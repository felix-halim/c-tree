#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <algorithm>

#include "random.h"

using namespace std;

int nLeaves, nInternals, nCap, locked;

namespace ctree {

#ifndef INTERNAL_BSIZE
  #define INTERNAL_BSIZE 64
#endif

#define LEAF_SIZE 64
#define LEAF_BUFFER_SIZE 4096

template <typename T>
class Allocator {
 public:

  priority_queue<int, vector<int>, greater<int>> free_indices;
  int cap;
  int N;
  T *D;

  void init(int initial_cap) {
    D = new T[cap = initial_cap];
    N = 0;
  }

  void save(string fn) {
    FILE *out = fopen(fn.c_str(), "wb");
    fwrite(&N, sizeof(int), 1, out);
    fwrite(D, sizeof(T), N, out);
    fclose(out);
  }

  void load(string fn) {
    FILE *f = fopen(fn.c_str(), "rb");
    fread(&N, sizeof(int), 1, f);
    D = new T[cap = N];
    fread(D, sizeof(T), N, f);
    fclose(f);
  }

  // TODO: mmap support
  int alloc() {
    if (free_indices.empty()) {
      if (N == cap) {
        fprintf(stderr, "double %d\n", cap);
        T *newD = new T[cap * 2];
        memcpy(newD, D, sizeof(T) * cap);
        cap *= 2;
        delete[] D;
        D = newD;
      }
      return N++;
    }
    int idx = free_indices.top();
    free_indices.pop();
    return idx;
  }

  void destroy(int idx) {
    free_indices.push(idx);
  }

  T* get(int idx) const {
    assert(idx >= 0);
    assert(idx < N);
    return &D[idx];
  }
};

template <typename T, int LEAF_BSIZE>
class Bucket {
  int P;              // Pending insert if positive or pending delete if negative.
  int N;              // Number of data elements in this bucket pointed by D.
  int nextp;          // Pointer to the next chained bucket in leaf_bucket_allocator.
  int tailp;          // Pointer to the last chained bucket in leaf_bucket_allocator.
  int parentp;        // Pointer to the parent bucket in leaf_bucket_allocator.
  T D[LEAF_BSIZE];  // Data values.

 public:

  int size() const { assert(is_valid()); return N; }
  T data(int i) const { assert(is_valid()); assert(i >= 0 && i < N); return D[i]; };
  T* data_arr() { assert(is_valid()); return D; }
  int next() const { assert(is_valid()); return nextp; }
  int tail() const { assert(is_valid()); return tailp; }
  int parent() const { assert(is_valid()); return parentp; }
  int slack() const { assert(is_valid()); return LEAF_BSIZE - size(); }
  bool is_full() const { assert(is_valid()); return slack() == 0; }
  bool is_valid() const { return N <= LEAF_BSIZE; }
  int last_data_is_at_least(T value) const { assert(is_valid()); return D[N - 1] >= value; }

  void set_parent(int parent) {
    assert(is_valid());
    parentp = parent;
  }

  void init(int parent) {
    parentp = parent;
    P = 0;
    N = 0;
    nextp = 0;
    tailp = 0;
    nCap += LEAF_BSIZE;
    nLeaves++;
  }

  void destroy() {
    invalidate();
    nCap -= LEAF_BSIZE;
    nLeaves--;
  }

  void bulk_insert(T *arr, int sz) {
    assert(is_valid());
    assert(!is_full());
    memcpy(D, arr, sizeof(T) * sz);
    N = sz;
    P = 1;
  }

  int copy_data_to(T *to) {
    assert(is_valid());
    for (int i = 0; i < N; i++) {
      to[i] = D[i];
    }
    return N;
  }

  void copy_data_from(T *from, int cnt) {
    assert(is_valid());
    for (int i = 0; i < cnt; i++) {
      D[i] = from[i];
    }
    N = cnt;
    P = 1;
  }

  T remove_random_data(Random &rng) {
    assert(is_valid());
    int j = rng.nextInt(N);
    T ret = D[j];
    D[j] = D[--N];
    return ret;
  }

  void swap_random_data_with(T &R, Random &rng) {
    assert(is_valid());
    assert(N > 0);
    swap(R, D[rng.nextInt(N)]);
  }

  void append(T value) {
    assert(is_valid());
    assert(!is_full());
    D[N++] = value;
    P = 1;
  }

  void leaf_sort() {
    assert(is_valid());
    if (P) {
      assert(!locked);
      sort(D, D + N);
      P = 0;
    }
  }

  T leaf_promote_last() {
    assert(is_valid());
    assert(N > 0);
    return D[--N];
  }

  int leaf_lower_pos(T value) {
    assert(is_valid());
    leaf_sort();
    if (LEAF_BSIZE >= 256) {
      return std::lower_bound(D, D + N, value) - D;
    }
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }

  T leaf_promote_first() {
    assert(is_valid());
    // TODO: optimize
    P = 1;
    int smallest_pos = 0;
    int pos = 1;
    while (pos < N) {
      if (D[pos] < D[smallest_pos]) smallest_pos = pos;
      pos++;
    }
    swap(D[smallest_pos], D[--N]);
    return D[N];
  }

  int detach_and_get_next() {
    assert(is_valid());
    int ret = nextp;
    nextp = tailp = 0;
    return ret;
  }

  void mark_hi(T P, int *hi, int &nhi) {
    for (int i = 0; i < N; i++) {
      hi[nhi] = i;
      nhi += D[i] >= P;
    }
  }

  void mark_lo(T P, int *lo, int &nlo) {
    for (int i = 0; i < N; i++) {
      lo[nlo] = i;
      nlo += D[i] < P;
    }
  }

  void set_next(int next) {
    assert(is_valid());
    nextp = next;
  }

  void set_tail(int tail) {
    assert(is_valid());
    tailp = tail;
  }

  T leaf_erase_pos(int pos) {
    assert(is_valid());
    assert(pos >= 0 && pos < N);
    swap(D[pos], D[--N]);
    P = 1;
    return D[N];
  }

  int leaf_largest_pos() {
    assert(nextp == 0);
    // fprintf(stderr, "leaf_erase_largest\n");
    int largest_pos = 0, pos = 1;
    while (pos < N) {
      if (D[pos] >= D[largest_pos])
        largest_pos = pos;
      pos++;
    }
    return largest_pos;
  }

  bool leaf_erase(T v) {
    for (int i = 0; i < N; i++) {
      if (D[i] == v) {
        D[i] = D[--N];
        P = 1;
        return true;
      }
    }
    return false;
  }

  void debug_data() {
    fprintf(stderr, "N = %d  %sLEAF  [ ", N, P ? "UNSORTED_" : "SORTED_");
    for (int i = 0; i < N; i++) fprintf(stderr, "%d ", D[i]);
    fprintf(stderr, "]");
  }

  bool leaf_check(T lo, bool useLo, T hi, bool useHi) {
    if (useLo) for (int i = 0; i < N; i++) if (D[i] < lo) {
      fprintf(stderr,"useLo failed: D[%d] = %lld, lo = %lld\n", i, D[i], lo);
      return false;
    }
    if (useHi) for (int i = 0; i < N; i++) if ((D[i] > hi)) {
      fprintf(stderr,"useHi failed: D[%d] = %lld, hi = %lld\n", i, D[i], hi);
      return false;
    }
    return true;
  }

  void invalidate() {
    assert(is_valid());
    N = LEAF_BSIZE + 1;
  }
};

template <typename T>
class IBucket {
  int N;                       // Number of data elements in this bucket pointed by D.
  int parentp;                 // Pointer to the parent bucket in leaf_bucket_allocator.
  T D[INTERNAL_BSIZE];       // Data values.
  int C[INTERNAL_BSIZE + 1];   // Pointer to children.
  unsigned long long is_buffer;

 public:

  int size() const { assert(is_valid()); return N; };
  int slack() const { assert(is_valid()); return INTERNAL_BSIZE - size(); }
  bool is_full() const { assert(is_valid()); return slack() == 0; }
  bool is_valid() const { return N <= INTERNAL_BSIZE; }
  int parent() const { assert(is_valid()); return parentp; }
  int child(int i) const { assert(is_valid()); assert(i >= 0 && i <= N); return C[i]; };
  T data(int i) const { assert(is_valid()); assert(i >= 0 && i < N); return D[i]; };

  void set_data(int i, T value) { assert(is_valid()); assert(i >= 0 && i < N); D[i] = value; };
  void set_parent(int parent) { assert(is_valid()); parentp = parent; }

  void init(int parent, int left_child) {
    parentp = parent;
    N = 0;
    C[0] = left_child;
    nCap += INTERNAL_BSIZE;
    nInternals++;
  }

  void destroy() {
    invalidate();
    nCap -= INTERNAL_BSIZE;
    nInternals--;
  }

  T internal_promote_first(int *C) {
    assert(is_valid());
    T ret = D[0];
    N--;
    for (int i = 0; i < N; i++) {
      D[i] = D[i + 1];
      C[i] = C[i + 1];
    }
    C[N] = C[N + 1];
    return ret;
  }

  T internal_promote_last() {
    assert(is_valid());
    return D[--N];
  }

  int internal_lower_pos(T value) const {
    assert(is_valid());
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }

  void internal_insert(T value, int nb, int left) {
    assert(is_valid());
    assert(!is_full());
    // fprintf(stderr, "ii %d\n", b);
    int i = N - 1;
    while (i >= 0 && D[i] > value) {
      D[i + 1] = D[i];
      C[i + 2] = C[i + 1];
      i--;
    }
    D[i + 1] = value;
    if (left == -1) {
      C[i + 2] = C[i + 1];
      C[i + 1] = nb;
    } else {
      C[i + 2] = nb;
    }
    N++;
  }

  void internal_erase(int pos, int stride) {
    assert(is_valid());
    N--;
    while (pos < N) {
      D[pos] = D[pos + 1];
      C[pos + stride] = C[pos + stride + 1];
      pos++;
    }
    if (!stride) C[pos] = C[pos + 1];
  }

  int mid_child() {
    assert(is_valid());
    return C[N / 2];
  }

  void move_half_to(IBucket *that) {
    assert(is_valid());
    int H = N / 2;
    for (int i = H, j = 0; i < N; i++) {
      that->D[j++] = D[i];
      that->C[j] = C[i + 1];
    }
    that->N = N - H;
    N = H;
  }

  void debug_data() {
    assert(is_valid());
    fprintf(stderr, "N = %d  [ ", N);
    for (int i = 0; i < N; i++) fprintf(stderr, "%d ", D[i]);
    fprintf(stderr, "]");
  }

  bool equal(int pos, T value) {
    assert(is_valid());
    return pos >= 0 && pos < N && D[pos] == value;
  }

  void invalidate() {
    assert(is_valid());
    N = INTERNAL_BSIZE + 1;
  }
};

template <typename T>
class CTree {
  Allocator<Bucket<T, LEAF_SIZE>> leaf_bucket_allocator;
  Allocator<Bucket<T, LEAF_BUFFER_SIZE>> leaf_buffer_allocator;
  Allocator<IBucket<T>> internal_bucket_allocator;
  int root;

  bool is_leaf(int b) {
    assert(b);
    return b > 0;
  }

  bool is_buffer(int b) {
    assert(b);
    return b >= 1000000000;
  }

  Bucket<T, LEAF_SIZE>* LEAF_BUCKET(int leafb) {
    assert(is_leaf(leafb) && !is_buffer(leafb));
    return leaf_bucket_allocator.get((leafb) - 1);
  }

  Bucket<T, LEAF_BUFFER_SIZE>* LEAF_BUFFER(int leafb) {
    assert(is_leaf(leafb) && is_buffer(leafb));
    return leaf_buffer_allocator.get((leafb) - 1000000000);
  }

  IBucket<T>* INTERNAL_BUCKET(int internalb) {
    assert(!is_leaf(internalb));
    return internal_bucket_allocator.get(-(internalb) - 1);
  }

  // Returns positive integer ID.
  int new_leaf_bucket(int parent) {
    int leafb = leaf_bucket_allocator.alloc() + 1;
    LEAF_BUCKET(leafb)->init(parent);
    return leafb;
  }

  // Returns positive integer ID.
  int new_buffer_bucket(int parent) {
    int leafb = leaf_buffer_allocator.alloc() + 1000000000;
    LEAF_BUFFER(leafb)->init(parent);
    return leafb;
  }

  // Returns negative integer ID.
  int new_internal_bucket(int parent, int left_child) {
    int internalb = -(internal_bucket_allocator.alloc() + 1);
    INTERNAL_BUCKET(internalb)->init(parent, left_child);
    set_parent(left_child, internalb);
    return internalb;
  }

  void delete_leaf(int &leafb) {
    if (is_buffer(leafb)) {
      LEAF_BUFFER(leafb)->destroy();
      leaf_buffer_allocator.destroy(leafb - 1000000000);
    } else {
      LEAF_BUCKET(leafb)->destroy();
      leaf_bucket_allocator.destroy(leafb - 1);
    }
    leafb = 0;
  }

  void delete_internal_bucket(int &internalb) {
    INTERNAL_BUCKET(internalb)->destroy();
    internal_bucket_allocator.destroy(-internalb - 1);
    internalb = 0;
  }

  void distribute_values(int leafb, T pivot, int chain[2]) {
    if (is_buffer(leafb)) {
      while (leaf_size(leafb)) {
        int i = LEAF_BUFFER(leafb)->last_data_is_at_least(pivot);
        leaf_insert(chain[i], LEAF_BUFFER(leafb)->leaf_promote_last());
      }
    } else {
      while (LEAF_BUCKET(leafb)->size()) {
        int i = LEAF_BUCKET(leafb)->last_data_is_at_least(pivot);
        leaf_insert(chain[i], LEAF_BUCKET(leafb)->leaf_promote_last());
      }
    }
  }

  // Returns the right boundary of the leftmost (not yet promoted) bucket.
  int rec_to_leaves(T *D, int lo, int hi, int parent, bool left_most) {
    assert(lo < hi);
    if (hi - lo <= LEAF_SIZE) {
      if (left_most) return hi;
      int b = new_leaf_bucket(parent);
      LEAF_BUCKET(b)->copy_data_from(D + lo, hi - lo);
      int promotedValue = LEAF_BUCKET(b)->leaf_promote_first();
      promote_internal(parent, promotedValue, b);
      return hi;
    }

    int mid = (lo + hi) / 2;
    nth_element(D + lo, D + mid, D + hi);

    int ret = rec_to_leaves(D, lo, mid, parent, left_most);
    rec_to_leaves(D, mid, hi, parent, false);
    return ret;
  }

  int leaf_parent(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->parent()
      : LEAF_BUCKET(leafb)->parent();
  }

  int copy_data_to(int leafb, T *D) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->copy_data_to(D)
      : LEAF_BUCKET(leafb)->copy_data_to(D);
  }

  int detach_and_get_next(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->detach_and_get_next()
      : LEAF_BUCKET(leafb)->detach_and_get_next();
  }

  int get_next(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->next()
      : LEAF_BUCKET(leafb)->next();
  }

  int get_tail(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->tail()
      : LEAF_BUCKET(leafb)->tail();
  }

  void leaf_split_one_chain(int &leafb) {
    assert(is_leaf(leafb));
    T D[LEAF_BUFFER_SIZE * 2];
    int N = copy_data_to(leafb, D);
    int parent = leaf_parent(leafb);
    int new_leafb = detach_and_get_next(leafb);
    assert(!get_next(new_leafb));
    assert(!get_tail(new_leafb));

    N += copy_data_to(new_leafb, D + N);
    assert(N <= LEAF_BUFFER_SIZE * 2);

    delete_leaf(leafb);
    delete_leaf(new_leafb);

    N = rec_to_leaves(D, 0, N, parent, true);

    leafb = new_leaf_bucket(parent);
    LEAF_BUCKET(leafb)->copy_data_from(D, N);
  }

  void append(int leafb, T value) {
    is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->append(value)
      : LEAF_BUCKET(leafb)->append(value);
  }

  int leaf_size(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->size()
      : LEAF_BUCKET(leafb)->size();
  }

  T remove_random_data(int leafb, Random &rng) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->remove_random_data(rng)
      : LEAF_BUCKET(leafb)->remove_random_data(rng);
  }

  void swap_random_data_with(int leafb, int &value, Random &rng) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->swap_random_data_with(value, rng)
      : LEAF_BUCKET(leafb)->swap_random_data_with(value, rng);
  }

  // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
  T pick_random_pivot(int leafb) {
    Bucket<T, LEAF_BUFFER_SIZE> *C = LEAF_BUFFER(get_tail(leafb));

    // fprintf(stderr, "Picking random 11 elements, sizes = %d + %d\n", B->size(), T->size());
    while (leaf_size(leafb) < 11) append(leafb, C->leaf_promote_last());

    T R[11]; // Randomly pick 11 elements from B.
    Random rng(140384); // TODO: use randomized seed.
    for (int i = 0; i < 11; i++) R[i] = remove_random_data(leafb, rng);

    // Replace R with the next buckets in the chain using reservoir sampling.
    for (int i = 1, Nb = get_next(leafb); Nb; i++) {
      int j = rng.nextInt(i);
      if (j < 11) swap_random_data_with(Nb, R[j], rng);
      Nb = get_next(Nb);
    }

    std::nth_element(R, R + 5, R + 11);
    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      if (i != 5) append(leafb, R[i]);
    }
    // fprintf(stderr, "pivot = %d\n", pivot);
    return R[5];
  }

  T data(int leafb, int i) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->data(i)
      : LEAF_BUCKET(leafb)->data(i);
  }

  T leaf_erase_pos(int leafb, int i) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->leaf_erase_pos(i)
      : LEAF_BUCKET(leafb)->leaf_erase_pos(i);
  }

  void move_data_at_least(int leafb, int &value, int to_leaf) {
    for (int i = 0; i < leaf_size(leafb); i++) {
      T d = data(leafb, i);
      if (d >= value) {
        append(to_leaf, d);
        leaf_erase_pos(leafb, i--);
      }
    }
  }

  int is_full(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->is_full()
      : LEAF_BUCKET(leafb)->is_full();
  }

  void mark_hi(int leafb, T &pivot, int *hi, int &nhi) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->mark_hi(pivot, hi, nhi)
      : LEAF_BUCKET(leafb)->mark_hi(pivot, hi, nhi);
  }

  void mark_lo(int leafb, T &pivot, int *lo, int &nlo) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->mark_lo(pivot, lo, nlo)
      : LEAF_BUCKET(leafb)->mark_lo(pivot, lo, nlo);
  }

  T* data_arr(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->data_arr()
      : LEAF_BUCKET(leafb)->data_arr();
  }

  // Only swaps as necessary.
  void fusion(T *Lp, T *Rp, int *hi, int *lo, int &nhi, int &nlo) {
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  void leaf_split_long_chain(int leafb) {
    // fprintf(stderr, "split long chain = %d\n", LEAF(leafb)->size());
    T pivot = pick_random_pivot(leafb);
    int new_leafb = new_leaf_bucket(leaf_parent(leafb));
    move_data_at_least(leafb, pivot, new_leafb);
    int chain[2] { leafb, new_leafb };

    // fprintf(stderr, "fusion %d\n", leafb);
    int Nb = detach_and_get_next(leafb);
    int Lb = 0, Rb = 0;
    // TODO: optimize locality.
    int hi[LEAF_BUFFER_SIZE], nhi = 0;
    int lo[LEAF_BUFFER_SIZE], nlo = 0;
    while (true) {
      if (nhi && nlo) {
        assert(Lb != 0 && Rb != 0);
        fusion(data_arr(Lb), data_arr(Rb), hi, lo, nhi, nlo);
        if (!nhi) { add_chain(chain[0], Lb); Lb = 0; }
        if (!nlo) { add_chain(chain[1], Rb); Rb = 0; }
      } else if (Lb == 0) {
        if (Nb == 0) break;
        Lb = Nb;
        Nb = detach_and_get_next(Nb);
        if (!is_full(Lb)) break;
      } else if (!nhi) {
        assert(Lb != 0);
        mark_hi(Lb, pivot, hi, nhi);
        if (!nhi){ add_chain(chain[0], Lb); Lb = 0; }
      } else if (Rb == 0) {
        if (Nb == 0) break;
        Rb = Nb;
        Nb = detach_and_get_next(Nb);
        if (!is_full(Rb)) break;
      } else if (!nlo) {
        assert(Rb != 0);
        mark_lo(Rb, pivot, lo, nlo);
        if (!nlo){ add_chain(chain[1], Rb); Rb = 0; }
      } else {
        assert(0);
      }
    }
    assert(Nb == 0);

    // fprintf(stderr, "fusioned\n");
    if (Lb != 0) distribute_values(Lb, pivot, chain), delete_leaf(Lb);
    if (Rb != 0) distribute_values(Rb, pivot, chain), delete_leaf(Rb);
    // fprintf(stderr, "done split long chain\n");
    // assert(leaf_check());

    promote_internal(leaf_parent(leafb), pivot, new_leafb);
  }

  void set_parent(int child, int parent) {
    if (is_leaf(child)) {
      is_buffer(child)
        ? LEAF_BUFFER(child)->set_parent(parent)
        : LEAF_BUCKET(child)->set_parent(parent);
    } else {
      INTERNAL_BUCKET(child)->set_parent(parent);
    }
  }

  void internal_insert(int internalb, T value, int newb, int left = 0) {
    INTERNAL_BUCKET(internalb)->internal_insert(value, newb, left);
    set_parent(newb, internalb);
    // assert(check());
  }

  void internal_erase(int internalb, int pos, int stride) {
    INTERNAL_BUCKET(internalb)->internal_erase(pos, stride);
  }

  int internal_split(int internalb) {
    int new_internalb = new_internal_bucket(
      INTERNAL_BUCKET(internalb)->parent(),
      INTERNAL_BUCKET(internalb)->mid_child());

    IBucket<T> *ib = INTERNAL_BUCKET(new_internalb);
    INTERNAL_BUCKET(internalb)->move_half_to(ib);
    for (int i = 0; i <= ib->size(); i++) {
      set_parent(ib->child(i), new_internalb);
    }
    return new_internalb;
  }

  void promote_internal(int parent, int &promotedValue, int &nb) {
    // fprintf(stderr, "parent = %d\n", parent);
    while (parent != 0 && nb != 0) {
      if (INTERNAL_BUCKET(parent)->is_full()) {
        // fprintf(stderr, "parful\n");

        // Optional optimization:
        // transfer_one_to_left_or_right();
        
        int inb = internal_split(parent);
        T promotedValueInternal = INTERNAL_BUCKET(parent)->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          internal_insert(inb, promotedValue, nb);
        } else {
          internal_insert(parent, promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = INTERNAL_BUCKET(nb)->parent();
      } else {
        // fprintf(stderr, "internal\n");
        internal_insert(parent, promotedValue, nb);
        nb = 0;
      }
    }
    assert(nb == 0 || parent == 0);

    // assert(check());
    // fprintf(stderr, "nb = %d\n", nb);
    if (nb != 0) {
      // fprintf(stderr, "OLD ROOT %d\n", root);
      assert(parent == 0);
      if (is_leaf(root)) {
        assert(is_leaf(root));
        assert(leaf_parent(root) == 0);
      } else {
        assert(INTERNAL_BUCKET(root)->parent() == 0);
      }
      root = new_internal_bucket(0, root);
      internal_insert(root, promotedValue, nb);
      // fprintf(stderr, "NEW ROOT %d\n", root);
    }
  }

  bool split_chain(int &leafb) {
    // assert(check());
    assert(is_leaf(leafb));
    if (get_next(leafb) == 0) return false;
    assert(!locked);

    // fprintf(stderr, "split_chain %d\n", leafb);
    // assert(check());

    if (!get_next(get_next(leafb))) {
      leaf_split_one_chain(leafb);
    } else {
      leaf_split_long_chain(leafb);
    }
    // fprintf(stderr, "leaf_split2\n");

    // fprintf(stderr, "split_chain2 %d\n", leafb);
    // assert(check());
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);


    // fprintf(stderr, "\n\n\ndone split %d\n", b);
    // debug();
    // assert(check());
    return true;
  }

  void debug_data(int b) {
    if (is_leaf(b)) {
      assert(is_leaf(b));
      // LEAF(b)->debug_data();
    } else {
      INTERNAL_BUCKET(b)->debug_data();
    }
  }

  int leaf_sizes(int leafb) {
    int ret = 0;
    while (leafb != 0) {
      ret += size(leafb);
      leafb = get_next(leafb);
    }
    return ret;
  }

  // Returns <bucket, pos> if found in internal node,
  // Otherwise returns <bucket, splitted> for leaf node.
  pair<int, int> find_bucket(T value, bool include_internal) {
    int *b = &root, splitted = 0;
    // fprintf(stderr, "find_bucket %d\n", b);
    while (true) {
      if (is_leaf(*b)) {
        // fprintf(stderr, "find_bucket2 %d\n", b);
        if (!split_chain(*b)) break;
        // fprintf(stderr, "find_bucket3 %d\n", b);
        assert(is_leaf(*b));
        *b = leaf_parent(*b);
        assert(*b != 0);
        splitted = 1;
      } else {
        int pos = INTERNAL_BUCKET(*b)->internal_lower_pos(value);
        if (include_internal && INTERNAL_BUCKET(*b)->equal(pos, value)) {
          // fprintf(stderr, "find_bucket2 %d\n", b);
          return make_pair(*b, pos); // Found in the internal bucket.
        }
        *b = INTERNAL_BUCKET(*b)->child(pos);    // Search the child.
      }
    }
    // fprintf(stderr, "find_bucket3 %d\n", b);
    return make_pair(*b, splitted);
  }

  void set_next(int leafb, int next) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->set_next(next)
      : LEAF_BUCKET(leafb)->set_next(next);
  }

  void set_tail(int leafb, int tail) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->set_tail(tail)
      : LEAF_BUCKET(leafb)->set_tail(tail);
  }

  void add_chain(int head, int next) {
    // assert(BUCKET(next)->is_leaf());
    assert(is_leaf(head));
    if (get_next(head) == 0) {
      set_next(head, next);
      set_tail(head, next);
    } else {
      int tail = get_tail(head);
      assert(is_leaf(tail));
      set_next(tail, next);
      set_tail(head, next);
    }
  }

  void leaf_insert(int leafb, T value) {
    assert(is_leaf(leafb));
    if (!is_full(leafb)) {
      append(leafb, value);
      return;
    }
    int tail = get_tail(leafb);
    assert(tail == 0 || get_next(tail) == 0);
    if (tail == 0 || is_full(tail)) {
      add_chain(leafb, tail = new_leaf_bucket(0));
    }
    append(tail, value);
  }

  bool check(int b, T lo, T hi) {
    assert(is_leaf(b));
    // return is_leaf(b) ? LEAF(b)->leaf_check(lo, true, hi, true) : check(b, lo);
  }

  int parent_of(int b) {
    return is_leaf(b) ? leaf_parent(b) : INTERNAL_BUCKET(b)->parent();
  }

 public:

  CTree() {
    int max_size = 100000000;
    leaf_bucket_allocator.init(max_size / LEAF_SIZE * 2);
    internal_bucket_allocator.init(max_size / LEAF_SIZE / INTERNAL_BSIZE * 3);
    root = new_leaf_bucket(0);
  }

  void alloc_sizes(int &ia_free, int &ia_size, int &la_free, int &la_size) {
    ia_free = internal_bucket_allocator.free_indices.size();
    ia_size = internal_bucket_allocator.N;
    la_free = leaf_bucket_allocator.free_indices.size();
    la_size = leaf_bucket_allocator.N;
  }

  void leaf_sort(int leafb) {
    is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->leaf_sort()
      : LEAF_BUCKET(leafb)->leaf_sort();
  }

  bool optimize(int &b) {
    if (b == 0) {
      if (root == 0) return false;
      while (optimize(root)) {
        if (!is_leaf(root) && INTERNAL_BUCKET(root)->size() == 0) {
          // fprintf(stderr, "PROMOTE ROOT\n");
          int c = INTERNAL_BUCKET(root)->child(0);
          delete_leaf(root);
          root = c;
          set_parent(root, 0);
        }
        // fprintf(stderr, "optimize\n");
      }
      return false;
    }

    if (is_leaf(b)) {
      bool splitted = split_chain(b);
      while (split_chain(b));
      assert(is_leaf(b));
      leaf_sort(b);
      return splitted;
    }

    bool changed = false;
    for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
      if (optimize(INTERNAL_BUCKET(b)->child(i))) {
        i = -1;
        changed = true;
      }
    }

    // OPTIONAL:
    // for (int i = 0; i < BUCKET(b)->size(); i++) {
    //   int L = child(b, i);
    //   int R = child(b, i + 1);
    //   assert(BUCKET(L)->is_leaf() == BUCKET(R)->is_leaf());
    //   if (BUCKET(L)->is_leaf()) {
    //     if (leaf_shift_left(b, i)) changed = 1;
    //   } else {
    //     if (internal_shift_left(b, i)) changed = 1;
    //   }
    // }

    // fprintf(stderr, "internal %d %d\n", b, changed);
    return changed;
  }

  int debug(int b = 0, int depth = 0) {
    for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
    if (b == 0) b = root;
    int sz = 0;
    if (is_leaf(b)) {
      while (true) {
        assert(is_leaf(b));
        sz += leaf_size(b);
        // LEAF(b)->debug_data();
        b = get_next(b);
        if (b == 0) break;
        fprintf(stderr, "------ ");
      }
      fprintf(stderr, "\n");
    } else {
      INTERNAL_BUCKET(b)->debug_data();
      fprintf(stderr, "\n");
      for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
        sz += debug(INTERNAL_BUCKET(b)->child(i), depth + 1);
      }
      for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
      fprintf(stderr, "size = %d\n", sz);
    }
    return sz;
  }

  int max_depth(int b = 0) {
    if (b == 0) b = root;
    if (is_leaf(b)) return 1;
    int ret = 0;
    for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
      int d = max_depth(INTERNAL_BUCKET(b)->child(i));
      assert(ret == 0 || ret == d);
      ret = d;
    }
    return ret + 1;
  }

  int slack(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->slack()
      : LEAF_BUCKET(leafb)->slack();
  }

  int slack() { return slack(0, 0); }
  int slack(int b, int last) {
    if (b == 0) b = root;
    int ret = is_leaf(b) ? slack(b) : INTERNAL_BUCKET(b)->slack();
    // if (ret > 10) fprintf(stderr, "slack = %d, for leaf = %d, last = %d\n", ret, b->is_leaf(), last);
    if (is_leaf(b)) return ret;
    for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
      ret += slack(INTERNAL_BUCKET(b)->child(i), i >= INTERNAL_BUCKET(b)->size());
    }
    return ret;
  }

  int size(int b = 0) {
    if (b == 0) b = root;
    if (is_leaf(b)) return leaf_size(b);
    int ret = 0;
    for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
      ret += size(INTERNAL_BUCKET(b)->child(i));
    }
    return ret;
  }

  int leaf_lower_pos(int leafb, int value) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->leaf_lower_pos(value)
      : LEAF_BUCKET(leafb)->leaf_lower_pos(value);
  }

  pair<bool, T> lower_bound(T value) {
    // TODO: optimize leaf slack

    // assert(check());
    // fprintf(stderr, "lower_bound %d\n", value);
    pair<int, int> p = find_bucket(value, true);
    // fprintf(stderr, "lower_bound1 %d\n", value);

    // Found in internal bucket.
    pair<bool, T> ret = make_pair(false, 0);
    if (!is_leaf(p.first)) {
      ret = make_pair(true, value);
    } else {
      assert(is_leaf(p.first));
      int pos = leaf_lower_pos(p.first, value);
      if (pos < leaf_size(p.first)) {
        ret = make_pair(true, data(p.first, pos));

        // OPTIONAL optimization:
        // int parent = BUCKET(p.first)->parent;
        // if (parent != 0) {
        //   leaf_compact(parent);
        // }
      } else {
        assert(is_leaf(p.first));
        int b = leaf_parent(p.first);
        while (b != 0) {
          IBucket<T> *ib = INTERNAL_BUCKET(b);
          pos = ib->internal_lower_pos(value);
          if (pos < ib->size()) {
            ret = make_pair(true, ib->data(pos));
            break;
          }
          b = ib->parent();
        }
      }
    }
    return ret;
  }

  void insert(T value) {
    // fprintf(stderr, "ins %d\n", value);
    int b = root;
    assert(b != 0);
    while (true) {
      if (is_leaf(b)) break;
      int pos = INTERNAL_BUCKET(b)->internal_lower_pos(value);
      b = INTERNAL_BUCKET(b)->child(pos);
    }
    leaf_insert(b, value);
  }

  void bulk_insert(int leafb, T *arr, int sz) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->bulk_insert(arr, sz)
      : LEAF_BUCKET(leafb)->bulk_insert(arr, sz);
  }

  void batch_insert(T *arr, int N) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0;
    while (i + LEAF_SIZE <= N) {
      int idx = new_leaf_bucket(0);
      assert(is_leaf(idx));
      bulk_insert(idx, arr + i, LEAF_SIZE);
      i += LEAF_SIZE;
      // fprintf(stderr, "chain %d, %d\n", i, b->size());
      add_chain(root, idx);
    }
    // fprintf(stderr, "done %d %d\n", i, N);
    while (i < N) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "inserted %d elements\n", size());
  }

  int leaf_erase(int leafb, T value) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->leaf_erase(value)
      : LEAF_BUCKET(leafb)->leaf_erase(value);
  }

  int leaf_largest_pos(int leafb) {
    return is_buffer(leafb)
      ? LEAF_BUFFER(leafb)->leaf_largest_pos()
      : LEAF_BUCKET(leafb)->leaf_largest_pos();
  }

  bool erase(T value) {
    // assert(check());
    // fprintf(stderr, "ERASE %d\n", value);
    // debug();

    pair<int, int> p = find_bucket(value, true);
    if (is_leaf(p.first)) {
      // fprintf(stderr, "ERASE1 %d\n", value);
      return leaf_erase(p.first, value);
    }

    // Found in an internal node, delete the largest node <= value.
    auto upper = find_bucket(value, false);

    // It is possible that finding upper invalidated p's references.
    if (upper.second) p = find_bucket(value, true); // Refresh.

    // fprintf(stderr, "ERASE2 %d\n", value);

    // upper.first is the leaf bucket containing the value.
    // upper.second signify whether a leaf_split happened when finding the bucket.
    int b = upper.first, parent;
    T next_largest = 0;
    assert(is_leaf(b));
    if (leaf_size(b)) {
      // fprintf(stderr, "ERASE3 %d\n", value);
      int pos = leaf_largest_pos(b);
      next_largest = leaf_erase_pos(b, pos);
    } else {
      // Bucket b is empty, search ancestors.
      while (true) {
        // fprintf(stderr, "ERASE4 %d\n", b);
        assert(b != p.first);
        if (is_leaf(b)) {
          parent = leaf_parent(b);
          delete_leaf(b);
        } else {
          parent = INTERNAL_BUCKET(b)->parent();
          delete_internal_bucket(b);
        }
        assert(parent);
        if (INTERNAL_BUCKET(parent)->size() > 0) {
          if (parent == p.first) {
            // Found in the same internal node as p.first, just delete it.
            INTERNAL_BUCKET(p.first)->internal_erase(p.second, 0);
            // fprintf(stderr, "yay\n");
            return true;
          }
          next_largest = INTERNAL_BUCKET(parent)->internal_promote_last();
          break;
        }
        b = parent;
      }
    }
    // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
    assert(INTERNAL_BUCKET(p.first)->data(p.second) == value);
    assert(next_largest <= value);
    INTERNAL_BUCKET(p.first)->set_data(p.second, next_largest);
    // assert(check());
    // fprintf(stderr, "ERASE5 %d\n", value);
    return true;
  }

  bool check(int b = 0, T lo = -2147483648) {
    if (b == 0) b = root;
    // fprintf(stderr, "check %d, leaf = %d\n", b, BUCKET(b)->is_leaf());
    if (is_leaf(b)) {
      if (b == root) assert(leaf_parent(b) == 0);
      return true; //LEAF(b)->leaf_check(lo, true, 0, false);
    }
    IBucket<T> *ib = INTERNAL_BUCKET(b);
    if (parent_of(ib->child(0)) != b) {
      fprintf(stderr, "parent mismatch %d != %d\n", parent_of(ib->child(0)), b);
      return false;
    }
    if (ib->size() && !check(ib->child(0), lo, ib->data(0))) return false;
    for (int i = 0; i < ib->size(); i++) {
      assert(i == 0 || ib->data(i - 1) <= ib->data(i));
      if (parent_of(ib->child(i + 1)) != b) {
        fprintf(stderr, "parent mismatch internal %d != %d\n", parent_of(ib->child(i + 1)), b);
        return false;
      }
      if (!check(ib->child(i + 1), ib->data(i), (i + 1 < ib->size()) ? ib->data(i + 1) : 2147483647)) return false;
    }
    return true;
  }

  void save(string fn) {
    leaf_bucket_allocator.save(fn + ".leaf");
    internal_bucket_allocator.save(fn + ".internal");
  }

  void load(string fn) {
    leaf_bucket_allocator.load(fn + ".leaf");
    internal_bucket_allocator.load(fn + ".internal");
  }
};


/*
class CTree {

  bool compact_leaves(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    while (!L->is_full() && M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->set_data(lpos, M->leaf_promote_first());
    }
    if (!L->is_full() && !M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->internal_erase(lpos, 1);
      delete_bucket(M);
      return true; // M is empty, compaction is done.
    }

    Bucket *R = (Bucket*) ib->child(rpos + 1);
    assert(!R->next_bucket());

    // Move from M to R as many as possible.
    while (M->size()) {
      assert(!R->is_full());
      R->leaf_insert(ib->data(rpos));
      ib->set_data(rpos, M->leaf_promote_last());
    }
    assert(!R->is_full());
    R->leaf_insert(ib->data(rpos));
    ib->internal_erase(rpos, 0);
    delete_bucket(M);
    return true;
  }

  bool compact_internals(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    while (!L->is_full() && M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->set_data(lpos, M->internal_promote_first());
    }
    if (!L->is_full() && !M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->internal_erase(lpos, 1);
      delete M;
      return true; // M is empty, compaction is done.
    }

    Bucket *R = (Bucket*) ib->child(rpos + 1);
    assert(!R->next_bucket());

    // Move from M to R as many as possible.
    while (M->size()) {
      assert(!R->is_full());
      R->internal_insert(ib->data(rpos), M->child(M->size()), -1);
      ib->set_data(rpos, M->internal_promote_last());
    }
    assert(!R->is_full());
    R->internal_insert(ib->data(rpos), M->child(M->size()), -1);
    ib->internal_erase(rpos, 0);
    delete M;
    return true;
  }

  int leaf_shift_left(int b, int pos) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(BUCKET(L)->is_leaf());
    assert(BUCKET(R)->is_leaf());
    if (BUCKET(L)->next != 0) return false;
    if (BUCKET(R)->next != 0) return false;

    // Move from M to L as many as possible.
    int changed = 0;
    while (!BUCKET(L)->is_full() && BUCKET(R)->size()) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      BUCKET(b)->D[pos] = BUCKET(R)->leaf_promote_first();
      changed = 1;
    }
    if (!BUCKET(L)->is_full() && !BUCKET(R)->size()) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      internal_erase(b, CHILDREN(b), pos, 1);
      delete_bucket(R);
      return 2; // R is empty, compaction is done.
    }
    return changed;
  }

  bool internal_shift_left(int b, int pos, int numMove = INTERNAL_BSIZE + 1) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(!BUCKET(L)->is_leaf());
    assert(!BUCKET(R)->is_leaf());
    assert(BUCKET(L)->next == 0);
    assert(BUCKET(R)->next == 0);

    // Move from R to L as many as possible.
    bool changed = false;
    while (!BUCKET(L)->is_full() && BUCKET(R)->size() && numMove-- > 0) {
      internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
      BUCKET(b)->D[pos] = BUCKET(R)->internal_promote_first(CHILDREN(R));
      changed = true;
    }
    // if (!BUCKET(L)->is_full() && !BUCKET(R)->size()) {
    //   internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
    //   internal_erase(b, CHILDREN(b), pos, 1);
    //   delete_bucket(R);
    //   return true; // R is empty, compaction is done.
    // }
    return changed;
  }

  bool internal_shift_right(int b, int pos, int numMove = INTERNAL_BSIZE + 1) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(!BUCKET(L)->is_leaf());
    assert(!BUCKET(R)->is_leaf());
    assert(BUCKET(L)->next == 0);
    assert(BUCKET(R)->next == 0);

    // Move from L to R as many as possible.
    bool changed = false;
    while (!BUCKET(R)->is_full() && BUCKET(L)->size() && numMove-- > 0) {
      internal_insert(R, BUCKET(b)->D[pos], child(L, BUCKET(L)->size()), -1);
      BUCKET(b)->D[pos] = BUCKET(L)->internal_promote_last();
      changed = true;
    }
    // if (!BUCKET(R)->is_full());
    // R->internal_insert(ib->data(rpos), L->child(L->size()), -1);
    // ib->internal_erase(rpos, 0);
    // delete L;
    return changed;
  }


  int internal_find_child_pos(int b, int c) {
    int *C = CHILDREN(b);
    for (int i = 0; i <= BUCKET(b)->size(); i++) {
      if (C[i] == c) return i;
    }
    assert(0);
    return 0;
  }


  bool leaf_compact(int b, int start, int end) {
    int last = 0;
    for (int i = start; i < end; i++) {
      last = leaf_shift_left(b, i);
    }
    assert(last == 2);
    fprintf(stderr, "saved 1 leaf %d\n", b);
    return true;
  }

  bool leaf_compact(int b) {
    int slack = 0, start = 1;
    Bucket *B = BUCKET(CHILDREN(b)[0]);
    if (B->next == 0) slack += B->slack(), start = 0;
    for (int i = 1; i <= BUCKET(b)->size(); i++) {
      B = BUCKET(CHILDREN(b)[i]);
      if (B->next != 0) {
        slack += B->slack();
        if (slack >= INTERNAL_BSIZE)
          return leaf_compact(b, start, i);
      } else {
        return false;
        // if (slack) fprintf(stderr, "gathered slack = %d, %d\n", slack, i - start);
        slack = 0;
        start = i + 1;
      }
    }
    fprintf(stderr, "gathered slack = %d\n", slack);
    return false;
  }

  transfer_one_to_left_or_right() {
    assert(!BUCKET(parent)->is_leaf());
    int pp = BUCKET(parent)->parent;
    if (pp != 0) {
      assert(!BUCKET(pp)->is_leaf());
      int pos = internal_find_child_pos(pp, parent);
      // assert(check());

      if (pos > 0 && BUCKET(parent)->D[0] < promotedValue && internal_shift_left(pp, pos - 1, 1)) {
        // fprintf(stderr, "shift left %d, p = %d, b = %d, root = %d, pos = %d / %d\n", pp, parent, b, root, pos, BUCKET(pp)->size());
        assert(!BUCKET(parent)->is_full());
        // assert(check());
        internal_insert(parent, promotedValue, nb);
        nb = 0;
        // assert(check());
        break;
      } else if (pos < BUCKET(pp)->size() && promotedValue < BUCKET(parent)->D[BUCKET(parent)->size() - 1] && internal_shift_right(pp, pos, 1)) {
        // assert(check());
        assert(!BUCKET(parent)->is_full());
        internal_insert(parent, promotedValue, nb);
        nb = 0;
        // assert(check());
        break;
      }
    }
  }
*/
}

#endif
