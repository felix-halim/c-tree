#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace chrono;

int nLeaves, nInternals, nCap, locked;

namespace ctree {

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

#ifndef INTERNAL_BSIZE
  #define INTERNAL_BSIZE 64
#endif

#ifndef LEAF_BSIZE
  #define LEAF_BSIZE 128
#endif

template<typename T>
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

class Bucket {
  int P;              // Pending insert if positive or pending delete if negative.
  int N;              // Number of data elements in this bucket pointed by D.
  int nextp;          // Pointer to the next chained bucket in leaf_bucket_allocator.
  int tailp;          // Pointer to the last chained bucket in leaf_bucket_allocator.
  int parentp;        // Pointer to the parent bucket in leaf_bucket_allocator.
  int D[LEAF_BSIZE];  // Data values.

 public:

  int size() const { assert(is_valid()); return N; }
  int data(int i) const { assert(is_valid()); assert(i >= 0 && i < N); return D[i]; };
  int next() const { assert(is_valid()); return nextp; }
  int tail() const { assert(is_valid()); return tailp; }
  int parent() const { assert(is_valid()); return parentp; }
  int slack() const { assert(is_valid()); return LEAF_BSIZE - size(); }
  bool is_full() const { assert(is_valid()); return slack() == 0; }
  bool is_valid() const { return N <= LEAF_BSIZE; }
  int last_data_is_at_least(int value) const { assert(is_valid()); return D[N - 1] >= value; }

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

  int copy_data_to(int *to) {
    assert(is_valid());
    for (int i = 0; i < N; i++) {
      to[i] = D[i];
    }
    return N;
  }

  void copy_data_from(int *from, int cnt) {
    assert(is_valid());
    for (int i = 0; i < cnt; i++) {
      D[i] = from[i];
    }
    N = cnt;
    P = 1;
  }

  bool move_all_data_to(Bucket *that) {
    if (!size() || size() + that->size() > LEAF_BSIZE) return false;
    while (N-- > 0) that->append(D[--N]);
    return true;
  }

  void move_data_at_least(int value, Bucket *that) {
    assert(is_valid());
    for (int i = 0; i < N; i++) {
      if (D[i] >= value) {
        that->append(D[i]);
        D[i--] = D[--N];
      }
    }
  }

  int remove_random_data(Random &rng) {
    assert(is_valid());
    int j = rng.nextInt(N);
    int ret = D[j];
    D[j] = D[--N];
    return ret;
  }

  void swap_random_data_with(int &R, Random &rng) {
    assert(is_valid());
    assert(N > 0);
    swap(R, D[rng.nextInt(N)]);
  }

  void append(int value) {
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

  int leaf_promote_last() {
    assert(is_valid());
    assert(N > 0);
    return D[--N];
  }

  int leaf_lower_pos(int value) {
    assert(is_valid());
    leaf_sort();
    if (LEAF_BSIZE >= 256) {
      return std::lower_bound(D, D + N, value) - D;
    }
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }

  int leaf_promote_first() {
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

  void mark_hi(int P, int *hi, int &nhi) {
    for (int i = 0; i < N; i++) {
      hi[nhi] = i;
      nhi += D[i] >= P;
    }
  }

  void mark_lo(int P, int *lo, int &nlo) {
    for (int i = 0; i < N; i++) {
      lo[nlo] = i;
      nlo += D[i] < P;
    }
  }

  // Only swaps as necessary.
  void fusion(Bucket *that, int *hi, int *lo, int &nhi, int &nlo) {
    assert(is_valid());
    int *Lp = D, *Rp = that->D;
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  void set_next(int next) {
    assert(is_valid());
    nextp = next;
  }

  void set_tail(int tail) {
    assert(is_valid());
    tailp = tail;
  }

  int leaf_erase_pos(int pos) {
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

  bool leaf_erase(int v) {
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

  bool leaf_check(int lo, bool useLo, int hi, bool useHi) {
    if (useLo) for (int i = 0; i < N; i++) if (D[i] < lo) {
      fprintf(stderr,"useLo failed: D[%d] = %d, lo = %d\n", i, D[i], lo);
      return false;
    }
    if (useHi) for (int i = 0; i < N; i++) if ((D[i] > hi)) {
      fprintf(stderr,"useHi failed: D[%d] = %d, hi = %d\n", i, D[i], hi);
      return false;
    }
    return true;
  }

  void invalidate() {
    assert(is_valid());
    N = LEAF_BSIZE + 1;
  }
};

class IBucket {
  int N;                       // Number of data elements in this bucket pointed by D.
  int parentp;                 // Pointer to the parent bucket in leaf_bucket_allocator.
  int D[INTERNAL_BSIZE];       // Data values.
  int C[INTERNAL_BSIZE + 1];   // Pointer to children.

 public:

  int size() const { assert(is_valid()); return N; };
  int slack() const { assert(is_valid()); return INTERNAL_BSIZE - size(); }
  bool is_full() const { assert(is_valid()); return slack() == 0; }
  bool is_valid() const { return N <= INTERNAL_BSIZE; }
  int parent() const { assert(is_valid()); return parentp; }
  int child(int i) const { assert(is_valid()); assert(i >= 0 && i <= N); return C[i]; };
  int data(int i) const { assert(is_valid()); assert(i >= 0 && i < N); return D[i]; };

  void set_data(int i, int value) { assert(is_valid()); assert(i >= 0 && i < N); D[i] = value; };
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

  int internal_promote_first(int *C) {
    assert(is_valid());
    int ret = D[0];
    N--;
    for (int i = 0; i < N; i++) {
      D[i] = D[i + 1];
      C[i] = C[i + 1];
    }
    C[N] = C[N + 1];
    return ret;
  }

  int internal_promote_last() {
    assert(is_valid());
    return D[--N];
  }

  int internal_lower_pos(int value) const {
    assert(is_valid());
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }

  void internal_insert(int value, int nb, int left) {
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

  bool equal(int pos, int value) {
    assert(is_valid());
    return pos >= 0 && pos < N && D[pos] == value;
  }

  void invalidate() {
    assert(is_valid());
    N = INTERNAL_BSIZE + 1;
  }
};


class CTree {
  Allocator<Bucket> leaf_bucket_allocator;
  Allocator<IBucket> internal_bucket_allocator;
  int root;

  bool is_leaf(int b) {
    assert(b);
    return b > 0;
  }

  Bucket* LEAF_BUCKET(int leafb) {
    assert(is_leaf(leafb));
    return leaf_bucket_allocator.get((leafb) - 1);
  }

  IBucket* INTERNAL_BUCKET(int internalb) {
    assert(!is_leaf(internalb));
    return internal_bucket_allocator.get(-(internalb) - 1);
  }

  // Returns positive integer ID.
  int new_leaf_bucket(int parent) {
    int leafb = leaf_bucket_allocator.alloc() + 1;
    LEAF_BUCKET(leafb)->init(parent);
    return leafb;
  }

  // Returns negative integer ID.
  int new_internal_bucket(int parent, int left_child) {
    int internalb = -(internal_bucket_allocator.alloc() + 1);
    INTERNAL_BUCKET(internalb)->init(parent, left_child);
    set_parent(left_child, internalb);
    return internalb;
  }

  void delete_leaf_bucket(int &leafb) {
    LEAF_BUCKET(leafb)->destroy();
    leaf_bucket_allocator.destroy(leafb - 1);
    leafb = 0;
  }

  void delete_internal_bucket(int &internalb) {
    INTERNAL_BUCKET(internalb)->destroy();
    internal_bucket_allocator.destroy(-internalb - 1);
    internalb = 0;
  }

  void distribute_values(int leafb, int pivot, int chain[2]) {
    while (LEAF_BUCKET(leafb)->size()) {
      int i = LEAF_BUCKET(leafb)->last_data_is_at_least(pivot);
      // fprintf(stderr, "distribute_values %d\n", LEAF_BUCKET(leafb)->size());
      leaf_insert(chain[i], LEAF_BUCKET(leafb)->leaf_promote_last());
    }
  }

  void leaf_split_one_chain(int leafb, int &promotedValue, int &new_leafb) {
    new_leafb = LEAF_BUCKET(leafb)->detach_and_get_next();
    assert(!LEAF_BUCKET(new_leafb)->next());
    assert(!LEAF_BUCKET(new_leafb)->tail());
    if (LEAF_BUCKET(new_leafb)->move_all_data_to(LEAF_BUCKET(leafb))) {
      delete_leaf_bucket(new_leafb);
    } else {
      int D[LEAF_BSIZE * 2];
      int N = LEAF_BUCKET(leafb)->copy_data_to(D);
      N += LEAF_BUCKET(new_leafb)->copy_data_to(D + N);
      assert(N <= LEAF_BSIZE * 2);
      int H = N / 2;
      nth_element(D, D + H, D + N);
      LEAF_BUCKET(leafb)->copy_data_from(D, H);
      promotedValue = D[H++];
      LEAF_BUCKET(new_leafb)->copy_data_from(D + H, N - H);
    }
  }

  // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
  int pick_random_pivot(int leafb) {
    Bucket *B = LEAF_BUCKET(leafb);
    Bucket *T = LEAF_BUCKET(B->tail());

    // fprintf(stderr, "Picking random 11 elements, sizes = %d + %d\n", B->size(), T->size());
    while (B->size() < 11) B->append(T->leaf_promote_last());

    int R[11]; // Randomly pick 11 elements from B.
    Random rng(140384); // TODO: use randomized seed.
    for (int i = 0; i < 11; i++) R[i] = B->remove_random_data(rng);

    // Replace R with the next buckets in the chain using reservoir sampling.
    for (int i = 1, Nb = B->next(); Nb && i < 11; i++) {
      int j = rng.nextInt(i);
      if (j < 11) LEAF_BUCKET(Nb)->swap_random_data_with(R[j], rng);
      Nb = LEAF_BUCKET(Nb)->next();
    }

    std::nth_element(R, R + 5, R + 11);
    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      if (i != 5) B->append(R[i]);
    }
    // fprintf(stderr, "pivot = %d\n", pivot);
    return R[5];
  }

  void leaf_split_long_chain(int leafb, int &pivot, int &new_leafb) {
    // fprintf(stderr, "split long chain = %d\n", LEAF_BUCKET(leafb)->size());
    pivot = pick_random_pivot(leafb);
    new_leafb = new_leaf_bucket(LEAF_BUCKET(leafb)->parent());
    LEAF_BUCKET(leafb)->move_data_at_least(pivot, LEAF_BUCKET(new_leafb));
    int chain[2] { leafb, new_leafb };

    // fprintf(stderr, "fusion\n");
    int Nb = LEAF_BUCKET(leafb)->detach_and_get_next();
    int Lb = 0, Rb = 0;
    // TODO: optimize locality.
    int hi[LEAF_BSIZE], nhi = 0;
    int lo[LEAF_BSIZE], nlo = 0;
    while (true) {
      if (nhi && nlo) {
        assert(Lb != 0 && Rb != 0);
        LEAF_BUCKET(Lb)->fusion(LEAF_BUCKET(Rb), hi, lo, nhi, nlo);
        if (!nhi) { add_chain(chain[0], Lb); Lb = 0; }
        if (!nlo) { add_chain(chain[1], Rb); Rb = 0; }
      } else if (Lb == 0) {
        if (Nb == 0) break;
        Lb = Nb;
        Nb = LEAF_BUCKET(Nb)->detach_and_get_next();
        if (!LEAF_BUCKET(Lb)->is_full()) break;
      } else if (!nhi) {
        assert(Lb != 0);
        LEAF_BUCKET(Lb)->mark_hi(pivot, hi, nhi);
        if (!nhi){ add_chain(chain[0], Lb); Lb = 0; }
      } else if (Rb == 0) {
        if (Nb == 0) break;
        Rb = Nb;
        Nb = LEAF_BUCKET(Nb)->detach_and_get_next();
        if (!LEAF_BUCKET(Rb)->is_full()) break;
      } else if (!nlo) {
        assert(Rb != 0);
        LEAF_BUCKET(Rb)->mark_lo(pivot, lo, nlo);
        if (!nlo){ add_chain(chain[1], Rb); Rb = 0; }
      } else {
        assert(0);
      }
    }
    assert(Nb == 0);

    // fprintf(stderr, "fusioned\n");
    if (Lb != 0) distribute_values(Lb, pivot, chain), delete_leaf_bucket(Lb);
    if (Rb != 0) distribute_values(Rb, pivot, chain), delete_leaf_bucket(Rb);
    // fprintf(stderr, "done split long chain\n");
    // assert(leaf_check());
  }

  void set_parent(int child, int parent) {
    if (is_leaf(child)) {
      LEAF_BUCKET(child)->set_parent(parent);
    } else {
      INTERNAL_BUCKET(child)->set_parent(parent);
    }
  }

  void leaf_split(int leafb, int &promotedValue, int &new_leafb) {
    assert(LEAF_BUCKET(leafb)->next() != 0);
    new_leafb = 0;
    if (!LEAF_BUCKET(LEAF_BUCKET(leafb)->next())->next()) {
      leaf_split_one_chain(leafb, promotedValue, new_leafb);
    } else {
      leaf_split_long_chain(leafb, promotedValue, new_leafb);
    }
  }

  void internal_insert(int internalb, int value, int newb, int left = 0) {
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

    IBucket *ib = INTERNAL_BUCKET(new_internalb);
    INTERNAL_BUCKET(internalb)->move_half_to(ib);
    for (int i = 0; i <= ib->size(); i++) {
      set_parent(ib->child(i), new_internalb);
    }
    return new_internalb;
  }

  bool split_chain(int leafb) {
    // assert(check());
    if (LEAF_BUCKET(leafb)->next() == 0) return false;
    assert(!locked);

    int promotedValue;
    int nb;
    // fprintf(stderr, "split_chain %d, %d\n", BUCKET(b)->size(), BUCKET(b)->next);
    // assert(check());
    leaf_split(leafb, promotedValue, nb);
    // assert(check());
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);

    int parent = LEAF_BUCKET(leafb)->parent();
    // fprintf(stderr, "parent = %d\n", parent);
    while (parent != 0 && nb != 0) {
      if (INTERNAL_BUCKET(parent)->is_full()) {
        // fprintf(stderr, "parful\n");

        // Optional optimization:
        // transfer_one_to_left_or_right();
        
        int inb = internal_split(parent);
        int promotedValueInternal = INTERNAL_BUCKET(parent)->internal_promote_last();
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
    // assert(check());
    // fprintf(stderr, "nb = %d\n", nb);
    if (nb != 0) {
      // fprintf(stderr, "OLD ROOT %d\n", root);
      assert(parent == 0);
      if (is_leaf(root)) {
        assert(LEAF_BUCKET(root)->parent() == 0);
      } else {
        assert(INTERNAL_BUCKET(root)->parent() == 0);
      }
      root = new_internal_bucket(0, root);
      internal_insert(root, promotedValue, nb);
      // fprintf(stderr, "NEW ROOT %d\n", root);
    }
    // fprintf(stderr, "\n\n\ndone split %d\n", b);
    // debug();
    // assert(check());
    return true;
  }

  void debug_data(int b) {
    if (is_leaf(b)) {
      LEAF_BUCKET(b)->debug_data();
    } else {
      INTERNAL_BUCKET(b)->debug_data();
    }
  }

  int leaf_size(int leafb) {
    int ret = 0;
    while (leafb != 0) {
      ret += LEAF_BUCKET(leafb)->size();
      leafb = LEAF_BUCKET(leafb)->next();
    }
    return ret;
  }

  // Returns <bucket, pos> if found in internal node,
  // Otherwise returns <bucket, splitted> for leaf node.
  pair<int, int> find_bucket(int value, bool include_internal) {
    int b = root, splitted = 0;
    // fprintf(stderr, "find_bucket %d\n", b);
    while (true) {
      if (is_leaf(b)) {
        if (!split_chain(b)) break;
        b = LEAF_BUCKET(b)->parent();
        assert(b != 0);
        splitted = 1;
      } else {
        int pos = INTERNAL_BUCKET(b)->internal_lower_pos(value);
        if (include_internal && INTERNAL_BUCKET(b)->equal(pos, value)) {
          return make_pair(b, pos); // Found in the internal bucket.
        }
        b = INTERNAL_BUCKET(b)->child(pos);    // Search the child.
      }
    }
    return make_pair(b, splitted);
  }

  void add_chain(int head, int next) {
    // assert(BUCKET(next)->is_leaf());
    assert(is_leaf(head));
    Bucket *B = LEAF_BUCKET(head);
    if (B->next() == 0) {
      B->set_next(next);
      B->set_tail(next);
    } else {
      LEAF_BUCKET(B->tail())->set_next(next);
      B->set_tail(next);
    }
  }

  void leaf_insert(int leafb, int value) {
    if (!LEAF_BUCKET(leafb)->is_full()) {
      LEAF_BUCKET(leafb)->append(value);
      return;
    }
    int tail = LEAF_BUCKET(leafb)->tail();
    assert(tail == 0 || LEAF_BUCKET(tail)->next() == 0);
    if (tail == 0 || LEAF_BUCKET(tail)->is_full()) {
      add_chain(leafb, tail = new_leaf_bucket(0));
    }
    LEAF_BUCKET(tail)->append(value);
  }

  bool check(int b, int lo, int hi) {
    return is_leaf(b) ? LEAF_BUCKET(b)->leaf_check(lo, true, hi, true) : check(b, lo);
  }

  int parent_of(int b) {
    return is_leaf(b) ? LEAF_BUCKET(b)->parent() : INTERNAL_BUCKET(b)->parent();
  }

 public:

  CTree() {
    int max_size = 100000000;
    leaf_bucket_allocator.init(max_size / LEAF_BSIZE * 2);
    internal_bucket_allocator.init(max_size / LEAF_BSIZE / INTERNAL_BSIZE * 3);
    root = new_leaf_bucket(0);
  }

  void alloc_sizes(int &ia_free, int &ia_size, int &la_free, int &la_size) {
    ia_free = internal_bucket_allocator.free_indices.size();
    ia_size = internal_bucket_allocator.N;
    la_free = leaf_bucket_allocator.free_indices.size();
    la_size = leaf_bucket_allocator.N;
  }

  bool optimize(int b = 0) {
    if (b == 0) {
      if (root == 0) return false;
      while (optimize(root)) {
        if (!is_leaf(root) && INTERNAL_BUCKET(root)->size() == 0) {
          // fprintf(stderr, "PROMOTE ROOT\n");
          int c = INTERNAL_BUCKET(root)->child(0);
          delete_leaf_bucket(root);
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
      LEAF_BUCKET(b)->leaf_sort();
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
        sz += LEAF_BUCKET(b)->size();
        LEAF_BUCKET(b)->debug_data();
        b = LEAF_BUCKET(b)->next();
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

  int slack(int b = 0, int last = 0) {
    if (b == 0) b = root;
    int ret = is_leaf(b) ? LEAF_BUCKET(b)->slack() : INTERNAL_BUCKET(b)->slack();
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

  pair<bool, int> lower_bound(int value) {
    // TODO: optimize leaf slack

    // assert(check());
    // fprintf(stderr, "lower_bound %d\n", value);
    pair<int, int> p = find_bucket(value, true);
    // fprintf(stderr, "lower_bound1 %d\n", value);

    // Found in internal bucket.
    pair<bool, int> ret = make_pair(false, 0);
    if (!is_leaf(p.first)) {
      ret = make_pair(true, value);
    } else {
      int pos = LEAF_BUCKET(p.first)->leaf_lower_pos(value);
      if (pos < LEAF_BUCKET(p.first)->size()) {
        ret = make_pair(true, LEAF_BUCKET(p.first)->data(pos));

        // OPTIONAL optimization:
        // int parent = BUCKET(p.first)->parent;
        // if (parent != 0) {
        //   leaf_compact(parent);
        // }
      } else {
        int b = LEAF_BUCKET(p.first)->parent();
        while (b != 0) {
          IBucket *ib = INTERNAL_BUCKET(b);
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

  void insert(int value) {
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

  void batch_insert(int *arr, int N) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0;
    while (i + LEAF_BSIZE <= N) {
      int idx = new_leaf_bucket(0);
      Bucket *b = LEAF_BUCKET(idx);
      for (int j = 0; j < LEAF_BSIZE; j++) {
        b->append(arr[i++]);
      }
      // fprintf(stderr, "chain %d, %d\n", i, b->size());
      add_chain(root, idx);
    }
    // fprintf(stderr, "done %d %d\n", i, N);
    while (i < N) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "inserted %d elements\n", size());
  }

  bool erase(int value) {
    // assert(check());
    // fprintf(stderr, "ERASE %d\n", value);
    // debug();

    pair<int, int> p = find_bucket(value, true);
    if (is_leaf(p.first)) return LEAF_BUCKET(p.first)->leaf_erase(value);

    // Found in an internal node, delete the largest node <= value.
    auto upper = find_bucket(value, false);

    // It is possible that finding upper invalidated p's references.
    if (upper.second) p = find_bucket(value, true); // Refresh.

    // upper.first is the leaf bucket containing the value.
    // upper.second signify whether a leaf_split happened when finding the bucket.
    int b = upper.first;
    int next_largest = 0;
    assert(is_leaf(b));
    if (LEAF_BUCKET(b)->size()) {
      int pos = LEAF_BUCKET(b)->leaf_largest_pos();
      next_largest = LEAF_BUCKET(b)->leaf_erase_pos(pos);
    } else {
      // Bucket b is empty, search ancestors.
      while (true) {
        assert(b != p.first);
        assert(LEAF_BUCKET(b)->size() == 0);
        int parent = LEAF_BUCKET(b)->parent();
        assert(parent);
        if (is_leaf(b)) {
          delete_leaf_bucket(b);
        } else {
          delete_internal_bucket(b);
        }
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
    return true;
  }

  bool check(int b = 0, int lo = -2147483648) {
    if (b == 0) b = root;
    // fprintf(stderr, "check %d, leaf = %d\n", b, BUCKET(b)->is_leaf());
    if (is_leaf(b)) {
      if (b == root) assert(LEAF_BUCKET(b)->parent() == 0);
      return LEAF_BUCKET(b)->leaf_check(lo, true, 0, false);
    }
    IBucket *ib = INTERNAL_BUCKET(b);
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
