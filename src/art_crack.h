#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <algorithm>
#include "art.h"

using namespace std;

int nLeaves, nInternals, nCap, locked;

namespace art_crack {

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

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

template <typename T>
class Bucket {
  int P;              // Pending insert if positive or pending delete if negative.
  int N;              // Number of data elements in this bucket pointed by D.
  int nextp;          // Pointer to the next chained bucket in leaf_bucket_allocator.
  int tailp;          // Pointer to the last chained bucket in leaf_bucket_allocator.
  T D[LEAF_BSIZE];  // Data values.

 public:

  int size() const { assert(is_valid()); return N; }
  T data(int i) const { assert(is_valid()); assert(i >= 0 && i < N); return D[i]; };
  int next() const { assert(is_valid()); return nextp; }
  int tail() const { assert(is_valid()); return tailp; }
  int slack() const { assert(is_valid()); return LEAF_BSIZE - size(); }
  bool is_full() const { assert(is_valid()); return slack() == 0; }
  bool is_valid() const { return N <= LEAF_BSIZE; }
  int last_data_is_at_least(T value) const { assert(is_valid()); return D[N - 1] >= value; }
  void set_data(int i, T v) { assert(is_valid()); assert(i >=0 && i < N); D[i] = v; }

  void init() {
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

  bool move_all_data_to(Bucket *that) {
    if (!size() || size() + that->size() > LEAF_BSIZE) return false;
    while (N-- > 0) that->append(D[--N]);
    return true;
  }

  void move_data_at_least(T value, Bucket *that) {
    assert(is_valid());
    for (int i = 0; i < N; i++) {
      if (D[i] >= value) {
        that->append(D[i]);
        D[i--] = D[--N];
      }
    }
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
    while (pos < N && D[pos] < value) {
      ART_DEBUG("leaf_lower_pos for value = %lld, D[%d/%d] = %lld\n", value, pos, N, D[pos]);
      pos++;
    }
    if (art_debug) {
      for (int i = 0; i < N; i++)
       ART_DEBUG("D[%d/%d] = %lld\n", i, N, D[i]);
    }
    ART_DEBUG("leaf_lower_pos for value = %lld, pos = %d/%d, v = %lld\n", value, pos, N, pos < N ? D[pos] : -1);
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

  // Only swaps as necessary.
  void fusion(Bucket *that, int *hi, int *lo, int &nhi, int &nlo) {
    assert(is_valid());
    T *Lp = D, *Rp = that->D;
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
class ArtCrack {
  Allocator<Bucket<T>> leaf_bucket_allocator;
  Node* tree = NULL;

  Bucket<T>* LEAF_BUCKET(int leafb) {
    // fprintf(stderr, "leafb = %d\n", leafb);
    return leaf_bucket_allocator.get((leafb) - 1);
  }

  // Returns positive integer ID.
  int new_leaf_bucket() {
    int leafb = leaf_bucket_allocator.alloc() + 1;
    LEAF_BUCKET(leafb)->init();
    return leafb;
  }

  void delete_leaf_bucket(int &leafb) {
    LEAF_BUCKET(leafb)->destroy();
    leaf_bucket_allocator.destroy(leafb - 1);
    leafb = 0;
  }

  void distribute_values(int leafb, T pivot, int chain[2]) {
    while (LEAF_BUCKET(leafb)->size()) {
      int i = LEAF_BUCKET(leafb)->last_data_is_at_least(pivot);
      // fprintf(stderr, "distribute_values %d\n", LEAF_BUCKET(leafb)->size());
      leaf_insert(chain[i], LEAF_BUCKET(leafb)->leaf_promote_last());
    }
  }

  void leaf_split_one_chain(int leafb, T &promotedValue, int &new_leafb) {
    new_leafb = LEAF_BUCKET(leafb)->detach_and_get_next();
    assert(!LEAF_BUCKET(new_leafb)->next());
    assert(!LEAF_BUCKET(new_leafb)->tail());
    if (LEAF_BUCKET(new_leafb)->move_all_data_to(LEAF_BUCKET(leafb))) {
      delete_leaf_bucket(new_leafb);
    } else {
      T D[LEAF_BSIZE * 2];
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
  T pick_random_pivot(int leafb) {
    Bucket<T> *B = LEAF_BUCKET(leafb);
    Bucket<T> *C = LEAF_BUCKET(B->tail());

    // fprintf(stderr, "Picking random 11 elements, sizes = %d + %d\n", B->size(), T->size());
    while (B->size() < 11) B->append(C->leaf_promote_last());

    T R[11]; // Randomly pick 11 elements from B.
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

  void leaf_split_long_chain(int leafb, T &pivot, int &new_leafb) {
    // fprintf(stderr, "split long chain = %d\n", LEAF_BUCKET(leafb)->size());
    pivot = pick_random_pivot(leafb);
    new_leafb = new_leaf_bucket();
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

  void leaf_split(int leafb, T &promotedValue, int &new_leafb) {
    assert(LEAF_BUCKET(leafb)->next() != 0);
    new_leafb = 0;
    if (!LEAF_BUCKET(LEAF_BUCKET(leafb)->next())->next()) {
      leaf_split_one_chain(leafb, promotedValue, new_leafb);
    } else {
      leaf_split_long_chain(leafb, promotedValue, new_leafb);
    }
  }

  void debug_data(int b) {
      LEAF_BUCKET(b)->debug_data();
  }

  int leaf_size(int leafb) {
    int ret = 0;
    while (leafb != 0) {
      ret += LEAF_BUCKET(leafb)->size();
      leafb = LEAF_BUCKET(leafb)->next();
    }
    return ret;
  }

  void add_chain(int head, int next) {
    Bucket<T> *B = LEAF_BUCKET(head);
    if (B->next() == 0) {
      B->set_next(next);
      B->set_tail(next);
    } else {
      LEAF_BUCKET(B->tail())->set_next(next);
      B->set_tail(next);
    }
  }

  void leaf_insert(int leafb, T value) {
    if (!LEAF_BUCKET(leafb)->is_full()) {
      LEAF_BUCKET(leafb)->append(value);
      return;
    }
    int tail = LEAF_BUCKET(leafb)->tail();
    assert(tail == 0 || LEAF_BUCKET(tail)->next() == 0);
    if (tail == 0 || LEAF_BUCKET(tail)->is_full()) {
      add_chain(leafb, tail = new_leaf_bucket());
    }
    LEAF_BUCKET(tail)->append(value);
  }

 public:

  ArtCrack() {
    int max_size = 100000000;
    leaf_bucket_allocator.init(max_size / LEAF_BSIZE * 2);
  }

  void alloc_sizes(int &ia_free, int &ia_size, int &la_free, int &la_size) {
    la_free = leaf_bucket_allocator.free_indices.size();
    la_size = leaf_bucket_allocator.N;
  }

  pair<bool, T> lower_bound_bucket(int &b, T value) {
    // fprintf(stderr, "split %d\n", b);
    T right_pivot = -1;
    while (LEAF_BUCKET(b)->next()) {
      T v1 = LEAF_BUCKET(b)->data(0);
      remove_root_bucket(v1, b);

      T pivot;
      int nb;
      leaf_split_long_chain(b, pivot, nb);
      T t = LEAF_BUCKET(nb)->data(0);
      LEAF_BUCKET(nb)->set_data(0, pivot);
      leaf_insert(nb, t);

      ART_DEBUG("pivot = %lld\n", pivot);

      assert(pivot == LEAF_BUCKET(nb)->data(0));
      insert_root(pivot, nb);

      T v2 = LEAF_BUCKET(b)->data(0);
      insert_root(v2, b);

      if (v1 != v2) {
      }


      if (value >= pivot) {
        ART_DEBUG("GO RIGHT\n");
        b = nb;
      } else {
        ART_DEBUG("GO LEFT\n");
        right_pivot = pivot;
      }
    }
    // fprintf(stderr, "splited %d\n", b);

    int pos = LEAF_BUCKET(b)->leaf_lower_pos(value);
    if (pos < LEAF_BUCKET(b)->size())
      return make_pair(true, LEAF_BUCKET(b)->data(pos));
    ART_DEBUG("right_pivot = %lld\n", right_pivot);
    return make_pair(false, right_pivot);
  }

  pair<bool, T> lower_bound(T value) {
    ART_DEBUG("\nquery = %lld\n", value);
    uint8_t key[8];
    loadKey(value, key);
    Node *leaf = ::lower_bound_prev(tree,key,8,0,8);
    if (leaf) {
      assert(isLeaf(leaf));
      int b = reinterpret_cast<uintptr_t>(leaf) >> 1;
      auto ret = lower_bound_bucket(b, value);
      if (ret.first) return ret;
      if (ret.second != -1) {
        ret.first = true;
        return ret;
      }
    }
    ART_DEBUG("second lower_bound\n");
    leaf = ::lower_bound(tree,key,8,0,8);
    if (!leaf) return make_pair(false, 0);
    assert(leaf);
    assert(isLeaf(leaf));
    int b = reinterpret_cast<uintptr_t>(leaf) >> 1;
    auto ret = lower_bound_bucket(b, value);
    // assert(ret.first);
    ret.first = true;
    return ret;
  }

  int get_root(T value) {
    uint8_t key[8];
    loadKey(value, key);
    Node *leaf = ::lower_bound_prev(tree,key,8,0,8);
    if (leaf) {
      assert(isLeaf(leaf));
      return reinterpret_cast<uintptr_t>(leaf) >> 1;
    }
    leaf = ::lower_bound(tree,key,8,0,8);
    assert(leaf);
    assert(isLeaf(leaf));
    return reinterpret_cast<uintptr_t>(leaf) >> 1;
  }

  void insert_root(T value, int bucket_number) {
    ART_DEBUG("insert root at %lld, b = %d\n", value, bucket_number);
    assert(value == LEAF_BUCKET(bucket_number)->data(0));
    uint8_t key[8];
    loadKey(value, key);
    ::insert(tree,&tree,key,0,bucket_number,8);

    // TODO: remove
    Node *leaf = ::lookup(tree, key, 8, 0, 8);
    assert(leaf);
    assert(isLeaf(leaf));
    assert(getLeafValue(leaf) == value);
  }

  void remove_root_bucket(T value, int bucket_number) {
    ART_DEBUG("remove root at %lld, b = %d\n", value, bucket_number);
    uint8_t key[8];
    loadKey(value, key);

    // TODO: remove
    Node *leaf = ::lookup(tree, key, 8, 0, 8);
    assert(leaf);
    assert(isLeaf(leaf));

    ::erase(tree,&tree,key,8,0,8);

    // TODO: remove
    assert(!::lookup(tree, key, 8, 0, 8));
  }

  T bucket_head_value(long long b) {
    // fprintf(stderr, "b = %lld\n", b);
    // fprintf(stderr, "data = %lld\n", LEAF_BUCKET(b)->data(0));
    return LEAF_BUCKET(b)->data(0);
  }

  void insert(T value) {
    // fprintf(stderr, "ins %d\n", value);
    int b = get_root(value);
    assert(b != 0);
    leaf_insert(b, value);
  }

  void batch_insert(T *arr, int N) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0;
    int root = -1;
    while (i + LEAF_BSIZE <= N) {
      int idx = new_leaf_bucket();
      Bucket<T> *b = LEAF_BUCKET(idx);
      for (int j = 0; j < LEAF_BSIZE; j++) {
        b->append(arr[i++]);
      }
      // fprintf(stderr, "chain %d, %d\n", i, b->size());
      if (root == -1) root = idx;
      else add_chain(root, idx);
    }

    insert_root(LEAF_BUCKET(root)->data(0), root);

    // fprintf(stderr, "done %d %d\n", i, N);
    while (i < N) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "inserted %d elements\n", size());
  }

  bool erase(T value) {
    // assert(check());
    // fprintf(stderr, "ERASE %d\n", value);
    // debug();

    uint8_t key[8];
    loadKey(value, key);
    Node *leaf = ::lower_bound_prev(tree,key,8,0,8);
    if (leaf) {
      assert(isLeaf(leaf));
      int b = reinterpret_cast<uintptr_t>(leaf) >> 1;
      auto ret = lower_bound_bucket(b, value);
      if (ret.first) {
        return LEAF_BUCKET(b)->leaf_erase(value);
      }
    }
    leaf = ::lower_bound(tree,key,8,0,8);
    if (!leaf) return false;
    assert(leaf);
    assert(isLeaf(leaf));
    int b = reinterpret_cast<uintptr_t>(leaf) >> 1;
    lower_bound_bucket(b, value);
    return LEAF_BUCKET(b)->leaf_erase(value);

    // assert(next_largest <= value);
    // assert(check());
  }
};


}


#endif
