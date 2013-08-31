/*
use leaf chained as separate allocator
*/

#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <chrono>
#include <algorithm>
#include "random.h"

using namespace std;
using namespace chrono;

int nLeaves, nInternals, nCap, nDes, locked;

namespace ctree {

#define INTERNAL_BSIZE        64  // Must be power of two.
#define LEAF_BSIZE            64  // Must be power of two.
#define LEAF_CHAINED_BSIZE  2048  // Must be power of two.

#define BUCKET(b) bucket_allocator.get(b)
#define CHILDREN(b) (*child_allocator.get(BUCKET(b)->C))

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

template<typename T>
class Allocator {
  priority_queue<int, vector<int>, greater<int>> free_indices;
  int N, cap;
  T *D;

 public:

  Allocator() {
    D = new T[cap = 100000000/INTERNAL_BSIZE*4];
    N = 0;
    for (int i = 0; i < 1000000; i++) {
      free_indices.push(N++);
    }
  }

  int alloc(bool priority) {
    if (priority) {
      assert(!free_indices.empty());
      int idx = free_indices.top();
      free_indices.pop();
      return idx;
    }
    if (N == cap) {
      fprintf(stderr, "double %d\n", cap);
      T *newD = new T[cap * 4];
      memcpy(newD, D, sizeof(T) * cap);
      cap *= 4;
      delete[] D;
      D = newD;
    }
    return N++;
  }

  void destroy(int idx) {
    free_indices.push(idx);
  }

  T* get(int idx) const {
    assert(idx != -1);
    return &D[idx];
  }
};


struct Bucket {
  int C;       // Pointer to children.
  int P;       // Pending insert if positive or pending delete if negative.
  int N;       // Number of data elements in this bucket pointed by D.
  int cap;     // The capacity of this bucket pointed by D.
  int parent;  // Pointer to the parent bucket in bucket_allocator.
  int next;    // Pointer to the next chained bucket in bucket_allocator.
  int tail;    // Pointer to the last chained bucket in bucket_allocator.
  int D[LEAF_BSIZE];       // Data values.

  void init(int parent, int cap, int C = -1) {
    this->C = C;
    this->cap = cap;
    this->parent = parent;
    P = 0;
    N = 0;
    next = -1;
    tail = -1;

    nCap += cap;
    if (is_leaf()) {
      nLeaves++;
    } else {
      nInternals++;
    }
    // C[0] = left_child;
    // left_child->set_parent(this);
  }

  bool is_full() const { return N == cap; }
  bool is_leaf() const { return C == -1; }

  bool append(int value) {
    assert(is_leaf());
    if (is_full()) return false;
    D[N++] = value;
    P++;
    return true;
  }

  void leaf_optimize() {
    assert(P >= 0);
    sort(D, D + N);
    P = 0;
  }

  int leaf_lower_pos(int value) {
    if (P) leaf_optimize();
    // assert(leaf_check());
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    // assert(leaf_check());
    // leaf_debug();
    // fprintf(stderr, "leaf_lower_pos = %d, pos = %d, %d %d\n", value, pos, L, R);
    return pos;
  }

  int internal_lower_pos(int value) {
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }


  int leaf_promote_first() {
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

  int internal_promote_first(int *C) {
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
   return D[--N];
  }

  int detach_and_get_next() {
    int ret = next;
    next = tail = -1;
    return ret;
  }

  void destroy() {
    assert(cap > 0);
    nCap -= cap;
    if (is_leaf()) {
      nLeaves--;
    } else {
      nInternals--;
    }
    nDes++;
    cap = 0;
  }
};


class CTree {
  Allocator<Bucket> bucket_allocator;
  Allocator<int[INTERNAL_BSIZE + 1]> child_allocator;
  int root;

 public:

  const char *version = "Exp LEAF_BSIZE 2048";

  CTree() {
    root = -1;
  }

  int child(int b, int i) {
    assert(BUCKET(b)->C != -1);
    return CHILDREN(b)[i];
  }


  void mark_hi(int *D, int N, int P, int *hi, int &nhi) {
    for (int i = 0; i < N; i++) {
      hi[nhi] = i;
      nhi += D[i] >= P;
    }
  }

  void mark_lo(int *D, int N, int P, int *lo, int &nlo) {
    for (int i = 0; i < N; i++) {
      lo[nlo] = i;
      nlo += D[i] < P;
    }
  }

  void fusion(int *Lp, int *Rp, int *hi, int *lo, int &nhi, int &nlo) {
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  void _add(int b, int nb) {
    assert(b != -1);
    add_chain(b, nb);
  }

  void distribute_values(int b, int pivot, int chain[2]) {
    // fprintf(stderr, "has left\n");
    while (BUCKET(b)->N) {
      int i = !(BUCKET(b)->D[--BUCKET(b)->N] < pivot);
      // fprintf(stderr, "proc left %d, i = %d\n", N, i);
      leaf_insert(chain[i], BUCKET(b)->D[BUCKET(b)->N]);
    }
  }

  void leaf_split(int b, int &promotedValue, int &new_bucket) {
    // assert(leaf_check());
    assert(BUCKET(b)->next != -1);
    assert(BUCKET(b)->cap == INTERNAL_BSIZE); // The first bucket must be the smallest capacity.
    new_bucket = -1;

    if (BUCKET(BUCKET(b)->next)->next == -1 && 0) {
      /*
      Bucket *b = detach_and_get_next(); b->detach_and_get_next();

      if (N + b->N <= INTERNAL_BSIZE) {
        for (int i = 0; i < b->N; i++)
          D[N++] = b->D[i];
      } else {
        // assert(b->cap == cap);

        // Ensure both have at least 5 elements.
        assert(N >= 5 || b->N >= 5);
        while (N < 5) D[N++] = b->D[--b->N];
        assert(N >= 5 || b->N >= 5);
        while (b->N < 5) b->D[b->N++] = D[--N];
        assert(N >= 5 && b->N >= 5);

        int R[5];
        Random rng(140384); // TODO: use randomized seed.
        for (int i = 0; i < 3; i++) {
          assert(N > 0);
          int j = rng.nextInt(N);
          R[i] = D[j];
          D[j] = D[--N];
        }
        for (int i = 0; i < 2; i++) {
          assert(b->N > 0);
          int j = rng.nextInt(b->N);
          R[i + 3] = b->D[j];
          b->D[j] = b->D[--b->N];
        }

        for (int i = 0; i < 5; i++) {
          // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
        }

        std::nth_element(R, R + 2, R + 5);
        int pivot = R[2];
        D[N++] = R[0];
        D[N++] = R[1];
        D[N++] = R[3];
        b->D[b->N++] = R[4];

        Bucket *nb = transfer_to(new_leaf(parent, INTERNAL_BSIZE), pivot);
        b->transfer_to(nb, pivot);
        for (int i = 0; i < b->N; i++) {
          leaf_insert(b->D[i]);
        }
        promotedValue = pivot;
        new_bucket = nb;
      }

      delete_leaf(b);
      */
    } else {
      // fprintf(stderr, "split N = %d\n", BUCKET(b)->N);
      // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
      assert(BUCKET(b)->N + BUCKET(BUCKET(b)->tail)->N >= 11);
      while (BUCKET(b)->N < 11) {
        BUCKET(b)->D[BUCKET(b)->N++] = BUCKET(BUCKET(b)->tail)->D[--BUCKET(BUCKET(b)->tail)->N];
      }
      // fprintf(stderr, "tail N = %d\n", BUCKET(b)->N);
      int R[11];
      Random rng(140384); // TODO: use randomized seed.
      for (int i = 0; i < 11; i++) {
        assert(BUCKET(b)->N > 0);
        int j = rng.nextInt(BUCKET(b)->N);
        R[i] = BUCKET(b)->D[j];
        BUCKET(b)->D[j] = BUCKET(b)->D[--BUCKET(b)->N];
      }
      assert(BUCKET(b)->N >= 0);

      // fprintf(stderr, "random N = %d\n", BUCKET(b)->N);

      int Nb = b;

      // Replace elements with gradually decreasing probability.
      for (int i = 1; BUCKET(Nb)->next != -1; i++) {
        Nb = BUCKET(Nb)->next;
        assert(i > 0);
        int j = rng.nextInt(i);
        if (j < 11) {
          assert(BUCKET(Nb)->N > 0);
          int k = rng.nextInt(BUCKET(Nb)->N);
          // fprintf(stderr, "swap %d  <>  %d,   %d %d\n", R[j], Nb->D[k], j, k);
          swap(R[j], BUCKET(Nb)->D[k]);
        }
      }

      for (int i = 0; i < 11; i++) {
        // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      }

      std::nth_element(R, R + 5, R + 11);
      int pivot = R[5];
      R[5] = R[10];
      for (int i = 0; i < 10; i++) {
        BUCKET(b)->D[BUCKET(b)->N++] = R[i];
      }
      BUCKET(b)->D[BUCKET(b)->N++] = BUCKET(Nb)->D[--BUCKET(Nb)->N];
      // fprintf(stderr, "pivot = %d\n", pivot);
      // fprintf(stderr, "split3 N = %d, pivot = %d\n", next->N, pivot);

      // debug(10);

      new_bucket = bucket_allocator.alloc(true);
      BUCKET(new_bucket)->init(BUCKET(b)->parent, INTERNAL_BSIZE);
      int chain[2] { b, new_bucket };

      // Split the first bucket (this bucket).
      for (int i = 0; i < BUCKET(b)->N; i++) {
        // fprintf(stderr, "hihi %d < %d\n", i, N);
        if (BUCKET(b)->D[i] >= pivot) {
          leaf_insert(chain[1], BUCKET(b)->D[i]);
          BUCKET(b)->D[i--] = BUCKET(b)->D[--BUCKET(b)->N];
        }
      }

      Nb = BUCKET(b)->detach_and_get_next();

      int Lb = -1, Rb = -1;
      // TODO: optimize locality.
      int hi[LEAF_BSIZE], nhi = 0;
      int lo[LEAF_BSIZE], nlo = 0;
      while (true) {
        if (nhi && nlo) {
          assert(Lb != -1 && Rb != -1);
          fusion(BUCKET(Lb)->D, BUCKET(Rb)->D, hi, lo, nhi, nlo);
          if (!nhi) { _add(chain[0], Lb); Lb = -1; }
          if (!nlo) { _add(chain[1], Rb); Rb = -1; }
        } else if (Lb == -1) {
          if (Nb == -1) break;
          Lb = Nb;
          Nb = BUCKET(Nb)->detach_and_get_next();
          if (!BUCKET(Lb)->is_full()) break;
        } else if (!nhi) {
          assert(Lb != -1);
          mark_hi(BUCKET(Lb)->D, BUCKET(Lb)->N, pivot, hi, nhi);
          if (!nhi){ _add(chain[0], Lb); Lb = -1; }
        } else if (Rb == -1) {
          if (Nb == -1) break;
          Rb = Nb;
          Nb = BUCKET(Nb)->detach_and_get_next();
          if (!BUCKET(Rb)->is_full()) break;
        } else if (!nlo) {
          assert(Rb != -1);
          mark_lo(BUCKET(Rb)->D, BUCKET(Rb)->N, pivot, lo, nlo);
          if (!nlo){ _add(chain[1], Rb); Rb = -1; }
        } else {
          assert(0);
        }
      }
      assert(Nb == -1);

      // fprintf(stderr, "splited\n");
      if (Lb != -1) distribute_values(Lb, pivot, chain), bucket_allocator.destroy(Lb);
      if (Rb != -1) distribute_values(Rb, pivot, chain), bucket_allocator.destroy(Rb);
      promotedValue = pivot;
    }
    // assert(leaf_check());
  }

  void internal_insert(int b, int value, int nb) {
    // fprintf(stderr, "ii %d\n", b);
    int *C = CHILDREN(b);
    assert(!BUCKET(b)->is_leaf());
    assert(!BUCKET(b)->is_full());
    int i = BUCKET(b)->N - 1;
    while (i >= 0 && BUCKET(b)->D[i] > value) {
      BUCKET(b)->D[i + 1] = BUCKET(b)->D[i];
      C[i + 2] = C[i + 1];
      i--;
    }
    BUCKET(b)->D[i + 1] = value;
    C[i + 2] = nb;
    BUCKET(b)->N++;
    BUCKET(nb)->parent = b;
    // assert(check());
  }

  void internal_erase(int b, int *C, int pos, int stride) {
    BUCKET(b)->N--;
    while (pos < BUCKET(b)->N) {
      BUCKET(b)->D[pos] = BUCKET(b)->D[pos + 1];
      C[pos + stride] = C[pos + stride + 1];
      pos++;
    }
    if (!stride) C[pos] = C[pos + 1];
  }

  bool shift_leaves(int b, int pos) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(BUCKET(L)->next == -1);
    assert(BUCKET(R)->next == -1);

    // Move from M to L as many as possible.
    bool changed = false;
    while (!BUCKET(L)->is_full() && BUCKET(R)->N) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      BUCKET(b)->D[pos] = BUCKET(R)->leaf_promote_first();
      changed = true;
    }
    if (!BUCKET(L)->is_full() && !BUCKET(R)->N) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      internal_erase(b, CHILDREN(b), pos, 1);
      bucket_allocator.destroy(R);
      return true; // R is empty, compaction is done.
    }
    return changed;
  }

  bool shift_internals(int b, int pos) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(BUCKET(L)->next == -1);
    assert(BUCKET(R)->next == -1);

    // Move from R to L as many as possible.
    bool changed = false;
    while (!BUCKET(L)->is_full() && BUCKET(R)->N) {
      internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
      BUCKET(b)->D[pos] = BUCKET(R)->internal_promote_first(CHILDREN(R));
      changed = true;
    }
    if (!BUCKET(L)->is_full() && !BUCKET(R)->N) {
      internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
      internal_erase(b, CHILDREN(b), pos, 1);
      bucket_allocator.destroy(R);
      return true; // R is empty, compaction is done.
    }
    return changed;
  }

  int internal_split(int b) {
    int nb = bucket_allocator.alloc(true);
    int nbc = child_allocator.alloc(true);
    BUCKET(nb)->init(BUCKET(b)->parent, LEAF_BSIZE, nbc);
    int *C = CHILDREN(b);
    int *nbC = CHILDREN(nb);
    CHILDREN(nb)[0] = C[BUCKET(b)->N / 2];
    for (int i = BUCKET(b)->N / 2, j = 0; i < BUCKET(b)->N; i++) {
      BUCKET(nb)->D[j++] = BUCKET(b)->D[i];
      nbC[j] = C[i + 1];
      BUCKET(nbC[j])->parent = nb;
    }
    BUCKET(b)->N /= 2;
    BUCKET(nb)->N = BUCKET(b)->N;
    // assert(check());
    return nb;
  }

  bool split_chain(int b) {
    if (BUCKET(b)->next == -1) return false;
    assert(!locked);

    int promotedValue;
    int nb;
    // fprintf(stderr, "split_chain %d, %d\n", BUCKET(b)->N, BUCKET(b)->next);
    leaf_split(b, promotedValue, nb);
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);

    int parent = BUCKET(b)->parent;
    while (parent != -1 && nb != -1) {
      if (BUCKET(parent)->is_full()) {
        // fprintf(stderr, "parful\n");
        int inb = internal_split(parent);
        int promotedValueInternal = BUCKET(parent)->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          internal_insert(inb, promotedValue, nb);
        } else {
          internal_insert(parent, promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = BUCKET(nb)->parent;
      } else {
        // fprintf(stderr, "internal\n");
        internal_insert(parent, promotedValue, nb);
        nb = -1;
      }
    }
    if (nb != -1) {
      // Replace root
      // fprintf(stderr, "OLD ROOT %d\n", root);
      assert(parent == -1);
      int new_root = bucket_allocator.alloc(true);
      int c = child_allocator.alloc(true);
      BUCKET(new_root)->init(-1, INTERNAL_BSIZE, c); // New internal bucket.
      CHILDREN(new_root)[0] = root;
      BUCKET(root)->parent = new_root;
      internal_insert(new_root, promotedValue, nb);
      root = new_root;
      // fprintf(stderr, "NEW ROOT %d\n", root);
    }
    // fprintf(stderr, "\n\n\ndone split %d\n", b);
    // debug();
    return true;
  }

  bool optimize(int b = -1) {
    if (b == -1) {
      if (root == -1) return false;
      while (optimize(root)) {
        if (BUCKET(root)->N == 0) {
          fprintf(stderr, "PROMOTE ROOT\n");
          assert(!BUCKET(root)->is_leaf());
          int c = child(root, 0);
          bucket_allocator.destroy(root);
          root = c;
          BUCKET(root)->parent = -1;
        }
        fprintf(stderr, "optimize\n");
      }
      return false;
    }

    if (BUCKET(b)->is_leaf()) {
      // fprintf(stderr, "split_chain %d\n", b);
      bool splitted = split_chain(b);
      // fprintf(stderr, "split_chain %d\n", b);
      while (split_chain(b));
      // fprintf(stderr, "split_chain %d %d\n", b, splitted);
      return splitted;
    }

    bool changed = false;
    // fprintf(stderr, "internal %d\n", b);
    for (int i = 0; i <= BUCKET(b)->N; i++) {
      if (optimize(child(b, i))) {
        i = -1;
        changed = true;
      }
    }
    for (int i = 0; i < BUCKET(b)->N; i++) {
      int L = child(b, i);
      int R = child(b, i + 1);
      assert(BUCKET(L)->is_leaf() == BUCKET(R)->is_leaf());
      if (BUCKET(L)->is_leaf()) {
        if (shift_leaves(b, i)) changed = 1;
      } else {
        if (shift_internals(b, i)) changed = 1;
      }
    }
    // fprintf(stderr, "internal %d %d\n", b, changed);
    return changed;
  }

  void debug_data(int b) {
    if (BUCKET(b)->is_leaf()) {
      fprintf(stderr, "N = %d (p=%d, LEAF), ", BUCKET(b)->N, BUCKET(b)->P);
    } else {
      fprintf(stderr, "N = %d, ", BUCKET(b)->N);
    }
    fprintf(stderr, "[");
    for (int i = 0; i < BUCKET(b)->N; i++) {
      fprintf(stderr, "%d ", BUCKET(b)->D[i]);
    }
    fprintf(stderr, "]");
  }

  int debug(int b = -1, int depth = 0) {
    for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
    if (b == -1) b = root;
    int sz = 0;
    if (BUCKET(b)->is_leaf()) {
      while (true) {
        sz += BUCKET(b)->N;
        debug_data(b);
        if (BUCKET(b)->next == -1) break;
        fprintf(stderr, "------ ");
        b = BUCKET(b)->next;
      }
      fprintf(stderr, "\n");
    } else {
      debug_data(b);
      fprintf(stderr, "\n");
      for (int i = 0; i <= BUCKET(b)->N; i++) {
        sz += debug(CHILDREN(b)[i], depth + 1);
      }
      for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
      fprintf(stderr, "size = %d\n", sz);
    }
    return sz;
  }

  int max_depth() {
    return 1;
  }

  int slack() {
    return 1;
  }

  int size() {
    int ret = 0;
    int b = root;
    while (b != -1) {
      Bucket *B = BUCKET(b);
      ret += B->N;
      b = B->next;
    }
    return ret;
  }

  pair<bool, int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    return lower_bound_rec(root, value);
  }

  pair<bool, int> lower_bound_rec(int b, int value) {
    if (BUCKET(b)->is_leaf()) {
      assert(!split_chain(b));
      int pos = BUCKET(b)->leaf_lower_pos(value);
      return  (pos < BUCKET(b)->N) ? make_pair(true, BUCKET(b)->D[pos]) : make_pair(false, 0);
    } else {
      int pos = BUCKET(b)->internal_lower_pos(value);
      if (pos < BUCKET(b)->N && BUCKET(b)->D[pos] == value) {
        return make_pair(true, value); // Found in the internal bucket.
      }
      auto res = lower_bound_rec(CHILDREN(b)[pos], value);
      if (!res.first && pos < BUCKET(b)->N) {
        res = make_pair(true, BUCKET(b)->D[pos]);
      }
      return res;
    }
  }

  void add_chain(int &head, int next) {
    assert(BUCKET(next)->is_leaf());
    if (head == -1) {
      head = next;
    } else {
      Bucket *b = BUCKET(head);
      assert(b->is_leaf());
      if (b->next == -1) {
        b->next = b->tail = next;
      } else {
        BUCKET(b->tail)->next = next;
        b->tail = next;
      }
    }
  }

  void batch_insert(int *arr, int N) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0;
    while (i + INTERNAL_BSIZE <= N) {
      int idx = bucket_allocator.alloc(false);
      BUCKET(idx)->init(-1, INTERNAL_BSIZE);
      for (int j = 0; j < INTERNAL_BSIZE; j++) {
        bool ok = BUCKET(idx)->append(arr[i++]);
        assert(ok);
      }
      // fprintf(stderr, "chain %d\n", i);
      add_chain(root, idx);
    }
    // fprintf(stderr, "done %d\n", i);
    while (i < N) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "done %d\n", i);
  }

  void insert(int value) {
    insert(root, value);
  }

  void insert(int &b, int value) {
    // fprintf(stderr, "ins %d\n", value);
    // if (value == 711)  debug();

    if (b == -1) {
      b = bucket_allocator.alloc(true);
      BUCKET(b)->init(-1, LEAF_BSIZE);
      BUCKET(b)->append(value);
      return;
    }

    while (true) {
      Bucket *B = BUCKET(b);
      if (B->is_leaf()) break;
      // int *C = *child_allocator.get(B->C);
      // b = B->child(value, child_allocator);
      assert(0);
    }

    leaf_insert(b, value);
    // root->debug(0);
  }

  void leaf_insert(int b, int value) {
    if (BUCKET(b)->append(value)) return;
    // assert(BUCKET(b)->cap == INTERNAL_BSIZE);
    // assert(BUCKET(BUCKET(b)->tail)->next == -1);
    int tail = BUCKET(b)->tail;
    if (tail == -1 || BUCKET(tail)->is_full()) {
      tail = bucket_allocator.alloc(false);
      Bucket *B = BUCKET(b);
      Bucket *nb = BUCKET(tail);
      nb->init(B->parent, INTERNAL_BSIZE);
      nb->append(value);

      if (B->tail == -1) {
        B->next = B->tail = tail;
      } else {
        BUCKET(B->tail)->next = tail;
        B->tail = tail;
      }
    } else {
      BUCKET(tail)->append(value);
    }
  }

  bool erase(int value) {
    return true;
  }

  bool check() {
    return true;
  }
};


/*

Random rng;



void optimize() {
  if (is_leaf()) {
    ((Bucket*) this)->leaf_optimize();
  } else {
    for (int i = 0; i <= N; i++) {
      ((Bucket*) this)->child(i)->optimize();
    }
  }
}

bool check(int lo, int hi) {
  if (is_leaf()) return ((Bucket*) this)->leaf_check(lo, true, hi, true);
  return check(lo);
}

bool check(int lo) {
  assert(N >=0 && N <= cap);
  assert(is_leaf() || ((Bucket*) this)->child(N));
  if (is_leaf()) return ((Bucket*) this)->leaf_check();
  Bucket *ib = (Bucket*) this;
  if (N) {
    if (ib->child(0)->get_parent() != ib) return false;
    ib->child(0)->check(lo, D[0]);
  }
  for (int i = 0; i < N; i++) {
    assert(ib->child(i));
    if (i > 0) assert(ib->data(i - 1) < ib->data(i));
    if (ib->child(i + 1)->get_parent() != ib) return false;
    if (!ib->child(i + 1)->check(D[i], (i + 1 < N) ? D[i + 1] : 2147483647))
      return false;
  }
  return true;
}

pair<bool,int> leaf_erase_largest() {
  assert(!next);
  if (N == 0) return make_pair(false, D[0]);
  int pos = 0;
  int largest_pos = pos++;
  while (pos < N) {
    if (D[pos] > D[largest_pos])
      largest_pos = pos;
    pos++;
  }
  // fprintf(stderr, "pos %d %d\n", pos, largest_pos);
  swap(D[largest_pos], D[--N]);
  // assert(leaf_check());
  return make_pair(true, D[N]);
}

bool leaf_erase(int &v) {
  // assert(leaf_check());
  int pos = 0;
  while (pos < N) {
    if (D[pos] == v) {
      D[pos] = D[--N];
      return true;
    }
  }
  return false;
}

bool leaf_debug(const char *msg, int i, int j) {
  return false;
}

bool leaf_check() {
  return leaf_check(D[0],false,D[0],false);
}

bool leaf_check(int lo, bool useLo, int hi, bool useHi) {
  if (useLo) for (int i = 0; i < size(); i++) if ((D[i] < lo)) {
    fprintf(stderr,"D[%d] = %d, lo = %d\n", i, D[i], lo);
    return leaf_debug("useLo failed", i, 0);
  }
  if (useHi) for (int i = 0; i < size(); i++) if ((D[i] > hi)) {
    fprintf(stderr,"D[%d] = %d, hi = %d\n", i, D[i], hi);
    return leaf_debug("useHi failed", i, 0);
  }
  return true;
}



Bucket* transfer_to(Bucket *b, int pivot) {
  int oldN = N;
  N = 0;
  Bucket *arr[2] { this, b };
  for (int i = 0; i < oldN; i++) {
    int j = !(D[i] < pivot);
    arr[j]->leaf_insert(D[i]);
  }
  return b;
}


int leaf_promote_last() {
  if (P > 0) leaf_optimize();
  return D[--N];
}

Bucket*& child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}


class CTree {
  Bucket *root;

 public:

  const char *version = "Exp LEAF_BSIZE 2048";

  CTree() {
    root = new_leaf(NULL, INTERNAL_BSIZE);
  }

  void debug() {
    root->debug(0);
    // fprintf(stderr, "\n");
  }


  int max_depth(Bucket *b = NULL) {
    if (!b) b = root;
    if (b->is_leaf()) return 1;
    int ret = -1;
    for (int i = 0; i <= b->size(); i++) {
      int d = max_depth(((Bucket*) b)->child(i));
      assert(ret == -1 || ret == d);
      ret = d;
    }
    return ret + 1;
  }

  int slack(Bucket *b = NULL, int last = 0) {
    if (!b) b = root;
    int ret = b->capacity() - b->size();
    // if (ret > 10) fprintf(stderr, "slack = %d, for leaf = %d, last = %d\n", ret, b->is_leaf(), last);
    if (b->is_leaf()) return ret;
    for (int i = 0; i <= b->size(); i++) {
      ret += slack(((Bucket*) b)->child(i), i >= b->size());
    }
    return ret;
  }

  double t1 = 0, t2 = 0, t3 = 0;

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
      delete_leaf(M);
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
    delete_leaf(M);
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

  pair<Bucket*, int> find_bucket(int value, bool include_internal) {
    Bucket *b = root;
    while (true) {
      if (b->is_leaf()) {
        Bucket *Lb = (Bucket*) b;
        if (!split_chain(Lb)) break;
        assert(Lb->get_parent());
        b = Lb->get_parent();
      } else {
        Bucket *ib = (Bucket*) b;
        int pos = ib->internal_lower_pos(value);
        if (include_internal && pos < ib->size() && ib->data(pos) == value) {
          return make_pair(ib, pos); // Found in the internal bucket.
        }
        // bool changed = 0;
        // for (int i = 1; i <= ib->size(); i++) {
        //   Bucket *L = ib->child(i - 1);
        //   Bucket *R = ib->child(i);
        //   assert(L->is_leaf() == R->is_leaf());
        //   if (L->is_leaf()) {
        //     if (!(((Bucket*) L)->next_bucket() || ((Bucket*) R)->next_bucket()))
        //       if (shift_leaves(ib, i)) changed = 1;
        //   } else {
        //     if (shift_internals(ib, i)) changed = 1;
        //   }
        // }
        // if (changed) b = ib;
        // else 
        b = ib->child(pos);    // Search the child.

        // if (pos > 0 && pos < ib->size()) {
        //   Bucket *L = ib->child(pos - 1);
        //   Bucket *R = ib->child(pos + 1);
        //   assert(b->is_leaf() == L->is_leaf());
        //   assert(b->is_leaf() == R->is_leaf());
        //   if (b->size() + L->size() + R->size() + 1 < INTERNAL_BSIZE * 2) {
        //     if (b->is_leaf()) {
        //       if (!(((Bucket*) L)->next_bucket() || ((Bucket*) b)->next_bucket() || ((Bucket*) R)->next_bucket()))
        //         if (compact_leaves(ib, pos)) b = ib;
        //     } else {
        //       if (compact_internals(ib, pos)) b = ib;
        //     }
        //   }
        // }
      }
    }
    return make_pair(b, 0);
  }

  pair<bool, int> lower_bound2(int value) {
    if (root->size() == 0) {
      // fprintf(stderr, "PROMOTE ROOT\n");
      assert(!root->is_leaf());
      Bucket *x = ((Bucket*) root)->child(0);
      delete ((Bucket*) root);
      root = x;
      root->set_parent(NULL);
    }

    // fprintf(stderr, "lower_bound %d\n", value);
    pair<Bucket*, int> p = find_bucket(value, true);

    // Found in internal bucket.
    pair<bool, int> ret = make_pair(false, 0);
    if (!p.first->is_leaf()) {
      ret = make_pair(true, value);
    } else {
      Bucket *b = (Bucket*) p.first;
      int pos = b->leaf_lower_pos(value);
      if (pos < b->size()) {
        ret = make_pair(true, b->data(pos)); 
      } else {
        Bucket *ib = (Bucket*) b->get_parent();
        while (ib) {
          pos = ib->internal_lower_pos(value);
          if (pos < ib->size()) {
            ret = make_pair(true, ib->data(pos));
            break;
          }
          ib = (Bucket*) ib->get_parent();
        }
      }
    }

    return ret;
  }


  pair<bool, int> erase_largest(Bucket *b) {
    assert(b->is_leaf());
    auto res = ((Bucket*) b)->leaf_erase_largest();
    if (res.first) return res;
    while (b) {
      assert(!b->size());
      Bucket *parent = (Bucket*) b->get_parent();
      if (b->is_leaf()) {
        delete_leaf((Bucket*) b);
      } else {
        delete b;
      }
      if (parent->size() > 0) return make_pair(true, parent->internal_promote_last());
      b = parent;
    }
    return res;
  }

  bool erase(int value) {
    // assert(check());
    // fprintf(stderr, "ERASE %d\n", value);
    // debug();

    pair<Bucket*, int> p = find_bucket(value, true);
    if (p.first->is_leaf()) {
      return ((Bucket*) p.first)->leaf_erase(value); // May return false.
    }

    Bucket* upper = find_bucket(value, false).first;

    // The buckets may have moved.
    p = find_bucket(value, true);

    Bucket *ib = (Bucket*) p.first;
    int pos = p.second;
    assert(ib->data(pos) == value);

    auto res = erase_largest(upper);
    // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
    if (!res.first) {
      ib->internal_erase(pos, 0);
    } else {
      ib->set_data(pos, res.second);
    }
    // assert(check());
    return true;
  }

  bool check() {
    return root->check(-2147483648);
  }
};
*/
}


#endif
