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

namespace ctree {

#define INTERNAL_BSIZE        256  // Must be power of two.
#define LEAF_BSIZE            256  // Must be power of two.
#define LEAF_CHAINED_BSIZE  2048  // Must be power of two.

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
    D = new T[cap = 1];
    N = 0;
  }

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
    } else {
      int idx = free_indices.top();
      free_indices.pop();
      return idx;
    }
  }

  void destroy(int idx) {
    free_indices.push(idx);
  }

  T* get(int idx) {
    return &D[idx];
  }
};


struct Bucket {
  int D;       // Pointer to data_allocator (if cap == LEAF_BSIZE), or chained_data_allocator otherwise.
  int C;       // Pointer to the child bucket in child_allocator.
  int P;       // Pending insert if positive or pending delete if negative.
  int N;       // Number of data elements in this bucket pointed by D.
  int cap;     // The capacity of this bucket pointed by D.
  int parent;  // Pointer to the parent bucket in bucket_allocator.
  int next;    // Pointer to the next chained bucket in bucket_allocator.
  int tail;    // Pointer to the last chained bucket in bucket_allocator.
};

Allocator<Bucket> bucket_allocator;
Allocator<int[LEAF_BSIZE]> data_allocator;
Allocator<int[INTERNAL_BSIZE + 1]> child_allocator;
int nLeaves, nInternals, nCap, nDes, locked;

void init(int b, int parent, int cap) {
  bucket_allocator.get(b)->P = 0;
  bucket_allocator.get(b)->N = 0;
  bucket_allocator.get(b)->cap = cap;
  bucket_allocator.get(b)->parent = parent;
  bucket_allocator.get(b)->next = -1;
  bucket_allocator.get(b)->tail = -1;

  nCap += cap;
  nLeaves++;
}

void init_leaf(int b, int parent, int cap) {
  if (bucket_allocator.get(b)->cap == LEAF_BSIZE) {
    bucket_allocator.get(b)->D = data_allocator.alloc();
  } else {
    bucket_allocator.get(b)->D = data_allocator.alloc();
    // D = chained_data_allocator.alloc();
  }
  bucket_allocator.get(b)->C = -1;
  init(b, parent, cap);
}

void init_internal(int b, int parent, int cap) {
  bucket_allocator.get(b)->D = data_allocator.alloc();
  bucket_allocator.get(b)->C = child_allocator.alloc();
  init(b, parent, cap);
}

void destroy(int b) {
  nCap -= bucket_allocator.get(b)->cap;
  nLeaves--;
  nDes++;
}


bool is_full(Bucket *b) { return b->N == b->cap; }
bool is_leaf(Bucket *b) { return b->C == -1; }

int next_bucket(Bucket* b) {
  assert(is_leaf(b));
  return b->next;
}

int child(Bucket *b, int i) {
  assert(!is_leaf(b));
  assert(i >= 0 && i <= b->N);
  return (*child_allocator.get(b->C))[i];
}

int data(Bucket *b, int i) {
  assert(i >= 0 && i < b->N);
  return (*data_allocator.get(b->D))[i];
}

void set_data(Bucket *b, int i, int value) {
  assert(i >= 0 && i < b->N);
  (*data_allocator.get(b->D))[i] = value;
}

void leaf_insert(int b, int value) {
  // assert(leaf_check());

  assert(is_leaf(bucket_allocator.get(b)));
  assert(bucket_allocator.get(b)->N >= 0);
  if (!is_full(bucket_allocator.get(b))) {
    (*data_allocator.get(bucket_allocator.get(b)->D))[bucket_allocator.get(b)->N++] = value;
    bucket_allocator.get(b)->P++;
  } else {
    if (bucket_allocator.get(b)->tail == -1) {
      assert(bucket_allocator.get(b)->cap == INTERNAL_BSIZE);
      int idx = bucket_allocator.alloc();
      init_leaf(idx, bucket_allocator.get(b)->parent, INTERNAL_BSIZE);
      bucket_allocator.get(b)->next =
      bucket_allocator.get(b)->tail = idx;
    } else {
      if (is_full(bucket_allocator.get(bucket_allocator.get(b)->tail))) {
        assert(bucket_allocator.get(bucket_allocator.get(b)->tail)->next == -1);
        int idx = bucket_allocator.alloc();
        init_leaf(idx, bucket_allocator.get(b)->parent, LEAF_BSIZE); // Doubling?
        bucket_allocator.get(bucket_allocator.get(b)->tail)->next = idx;
        bucket_allocator.get(b)->tail = idx;
      }
    }
    int tail = bucket_allocator.get(b)->tail;
    Bucket *b = bucket_allocator.get(tail);
    (*data_allocator.get(b->D))[b->N++] = value;
  }
  // assert(leaf_check());
}


class CTree {
  int root;

 public:

  const char *version = "Exp LEAF_BSIZE 2048";

  CTree() {
    root = bucket_allocator.alloc();
    init_leaf(root, -1, LEAF_BSIZE);
  }

  bool optimize() {
    return 0;
  }

  int max_depth() {
    return 1;
  }

  int slack() {
    return 1;
  }

  pair<bool, int> lower_bound(int value) {
    // Bucket *b = bucket_allocator.get(root);
    return make_pair(0,0);
  }

  void insert(int value) {
    // fprintf(stderr, "ins %d\n", value);
    // if (value == 711)  debug();
    int b = root;
    while (!is_leaf(bucket_allocator.get(b)))
      b = child(bucket_allocator.get(b), value);
    leaf_insert(b, value);
    // root->debug(0);
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

int debug(int depth) {
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  Bucket *b = (Bucket *) this;
  int sz = 0;
  while (b) {
    sz += b->size();
    ((Bucket*) b)->leaf_debug();
    if ((b = b->next_bucket())) {
      fprintf(stderr, "------ ");
    }
  }
  fprintf(stderr, "\n");
  if (!is_leaf()) {
    for (int i = 0; i <= ((Bucket*) this)->size(); i++) {
      sz += ((Bucket*) this)->child(i)->debug(depth + 1);
    }
  }
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  fprintf(stderr, "size = %d\n", sz);
  return sz;
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


void leaf_debug() {
  fprintf(stderr, "N = %d (p=%d, LEAF), ", N, P);
  for (int i = 0; i < N; i++) {
    fprintf(stderr, "%d ", D[i]);
  }
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

Bucket* detach_and_get_next() {
  Bucket* ret = next;
  next = tail = NULL;
  return ret;
}

void _add(Bucket *&b, Bucket *nb) {
  if (!b) b = nb;
  else b->add_chain(nb);
}

void distribute_values(int pivot, Bucket *chain[2]) {
  // fprintf(stderr, "has left\n");
  while (N) {
    int i = !(D[--N] < pivot);
    // fprintf(stderr, "proc left %d, i = %d\n", N, i);
    chain[i]->leaf_insert(D[N]);
  }
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

void leaf_split(int &promotedValue, Bucket *&new_bucket) {
  // assert(leaf_check());
  assert(next);
  assert(cap == INTERNAL_BSIZE); // The first bucket must be the smallest capacity.
  new_bucket = NULL;

  if (!next->next) {
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

  } else {
    // fprintf(stderr, "split N = %d\n", N);
    // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
    assert(N + tail->N >= 11);
    while (N < 11) {
      D[N++] = tail->D[--tail->N];
    }
    int R[11];
    Random rng(140384); // TODO: use randomized seed.
    for (int i = 0; i < 11; i++) {
      assert(N > 0);
      int j = rng.nextInt(N);
      R[i] = D[j];
      D[j] = D[--N];
    }
    assert(N >= 0);

    Bucket *Nb = this;

    // Replace elements with gradually decreasing probability.
    for (int i = 1; Nb->next; i++) {
      Nb = Nb->next;
      assert(i > 0);
      int j = rng.nextInt(i);
      if (j < 11) {
        assert(Nb->N > 0);
        int k = rng.nextInt(Nb->N);
        // fprintf(stderr, "swap %d  <>  %d,   %d %d\n", R[j], Nb->D[k], j, k);
        swap(R[j], Nb->D[k]);
      }
    }
    // fprintf(stderr, "split2 N = %d\n", q.size());

    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
    }

    std::nth_element(R, R + 5, R + 11);
    int pivot = R[5];
    R[5] = R[10];
    for (int i = 0; i < 10; i++) {
      D[N++] = R[i];
    }
    D[N++] = Nb->D[--Nb->N];
    // fprintf(stderr, "queue size = %lu\n", q.size());
    // fprintf(stderr, "split3 N = %d, pivot = %d\n", next->N, pivot);

    // debug(10);

    Bucket *chain[2] { this, new_leaf(parent, INTERNAL_BSIZE) };

    // Split the first bucket (this bucket).
    for (int i = 0; i < N; i++) {
      // fprintf(stderr, "hihi %d < %d\n", i, N);
      if (D[i] >= pivot) {
        chain[1]->leaf_insert(D[i]);
        D[i--] = D[--N];
      }
    }

    Nb = detach_and_get_next();

    Bucket *Lb = NULL, *Rb = NULL;
    // TODO: optimize locality.
    int hi[LEAF_BSIZE], nhi = 0;
    int lo[LEAF_BSIZE], nlo = 0;
    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        fusion(Lb->D, Rb->D, hi, lo, nhi, nlo);
        if (!nhi) { _add(chain[0], Lb); Lb = NULL; }
        if (!nlo) { _add(chain[1], Rb); Rb = NULL; }
      } else if (!Lb) {
        if (!Nb) break;
        Lb = Nb;
        Nb = Nb->detach_and_get_next();
        if (!Lb->is_full()) break;
      } else if (!nhi) {
        assert(Lb);
        mark_hi(Lb->D, Lb->N, pivot, hi, nhi);
        if (!nhi){ _add(chain[0], Lb); Lb = NULL; }
      } else if (!Rb) {
        if (!Nb) break;
        Rb = Nb;
        Nb = Nb->detach_and_get_next();
        if (!Rb->is_full()) break;
      } else if (!nlo) {
        assert(Rb);
        mark_lo(Rb->D, Rb->N, pivot, lo, nlo);
        if (!nlo){ _add(chain[1], Rb); Rb = NULL; }
      } else {
        assert(0);
      }
    }
    assert(!Nb);

    // fprintf(stderr, "splited\n");
    if (Lb) Lb->distribute_values(pivot, chain), delete_leaf(Lb);
    if (Rb) Rb->distribute_values(pivot, chain), delete_leaf(Rb);
    promotedValue = pivot;
    new_bucket = chain[1];
  }
  // assert(leaf_check());
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

int leaf_promote_last() {
  if (P > 0) leaf_optimize();
  return D[--N];
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

Bucket(Bucket *parent, Bucket *left_child) : Bucket(parent, INTERNAL_BSIZE) {
  nInternals++;
  P = -1;
  C[0] = left_child;
  left_child->set_parent(this);
}

~Bucket() {
  nInternals--;
}

Bucket* internal_split() {
  Bucket *nb = new Bucket(parent, C[N / 2]);
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
    nb->C[j] = C[i + 1];
    nb->C[j]->set_parent(nb);
  }
  N /= 2;
  nb->N = N;

  // for (int i = N + 1; i <= nb->capacity(); i++) nb->C[i] = NULL;
  // for (int i = N + 1; i <= cap; i++) C[i] = NULL;

  // assert(check());
  return nb;
}


void internal_insert(int value, Bucket *b, int left) {
  // assert(check());
  assert(!is_full());
  int i = N - 1;
  while (i >= 0 && D[i] > value) {
    D[i + 1] = D[i];
    C[i + 2] = C[i + 1];
    i--;
  }
  D[i + 1] = value;
  if (left == -1) {
    C[i + 2] = C[i + 1];
    C[i + 1] = b;
  } else {
    C[i + 2] = b;
  }
  N++;
  b->set_parent(this);
  // assert(check());
}

int internal_lower_pos(int value) {
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

int internal_promote_first() {
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

Bucket*& child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}

void internal_erase(int pos, int stride) {
  N--;
  while (pos < N) {
    D[pos] = D[pos + 1];
    C[pos + stride] = C[pos + stride + 1];
    pos++;
  }
  if (!stride) C[pos] = C[pos + 1];
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

  bool optimize(Bucket *b = NULL) {
    if (root->size() == 0) {
      // fprintf(stderr, "PROMOTE ROOT\n");
      assert(!root->is_leaf());
      Bucket *x = ((Bucket*) root)->child(0);
      delete ((Bucket*) root);
      root = x;
      root->set_parent(NULL);
    }

    if (!b) b = root;
    if (b->is_leaf()) {
      bool ok = split_chain((Bucket*) b);
      while (split_chain((Bucket*) b));
      return ok;
    }

    bool changed = false;
    Bucket *ib = (Bucket*) b;
    for (int i = 0; i <= ib->size(); i++) {
      if (optimize(ib->child(i))) {
        i = -1;
        changed = true;
      }
    }
    for (int i = 1; i <= ib->size(); i++) {
      Bucket *L = ib->child(i - 1);
      Bucket *R = ib->child(i);
      assert(L->is_leaf() == R->is_leaf());
      if (L->is_leaf()) {
        assert(!((Bucket*) L)->next_bucket());
        assert(!((Bucket*) R)->next_bucket());
        if (!(((Bucket*) L)->next_bucket() || ((Bucket*) R)->next_bucket())) {
          if (shift_leaves(ib, i)) changed = 1, i = 0;
        }
      } else {
        if (shift_internals(ib, i)) changed = 1, i = 0;
      }
    }
    return changed;
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

  bool split_chain(Bucket *b) {
    // fprintf(stderr, "split_chain %d, %d\n", b->size(), b->next_bucket());
    if (!b->next_bucket()) return false;
    assert(!locked);

    int promotedValue;
    Bucket *nb;
    b->leaf_split(promotedValue, nb);
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);
    Bucket *parent = (Bucket*) b->get_parent();

    while (parent && nb) {
      if (parent->is_full()) {
        // fprintf(stderr, "parful\n");
        Bucket *inb = parent->internal_split();
        int promotedValueInternal = parent->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          inb->internal_insert(promotedValue, nb);
        } else {
          parent->internal_insert(promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = (Bucket*) nb->get_parent();
      } else {
        // fprintf(stderr, "internal\n");
        parent->internal_insert(promotedValue, nb);
        nb = NULL;
        break;
      }
    }
    if (nb) {
      // Replace root
        // fprintf(stderr, "replace root\n");
      assert(parent == NULL);
      root = new Bucket(NULL, root);
      ((Bucket*) root)->internal_insert(promotedValue, nb);
    }
    return true;
  }

  bool shift_leaves(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    bool changed = false;
    while (!L->is_full() && M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->set_data(lpos, M->leaf_promote_first());
      changed = true;
    }
    if (!L->is_full() && !M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->internal_erase(lpos, 1);
      delete_leaf(M);
      return true; // M is empty, compaction is done.
    }
    return changed;
  }

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

  bool shift_internals(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    bool changed = false;
    while (!L->is_full() && M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->set_data(lpos, M->internal_promote_first());
      changed = true;
    }
    if (!L->is_full() && !M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->internal_erase(lpos, 1);
      delete M;
      return true; // M is empty, compaction is done.
    }
    return changed;
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

  pair<bool, int> lower_bound(int value) {
    return lower_bound_rec(root, value);
  }

  pair<bool, int> lower_bound_rec(Bucket *b, int value) {
    if (b->is_leaf()) {
      Bucket *Lb = (Bucket*) b;
      assert(!split_chain(Lb));

      int pos = Lb->leaf_lower_pos(value);
      return  (pos < b->size()) ? make_pair(true, b->data(pos)) : make_pair(false, 0);

    } else {
      Bucket *ib = (Bucket*) b;
      int pos = ib->internal_lower_pos(value);
      if (pos < ib->size() && ib->data(pos) == value) {
        return make_pair(true, value); // Found in the internal bucket.
      }
      auto res = lower_bound_rec(ib->child(pos), value);
      if (!res.first && pos < ib->size()) {
        res = make_pair(true, ib->data(pos));
      }
      return res;
    }
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
