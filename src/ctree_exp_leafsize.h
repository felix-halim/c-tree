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

#define INTERNAL_BSIZE 64   // Must be power of two.
#define MIN_LEAF_BSIZE 64     // Must be power of two.
#define MAX_LEAF_BSIZE 4096     // Must be power of two.

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

int nLeaves, nInternals, nCap, nDes, locked;
Random rng;

class Bucket {
 protected:
  int *D;
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N, cap;
  Bucket *parent;

 public:
  int get_cap() const { return cap; }
  int size() const { return N; }
  bool is_full() const { return size() == cap; }
  bool is_leaf() const { return pending_insert >= 0; }
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }
  void set_data(int i, int v) { D[i] = v; }
  Bucket* get_parent() const { return parent; }
  void set_parent(Bucket *parent) { this->parent = parent; }
  int debug(int depth);
  bool check(int lo) const;
  bool check(int lo, int hi) const;
  pair<bool,int> erase_largest();
};

class LeafBucket : public Bucket {
 protected:
  LeafBucket *next, *tail;  // Store pending inserts in a linked list.

 public:
  ~LeafBucket();
  LeafBucket(Bucket *parent, int cap);
  void init(int cap);
  LeafBucket* next_bucket() { return next; }

  // T nth(int idx) {
  //   assert(next() == -1);      // it doesn't make sense crack a chained bucket!
  //   int L, R, i = get_piece_by_index(idx, L, R);
  //   assert(L >= 0 && L <= R && R <= size());
  //   if (!piece_is_sorted(i)){    // sort the piece if it isn't
  //     std::sort(D+L,D+R);
  //     piece_set_sorted(i,true);
  //   }
  //   return D[idx];          // this is the correct nth element
  // }

  // T randomValue(Random &rng) const { return D[rng.nextInt(size())]; }

  bool leaf_erase(int &v);
  pair<bool,int> leaf_erase_largest();

  bool leaf_debug(const char *msg, int i, int j) const;

  // call this function to check the consistency of this Bucket structure
  bool leaf_check() const;
  bool leaf_check(int lo, bool useLo, int hi, bool useHi) const;

  void leaf_debug();
  void leaf_insert(int v);
  void leaf_split(int &promotedValue, LeafBucket *&nb);
  void leaf_optimize();
  int promote_last();
  int leaf_lower_pos(int value);
  LeafBucket* detach_and_get_next();
  void add_chain(LeafBucket *b);
  void distribute_values(int pivot, LeafBucket* chain[2]);
  LeafBucket* transfer_to(LeafBucket *b, int pivot);

  bool leaf_promote_first(int &pv) {
    pv = D[0];
    N--;
    for (int i = 0; i < N; i++)
      D[i] = D[i + 1];
    return N;
  }
};

class InternalBucket : public LeafBucket {
 protected:
  Bucket *C[INTERNAL_BSIZE + 1];

 public:
  ~InternalBucket();
  InternalBucket(Bucket *parent, Bucket *left_child);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  void set_child(int i, Bucket *c) { assert(i >= 0 && i <= N); C[i] = c; }
  int child_pos(int value);
  Bucket*& child_bucket(int value);
  int internal_lower_pos(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
  bool internal_erase(int &v);
  pair<bool,int> internal_erase_largest();
  void internal_erase_pos(int pos);
  bool internal_promote_first(int &pv, Bucket *&nb) {
    pv = D[0];
    nb = C[0];
    N--;
    for (int i = 0; i < N; i++) {
      D[i] = D[i + 1];
      C[i] = C[i + 1];
    }
    C[N] = C[N + 1];
    return N;
  }
};



vector<LeafBucket*> free_leaves[30];

LeafBucket* new_leaf(Bucket *parent, int cap) {
  // return new LeafBucket(parent, cap);
  for (int i = 2; ; i++) {
    if ((1 << i) == cap) {
      if (free_leaves[i].empty()) {
        free_leaves[i].push_back(new LeafBucket(parent, cap));
      }
      LeafBucket* leaf = free_leaves[i].back();
      leaf->init(cap);
      free_leaves[i].pop_back();
      return leaf;
    }
  }
}

void delete_leaf(LeafBucket *b) {
  // delete b; return;
  for (int i = 2; ; i++) {
    if ((1 << i) == b->get_cap()) {
      free_leaves[i].push_back(b);
      break;
    }
  }
}

pair<bool,int> Bucket::erase_largest() {
  return is_leaf()
    ? ((LeafBucket*) this)->leaf_erase_largest()
    : ((InternalBucket*) this)->internal_erase_largest();
}

int Bucket::debug(int depth) {
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  LeafBucket *b = (LeafBucket *) this;
  int sz = 0;
  while (b) {
    sz += b->size();
    ((LeafBucket*) b)->leaf_debug();
    if ((b = b->next_bucket())) {
      fprintf(stderr, "------ ");
    }
  }
  fprintf(stderr, "\n");
  if (!is_leaf()) {
    for (int i = 0; i <= ((LeafBucket*) this)->size(); i++) {
      sz += ((InternalBucket*) this)->child(i)->debug(depth + 1);
    }
  }
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  fprintf(stderr, "size = %d\n", sz);
  return sz;
}

bool Bucket::check(int lo, int hi) const {
  if (is_leaf()) return ((LeafBucket*) this)->leaf_check(lo, true, hi, true);
  return check(lo);
}

bool Bucket::check(int lo) const {
  assert(N >=0 && N <= cap);
  assert(is_leaf() || ((InternalBucket*) this)->child(N));
  if (is_leaf()) return ((LeafBucket*) this)->leaf_check();
  InternalBucket *ib = (InternalBucket*) this;
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


LeafBucket::LeafBucket(Bucket *parent, int cap) {
  this->parent = parent;
  this->cap = cap;
  D = new int[cap];
  init(cap);
}

void LeafBucket::init(int cap) {
  assert(this->cap == cap);
  N = 0;
  pending_insert = 0;
  next = tail = NULL;
  nCap += cap;
  nLeaves++;
}

LeafBucket::~LeafBucket() {
  nCap -= cap;
  nLeaves--;
  nDes++;
}

void LeafBucket::leaf_debug() {
  fprintf(stderr, "N = %d (p=%d, LEAF), ", N, pending_insert);
  for (int i = 0; i < N; i++) {
    fprintf(stderr, "%d ", D[i]);
  }
}

void LeafBucket::leaf_insert(int value) {
  // assert(leaf_check());
  assert(is_leaf());
  assert(N >= 0);
  pending_insert = 1;
  if (!is_full()) {
    D[N++] = value;
  } else {
    if (!tail) {
      assert(cap == MIN_LEAF_BSIZE);
      add_chain(new_leaf(parent, MIN_LEAF_BSIZE));
    } else if (tail->is_full()) {
      add_chain(new_leaf(parent, std::min(MAX_LEAF_BSIZE, MAX_LEAF_BSIZE)));
    }
    tail->D[tail->N++] = value;
  }
  // assert(leaf_check());
}


pair<bool,int> LeafBucket::leaf_erase_largest() {
  assert(!next);
  int pos = 0;
  int largest_pos = pos++;
  while (pos < N) {
    if (D[pos] > D[largest_pos])
      largest_pos = pos;
    pos++;
  }
  // fprintf(stderr, "pos %d %d\n", pos, largest_pos);
  swap(D[largest_pos], D[--N]);
  pending_insert = 1;
  // assert(leaf_check());
  return make_pair(true, D[N]);
}

bool LeafBucket::leaf_erase(int &v) {
  for (int i = 0; i < N; i++)
    if (D[i] == v) {
      swap(D[i], D[--N]);
      break;
    }
  pending_insert = 1;   // Adjust the pending index.
  // assert(leaf_check());
  return true;
}

bool LeafBucket::leaf_debug(const char *msg, int i, int j) const {
  return false;
}

bool LeafBucket::leaf_check() const {
  return leaf_check(D[0],false,D[0],false);
}

bool LeafBucket::leaf_check(int lo, bool useLo, int hi, bool useHi) const {
  if (useLo) for (int i = 0; i < size(); i++) if ((D[i] < lo)) {
    fprintf(stderr,"D[%d] = %d, lo = %d\n", i, D[i], lo);
    return leaf_debug("useLo failed", i, 0);
  }
  if (useHi) for (int i = 0; i < size(); i++) if (!(D[i] < hi)) {
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

LeafBucket* LeafBucket::detach_and_get_next() {
  LeafBucket* ret = next;
  next = tail = NULL;
  pending_insert = 1;
  return ret;
}

void LeafBucket::add_chain(LeafBucket *b) {
  if (!next || !tail) next = tail = b;
  else {
    tail->next = b;
    tail = b;
  }
}

void _add(LeafBucket *&b, LeafBucket *nb) {
  if (!b) b = nb;
  else b->add_chain(nb);
}

void LeafBucket::distribute_values(int pivot, LeafBucket *chain[2]) {
  // fprintf(stderr, "has left\n");
  while (N) {
    int i = !(D[--N] < pivot);
    // fprintf(stderr, "proc left %d, i = %d\n", N, i);
    chain[i]->leaf_insert(D[N]);
  }
}

LeafBucket* LeafBucket::transfer_to(LeafBucket *b, int pivot) {
  int oldN = N;
  N = 0;
  LeafBucket *arr[2] { this, b };
  for (int i = 0; i < oldN; i++) {
    int j = !(D[i] < pivot);
    arr[j]->leaf_insert(D[i]);
  }
  return b;
}

void LeafBucket::leaf_split(int &promotedValue, LeafBucket *&new_bucket) {
  // assert(leaf_check());
  assert(next);
  assert(cap == MIN_LEAF_BSIZE); // The first bucket must be the smallest capacity.
  new_bucket = NULL;

  if (!next->next) {
    LeafBucket *b = detach_and_get_next(); b->detach_and_get_next();

    if (N + b->N <= MIN_LEAF_BSIZE) {
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

      LeafBucket *nb = transfer_to(new_leaf(parent, MIN_LEAF_BSIZE), pivot);
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

    LeafBucket *Nb = this;

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

    LeafBucket *chain[2] { this, new_leaf(parent, MIN_LEAF_BSIZE) };

    // Split the first bucket (this bucket).
    for (int i = 0; i < N; i++) {
      // fprintf(stderr, "hihi %d < %d\n", i, N);
      if (D[i] >= pivot) {
        chain[1]->leaf_insert(D[i]);
        D[i--] = D[--N];
      }
    }

    Nb = detach_and_get_next();

    LeafBucket *Lb = NULL, *Rb = NULL;
    // TODO: optimize locality.
    int hi[MAX_LEAF_BSIZE], nhi = 0;
    int lo[MAX_LEAF_BSIZE], nlo = 0;
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

int LeafBucket::promote_last() {
  nth_element(D, D + N - 1, D + N);
  return D[--N];
}

void LeafBucket::leaf_optimize() {
  assert(pending_insert >= 0);
  // assert(!locked);
  sort(D, D + N);
  pending_insert = 0;
}

int LeafBucket::leaf_lower_pos(int value) {
  // assert(leaf_check());
  if (pending_insert) leaf_optimize();
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

InternalBucket::InternalBucket(Bucket *parent, Bucket *left_child) : LeafBucket(parent, INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
  C[0] = left_child;
  left_child->set_parent(this);
}

InternalBucket::~InternalBucket() {
  nInternals--;
}

InternalBucket* InternalBucket::internal_split() {
  InternalBucket *nb = new InternalBucket(parent, C[N / 2]);
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
    nb->C[j] = C[i + 1];
    nb->C[j]->set_parent(nb);
  }
  N /= 2;
  nb->N = N;
  // assert(check());
  return nb;
}


void InternalBucket::internal_insert(int value, Bucket *b) {
  // assert(check());
  assert(!is_full());
  int i = N - 1;
  while (i >= 0 && D[i] > value) {
    D[i + 1] = D[i];
    C[i + 2] = C[i + 1];
    i--;
  }
  D[i + 1] = value;
  C[i + 2] = b;
  N++;
  b->set_parent(this);
  // assert(check());
}

int InternalBucket::internal_lower_pos(int value) {
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

int InternalBucket::internal_promote_last() {
 return D[--N];
}

void InternalBucket::internal_erase_pos(int pos) {
  N--;
  while (pos < N) {
    D[pos] = D[pos + 1];
    C[pos] = C[pos + 1];
    pos++;
  }
  C[pos] = C[pos + 1];
}

bool InternalBucket::internal_erase(int &v) {
  int pos = internal_lower_pos(v);
  assert(D[pos] == v);
  auto res = C[pos]->erase_largest();
  // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
  if (!res.first) {
    N--;
    while (pos < N) {
      D[pos] = D[pos + 1];
      C[pos] = C[pos + 1];
      pos++;
    }
    C[pos] = C[pos + 1];
  } else {
    D[pos] = res.second;
  }
  // fprintf(stderr, "ii4");
  return true;
}

pair<bool,int> InternalBucket::internal_erase_largest() {
  auto res = C[N]->erase_largest();
  if (res.first) return res;
  if (C[N]->is_leaf()) {
    LeafBucket *b = (LeafBucket*) C[N];
    assert(!b->size());
    delete_leaf(b);
  } else {
    InternalBucket *b = (InternalBucket*) C[N];
    assert(!b->size());
    delete b;
  }
  if (N > 0) return make_pair(true, D[--N]);
  return make_pair(false, D[0]);
}

Bucket*& InternalBucket::child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}


class CTree {
  Bucket *root;

 public:

  const char *version = "Crack 2048";

  CTree() {
    root = new_leaf(NULL, MIN_LEAF_BSIZE);
  }

  void debug() {
    root->debug(0);
    // fprintf(stderr, "\n");
  }

  double t1 = 0, t2 = 0, t3 = 0;

  bool split_chain(LeafBucket *b) {
    // fprintf(stderr, "split_chain %d\n", b->size());
    if (!b->next_bucket()) return false;

    // assert(!locked);

    int promotedValue;
    LeafBucket *nb;
    b->leaf_split(promotedValue, nb);
    InternalBucket *parent = (InternalBucket*) b->get_parent();

    while (parent && nb) {
      if (parent->is_full()) {
        InternalBucket *inb = parent->internal_split();
        int promotedValueInternal = parent->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          inb->internal_insert(promotedValue, nb);
        } else {
          parent->internal_insert(promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = (InternalBucket*) nb->get_parent();
      } else {
        parent->internal_insert(promotedValue, nb);
        nb = NULL;
        break;
      }
    }
    if (nb) {
      // Replace root
      assert(parent == NULL);
      root = new InternalBucket(NULL, root);
      ((InternalBucket*) root)->internal_insert(promotedValue, nb);
    }
    return true;
  }

  bool compact(Bucket *b = NULL) {
    if (!b) b = root;
    if (b->is_leaf()) return true;
    InternalBucket *ib = (InternalBucket*) b;
    for (int i = 0; i <= b->size(); i++) {
      Bucket *c = ib->child(i);
      compact(c);
      if (!c->size()) {
        assert(!c->is_leaf()); // only internal can be empty.
        ib->set_child(i, ((InternalBucket*) c)->child(0));
        delete ((InternalBucket *) c);
      }
    }

    for (int i = 0; i < b->size(); i++) {
      Bucket *L = ((InternalBucket*) b)->child(i);
      Bucket *R = ((InternalBucket*) b)->child(i + 1);
      if (L->is_leaf() != R->is_leaf()) continue;
      if (L->is_leaf()) {
        while (!L->is_full() && R->size()) {
          int promotedValue;
          ((LeafBucket*) R)->leaf_promote_first(promotedValue);
          ((LeafBucket*) L)->leaf_insert(b->data(i));
          b->set_data(i, promotedValue);
        }
        if (!L->is_full() && !R->size()) {
          ((LeafBucket*) L)->leaf_insert(b->data(i));
          ((InternalBucket*) b)->internal_erase_pos(i--);
          delete_leaf((LeafBucket*) R);
        }
      } else {
        while (!L->is_full() && R->size()) {
          int promotedValue;
          Bucket *nb;
          ((InternalBucket*) R)->internal_promote_first(promotedValue, nb);
          ((InternalBucket*) L)->internal_insert(b->data(i), nb);
          b->set_data(i, promotedValue);
        }
        if (!L->is_full() && !R->size()) {
          ((InternalBucket*) L)->internal_insert(b->data(i), ((InternalBucket*) R)->child(0));
          ((InternalBucket*) b)->internal_erase_pos(i--);
          delete_leaf((InternalBucket*) R);
        }
      }
    }
    return true;
  }

  bool optimize(Bucket *b = NULL) {
    // fprintf(stderr, "optimize\n");
    if (!b) b = root;
    if (b->is_leaf()) {
      // fprintf(stderr, "optimize leaf\n");
      if (split_chain((LeafBucket*) b)) {
        while (split_chain((LeafBucket*) b));
        return true;
      }
      ((LeafBucket*) b)->leaf_optimize();
      return false;
    }
    // fprintf(stderr, "optimize internal\n");
    bool retry = false;
    for (int i = 0; i <= b->size(); i++) {
      if (optimize(((InternalBucket*) b)->child(i)))
        i = -1, retry = true;
      // fprintf(stderr, "optimize internal i = %d / %d\n", i, b->size());
    }
    return retry;
  }

  int max_depth(Bucket *b = NULL) {
    if (!b) b = root;
    if (b->is_leaf()) return 1;
    int ret = 1;
    for (int i = 0; i <= b->size(); i++) {
      ret = max(ret, 1 + max_depth(((InternalBucket*) b)->child(i)));
    }
    return ret;
  }

  int slack(Bucket *b = NULL) {
    if (!b) b = root;
    int ret = b->get_cap() - b->size();
    if (b->is_leaf()) return ret;
    for (int i = 0; i <= b->size(); i++) {
      ret += slack(((InternalBucket*) b)->child(i));
    }
    return ret;
  }

  pair<Bucket*, int> find_leaf_bucket(int value) {
    Bucket *b = root;
    while (true) {
      if (b->is_leaf()) {
        LeafBucket *Lb = (LeafBucket*) b;
        if (!split_chain(Lb)) break;
        assert(Lb->get_parent());
        b = Lb->get_parent();
      } else {
        InternalBucket *ib = (InternalBucket*) b;
        int pos = ib->internal_lower_pos(value);
        if (pos < ib->size() && ib->data(pos) == value) {
          return make_pair(ib, pos); // Found in the internal bucket.
        }
        b = ib->child(pos);    // Search the child.
      }
    }
    return make_pair(b, 0);
  }

  pair<bool, int> lower_bound(int value) {
    pair<bool, int> ret = make_pair(false, 0);
    pair<Bucket*, int> p;
    // t1 += time_it([&] {
      // fprintf(stderr, "lower_bound %d\n", value);
      p = find_leaf_bucket(value);
    // });
    // t2 += time_it([&] {
      // Found in internal bucket.
      if (!p.first->is_leaf()) {
        ret = make_pair(true, value);
      } else {
          LeafBucket *b = (LeafBucket*) p.first;
          int pos = b->leaf_lower_pos(value);
          if (pos < b->size()) {
            ret = make_pair(true, b->data(pos));
          } else {
              InternalBucket *ib = (InternalBucket*) b->get_parent();
              while (ib) {
                pos = ib->internal_lower_pos(value);
                if (pos < ib->size()) {
                  ret = make_pair(true, ib->data(pos));
                  break;
                }
                ib = (InternalBucket*) ib->get_parent();
              }
          }
      }
    // });
    return ret;
  }

  void insert(int value) {
    // if (t1 > 0) fprintf(stderr, "ins %d\n", value);
    // if (value == 711)  debug();
    Bucket *b = root;
    while (!b->is_leaf()) b = ((InternalBucket*) b)->child_bucket(value);
    ((LeafBucket*) b)->leaf_insert(value);
    // root->debug(0);
  }

  bool erase(int value) {
    // fprintf(stderr, "erase %d, is leaf = %d\n", value, L.first->is_leaf());
    pair<Bucket*, int> p = find_leaf_bucket(value);
    if (p.first->is_leaf()) {
      // fprintf(stderr, "leaf_erase\n");
      return ((LeafBucket*) p.first)->leaf_erase(value); // May return false.
    }
    // fprintf(stderr, "internal_erase\n");
    return ((InternalBucket*) p.first)->internal_erase(value); // Always return true.
  }

  bool check() {
    return root->check(-2147483648);
  }
};
}

#endif
