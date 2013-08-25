// TODO: why transitioning from exponential buckets to 64 sized buckets is so costly?
// If this problem is solved, then a cleaner solution can be achieved. Crack is not needed, use btree internals.
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

#define INTERNAL_BSIZE 60   // Must be power of two.
#define MAX_LEAF_BSIZE 2048     // Must be power of two.

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

int nLeaves, nInternals, nCap, nDes;
Random rng;

class Bucket {
 protected:
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N, cap;

 public:
  int get_cap() { return cap; }
  int size() { return N; }
  bool is_full() { return size() == cap; }
  bool is_leaf() { return pending_insert >= 0; }
  void optimize();
  int debug(int depth);
  bool check();
};

class LeafBucket : public Bucket {
 protected:
  int *D;
  LeafBucket *next, *tail;  // Store pending inserts in a linked list.

 public:
  ~LeafBucket();
  LeafBucket(int cap);
  void init(int cap);
  LeafBucket* next_bucket() { return next; }
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }
  void leaf_insert(int v);
  void leaf_split(vector<pair<int, LeafBucket*>> &nbs);
  void leaf_optimize();
  int promote_last();
  int lower_pos(int value);
  pair<bool,int> leaf_lower_bound(int value);
  LeafBucket* detach_and_get_next();
  void add_chain(LeafBucket *b);
  void distribute_values(int pivot, LeafBucket* chain[2]);
  LeafBucket* transfer_to(LeafBucket *b, int pivot);
};

class InternalBucket : public LeafBucket {
 protected:
  Bucket *C[INTERNAL_BSIZE + 1];

 public:
  ~InternalBucket();
  InternalBucket();
  InternalBucket(int promotedValue, Bucket *left, Bucket *right);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  int child_pos(int value);
  Bucket*& child_bucket(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
};



vector<LeafBucket*> free_leaves[30];

LeafBucket* new_leaf(int cap) {
  return new LeafBucket(cap);

  for (int i = 2; ; i++) {
    if ((1 << i) == cap) {
      if (free_leaves[i].empty()) {
        free_leaves[i].push_back(new LeafBucket(cap));
      }
      LeafBucket* leaf = free_leaves[i].back();
      leaf->init(cap);
      free_leaves[i].pop_back();
      return leaf;
    }
  }
}

void delete_leaf(LeafBucket *b) {
  delete b;
  return;
  for (int i = 2; ; i++) {
    if ((1 << i) == b->get_cap()) {
      free_leaves[i].push_back(b);
      break;
    }
  }
}

void Bucket::optimize() {
  if (is_leaf()) {
    ((LeafBucket*) this)->leaf_optimize();
  } else {
    for (int i = 0; i <= N; i++) {
      ((InternalBucket*) this)->child(i)->optimize();
    }
  }
}

int Bucket::debug(int depth) {
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  LeafBucket *b = (LeafBucket *) this;
  int sz = 0;
  while (b) {
    sz += b->size();
    fprintf(stderr, "N = %d (%d%s), ", b->size(), b->pending_insert, b->is_leaf() ? ", leaf" : "");
    for (int i = 0; i < b->size(); i++) {
      fprintf(stderr, "%d ", b->data(i));
    }
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

bool Bucket::check() {
  assert(N >=0 && N <= cap);
  assert(is_leaf() || ((InternalBucket*) this)->child(N));
  for (int i = 0; i < N; i++) {
    if (is_leaf()) {
    } else if (i > 0) {
      assert(((InternalBucket*) this)->child(i));
      assert(((InternalBucket*) this)->data(i - 1) < ((InternalBucket*) this)->data(i));
    }
  }
  return true;
}


LeafBucket::LeafBucket(int cap) {
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

void LeafBucket::leaf_insert(int value) {
  assert(is_leaf());
  assert(N >= 0);
  pending_insert = 1;
  if (!is_full()) {
    D[N++] = value;
  } else {
    if (!tail) {
      assert(cap == INTERNAL_BSIZE);
      add_chain(new_leaf(INTERNAL_BSIZE));
    } else if (tail->is_full()) {
      add_chain(new_leaf(std::min(tail->cap * 2, MAX_LEAF_BSIZE)));
    }
    tail->D[tail->N++] = value;
    tail->pending_insert = 1;
  }
  // assert(check());
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

void LeafBucket::leaf_split(vector<pair<int, LeafBucket*>> &ret) {
  ret.clear();
  assert(next);
  assert(cap == INTERNAL_BSIZE); // The first bucket must be the smallest capacity.

  if (!next->next && false) {
    LeafBucket *b = detach_and_get_next(); b->detach_and_get_next();
    int T[INTERNAL_BSIZE];
    for (int i = 0; i < N; i++)
      T[i] = D[i];

    int nT = N, nB = b->N, i = 0, j = 0;
    N = b->N = 0;
    sort(T, T + nT);
    sort(b->D, b->D + nB);
    while (N < INTERNAL_BSIZE && i < nT && j < nB) {
      if (T[i] <= b->D[j]) {
        D[N++] = T[i++];
      } else {
        D[N++] = b->D[j++];
      }
    }
    while (N < INTERNAL_BSIZE && i < nT) D[N++] = T[i++];
    while (N < INTERNAL_BSIZE && j < nB) D[N++] = b->D[j++];
    pending_insert = 0;

    LeafBucket *pnb = this;
    while (i < nT || j < nB) {
      LeafBucket *nb = new_leaf(INTERNAL_BSIZE);
      while (nb->N < INTERNAL_BSIZE && i < nT && j < nB) {
        if (T[i] <= b->D[j]) {
          nb->D[nb->N++] = T[i++];
        } else {
          nb->D[nb->N++] = b->D[j++];
        }
      }
      while (nb->N < INTERNAL_BSIZE && j < nB) nb->D[nb->N++] = b->D[j++];
      while (nb->N < INTERNAL_BSIZE && i < nT) nb->D[nb->N++] = T[i++];
      nb->pending_insert = 0;
      ret.push_back(make_pair(pnb->D[--pnb->N], nb));
      pnb = nb;
    }

    delete_leaf(b);

  } else if (!next->next->next && false) {
    LeafBucket *b1 = detach_and_get_next();
    LeafBucket *b2 = b1->detach_and_get_next(); b2->detach_and_get_next();

    int *T = new int[N + b1->N + b2->N];
    int nT = 0;
    for (int i = 0; i < N; i++) T[nT++] = D[i];
    for (int i = 0; i < b1->N; i++) T[nT++] = b1->D[i];
    for (int i = 0; i < b2->N; i++) T[nT++] = b2->D[i];

    sort(T, T + nT);
    int i = 0;
    while (i < nT && i < INTERNAL_BSIZE) {
      D[i] = T[i];
      i++;
    }
    N = i;
    pending_insert = 0;

    LeafBucket *pnb = this;
    while (i < nT) {
      LeafBucket *nb = new_leaf(INTERNAL_BSIZE);
      nb->pending_insert = 0;
      while (i < nT && !nb->is_full()) nb->D[nb->N++] = T[i++];
      ret.push_back(make_pair(pnb->D[--pnb->N], nb));
      pnb = nb;
    }

    delete_leaf(b1);
    delete_leaf(b2);

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

    LeafBucket *chain[2] { this, new_leaf(INTERNAL_BSIZE) };

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
    ret.push_back(make_pair(pivot, chain[1]));
  }
}

int LeafBucket::promote_last() {
  nth_element(D, D + N - 1, D + N);
  return D[--N];
}

void LeafBucket::leaf_optimize() {
  assert(pending_insert >= 0);
  sort(D, D + N);
  pending_insert = 0;
}

int LeafBucket::lower_pos(int value) {
  // return std::lower_bound(D, D + N, value) - D;
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

pair<bool,int> LeafBucket::leaf_lower_bound(int value) {
  if (pending_insert) leaf_optimize();
  int pos = lower_pos(value);
  return (pos == N) ? make_pair(false, 0) : make_pair(true, D[pos]);
}

InternalBucket::InternalBucket() : LeafBucket(INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
}

InternalBucket::InternalBucket(int promotedValue, Bucket *left, Bucket *right) : LeafBucket(INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
  D[0] = promotedValue;
  C[0] = left;
  C[1] = right;
  N = 1;
}

InternalBucket::~InternalBucket() {
  nInternals--;
}

InternalBucket* InternalBucket::internal_split() {
  InternalBucket *nb = new InternalBucket();
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
    nb->C[j] = C[i + 1];
  }
  N /= 2;
  nb->N = N;
  nb->C[0] = C[N];
  // assert(check());
  return nb;
}


void InternalBucket::internal_insert(int value, Bucket *b) {
  // assert(check());
  assert(!is_full());
  int i = N - 1;
  assert(i >= 0);
  while (i >= 0 && D[i] > value) {
    D[i + 1] = D[i];
    C[i + 2] = C[i + 1];
    i--;
  }
  D[i + 1] = value;
  C[i + 2] = b;
  N++;
  // assert(check());
}

int InternalBucket::internal_promote_last() {
 return D[--N];
}

Bucket*& InternalBucket::child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}


class CTree {
  Bucket *root;

 public:

  const char *version = "exp leaf size";

  CTree() {
    root = new_leaf(INTERNAL_BSIZE);
  }

  void optimize() { root->optimize(); }
  void debug() {
    root->debug(0);
    // fprintf(stderr, "\n");
  }

  void insert(int value) {
    // fprintf(stderr, "ins %d\n", value);
    Bucket *b = root;
    while (!b->is_leaf()) b = ((InternalBucket*) b)->child_bucket(value);
    ((LeafBucket*) b)->leaf_insert(value);
    // root->debug(0);
  }

  bool erase(int value) {
    return true;
  }

  pair<InternalBucket*, int> parents[10];
  vector<pair<int, LeafBucket*>> nbs;
  double t1 = 0, t2 = 0, t3 = 0;

  pair<bool, int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    int i = 0;
    Bucket *b = root;
    // if (!root->is_leaf()) parents[i++] = make_pair((InternalBucket*) root, ((InternalBucket*) root)->child_pos(value));
      // fprintf(stderr, "leaf_split iii = %d\n", i);
    while (true) {
      // fprintf(stderr, "leaf_split ii = %d\n", i);
      if (!b->is_leaf()) {
        // t1 += time_it([&] {
          assert(i < 10);
          parents[i].first = (InternalBucket*) b;
          parents[i].second = ((InternalBucket*) b)->lower_pos(value);
          // fprintf(stderr, "child pos %d <= %d\n", parents[i].second, parents[i].first->size());
          b = parents[i].first->child(parents[i].second);
          i++;
        // });

        if (parents[i-1].second < parents[i-1].first->size() && parents[i-1].first->data(parents[i-1].second) == value)
          return make_pair(true, value);

      } else if (((LeafBucket*) b)->next_bucket()) {
        // t2 += time_it([&] {
          ((InternalBucket*) b)->leaf_split(nbs);
          while (!nbs.empty()) {
            if (i == 0) {
              b = root = new InternalBucket(nbs[0].first, b, nbs[0].second);
              for (int j = 1; j < (int) nbs.size(); j++) {
                ((InternalBucket*) b)->internal_insert(nbs[j].first, nbs[j].second);
              }
              break;
            }
            b = parents[--i].first;
            while (!b->is_full() && !nbs.empty()) {
              ((InternalBucket*) b)->internal_insert(nbs.back().first, nbs.back().second);
              nbs.pop_back();
            }
            if (nbs.empty()) break;
            InternalBucket *inb = ((InternalBucket*) b)->internal_split();
            int promotedValueInternal = ((InternalBucket*) b)->internal_promote_last();
            for (int j = 0; j < (int) nbs.size(); j++) {
              if (nbs[j].first >= promotedValueInternal) {
                assert(!inb->is_full());
                inb->internal_insert(nbs[j].first, nbs[j].second);
              } else {
                assert(!b->is_full());
                ((InternalBucket*) b)->internal_insert(nbs[j].first, nbs[j].second);
              }
            }
            nbs.clear();
            nbs.push_back(make_pair(promotedValueInternal, inb));
          }
        // });
      } else {
        break;
      }
      // debug();
    }

    assert(b->get_cap() == INTERNAL_BSIZE);
    pair<bool, int> ret;
    // t3 += time_it([&] {
      ret = ((LeafBucket*) b)->leaf_lower_bound(value);
      while (!ret.first && i > 0) {
        auto &p = parents[--i];
        if (p.second < p.first->size()) {
          // fprintf(stderr, "here %d\n", p.first->data(p.second));
          ret = make_pair(true, p.first->data(p.second));
          break;
        }
      }
    // });
    // fprintf(stderr, "awww %d %d\n", ret.first, ret.second);
      // debug();
    return ret;
  }
};

}

#endif
