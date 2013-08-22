#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cassert>
#include <queue>
#include <algorithm>
#include "random.h"

using namespace std;

namespace ctree {

#define BSIZE 128 // Must be divisible by two.
// #define BSIZE 20 // Must be divisible by two.

class Bucket {
 protected:
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N;

 public:
  int size() { return N; }
  bool is_full() { return size() == BSIZE; }
  bool is_leaf() { return pending_insert >= 0; }
  void optimize();
  int debug(int depth);
  bool check();
};

class LeafBucket : public Bucket {
 protected:
  int D[BSIZE];
  LeafBucket *next, *tail;  // Store pending inserts in a linked list.

 public:
  LeafBucket();
  LeafBucket* next_bucket() { return next; }
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }
  void leaf_insert(int v);
  template<typename Func> void leaf_split(Func f);
  void leaf_optimize();
  int promote_last();
  int lower_pos(int value);
  pair<bool,int> leaf_lower_bound(int value);
  void set_standalone();
  void add_chain(LeafBucket *b);
};

class InternalBucket : public LeafBucket {
 protected:
  Bucket *C[BSIZE + 1];

 public:
  InternalBucket();
  InternalBucket(int promotedValue, Bucket *left, Bucket *right);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  int child_pos(int value);
  Bucket*& child_bucket(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
};





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
  assert(N >=0 && N <= BSIZE);
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


LeafBucket::LeafBucket() {
  N = 0;
  pending_insert = 0;
  next = tail = NULL;
}

void LeafBucket::leaf_insert(int value) {
  assert(is_leaf());
  assert(N >= 0);
  pending_insert++;
  if (!is_full()) {
    D[N++] = value;
  } else {
    if (!tail) next = tail = new LeafBucket();
    if (tail->is_full()) add_chain(new LeafBucket());
    tail->D[tail->N++] = value;
    tail->pending_insert++;
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

void LeafBucket::set_standalone() {
  next = tail = NULL;
  pending_insert = 1;
}

void LeafBucket::add_chain(LeafBucket *b) {
  if (!next || !tail) next = tail = b;
  else {
    tail->next = b;
    tail = b;
  }
}

void _add(LeafBucket *&b, LeafBucket *nb) {
  nb->set_standalone();
  if (!b) b = nb;
  else b->add_chain(nb);
}

template<typename Func>
void LeafBucket::leaf_split(Func f) {
  if (!next) return;

  // fprintf(stderr, "split N = %d\n", N);
  // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
  int R[11];
  Random rng(140384); // TODO: use randomized seed.
  for (int i = 0; i < 11; i++) {
    int j = rng.nextInt(N);
    R[i] = D[j];
    D[j] = D[--N];
  }
  assert(N >= 0);
  // fprintf(stderr, "split2 N = %d\n", N);

  LeafBucket *Lb, *Rb;
  queue<LeafBucket*> q;
  q.push(Lb = this);

  // Replace elements with gradually decreasing probability.
  for (int i = 1; Lb->next; i++) {
    q.push(Lb = Lb->next);
    int j = rng.nextInt(i);
    if (j < 11) {
      int k = rng.nextInt(Lb->N);
      // fprintf(stderr, "swap %d  <>  %d,   %d %d\n", R[j], Lb->D[k], j, k);
      swap(R[j], Lb->D[k]);
    }
  }

  for (int i = 0; i < 11; i++) {
    // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
  }

  std::nth_element(R, R + 5, R + 11);
  int pivot = R[5];
  R[5] = R[10];
  for (int i = 0; i < 10; i++) {
    D[N++] = R[i];
  }
  D[N++] = Lb->D[--Lb->N];
  // fprintf(stderr, "split3 N = %d\n", N);
  // fprintf(stderr, "split3 N = %d, pivot = %d\n", next->N, pivot);

  // debug(10);

  LeafBucket *chain[2] { NULL, NULL };
  Lb = Rb = NULL;
  int hi[BSIZE], lo[BSIZE];
  int nhi = 0, nlo = 0;
  while (true) {
    if (nhi && nlo) {
      assert(Lb && Rb);
      // fprintf(stderr, "fusion\n");
      fusion(Lb->D, Rb->D, hi, lo, nhi, nlo);
      // fprintf(stderr, "fusiond %d %d, %p %p, %p %p\n", nhi, nlo, Lb, Rb, chain[0], chain[1]);
      if (!nhi) { _add(chain[0], Lb); Lb = NULL; }
      // fprintf(stderr, "fusiondd1\n");
      if (!nlo) { _add(chain[1], Rb); Rb = NULL; }
      // fprintf(stderr, "fusiondd2\n");
    } else if (!Lb) {
      if (q.empty()) break;
      Lb = q.front(); q.pop();
      Lb->set_standalone();
      if (Lb->size() < BSIZE) break;
    } else if (!nhi) {
      assert(Lb);
      mark_hi(Lb->D, BSIZE, pivot, hi, nhi);
      if (!nhi){ _add(chain[0], Lb); Lb = NULL; }
    } else if (!Rb) {
      if (q.empty()) break;
      Rb = q.front(); q.pop();
      Rb->set_standalone();
      if (Rb->size() < BSIZE) break;
    } else if (!nlo) {
      assert(Rb);
      mark_lo(Rb->D, BSIZE, pivot, lo, nlo);
      if (!nlo){ _add(chain[1], Rb); Rb = NULL; }
    } else {
      assert(0);
    }
  }
  assert(q.empty());

  // fprintf(stderr, "splited\n");

  if (!chain[0]) {
    // fprintf(stderr, "no left\n");
    assert(Lb == this);
    chain[0] = Lb;
    for (int i = 0; i < Lb->N; i++) {
      // fprintf(stderr, "hihi %d < %d\n", i, Lb->N);
      if (Lb->D[i] >= pivot) {
        if (!chain[1]) chain[1] = new LeafBucket();
        chain[1]->leaf_insert(Lb->D[i]);
        Lb->D[i--] = Lb->D[--Lb->N];
      }
    }
  } else if (Lb) {
    // fprintf(stderr, "has left\n");
    assert(Lb != this);
    while (Lb->N) {
      int j = pivot <= Lb->D[--Lb->N];
      // fprintf(stderr, "proc left %d, j = %d\n", Lb->N, j);
      if (!chain[j]) chain[j] = new LeafBucket();
      chain[j]->leaf_insert(Lb->D[Lb->N]);
    }
    delete Lb;
  }

  if (Rb) {
    // fprintf(stderr, "no right\n");
    assert(chain[0]);
    while (Rb->N) {
      int j = pivot <= Rb->D[--Rb->N];
      if (!chain[j]) chain[j] = new LeafBucket();
      chain[j]->leaf_insert(Rb->D[Rb->N]);
    }
    delete Rb;
  }

  f(pivot, chain[1]);
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
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  // std::lower_bound(D, D + N, value) - D;
  return pos;
}

pair<bool,int> LeafBucket::leaf_lower_bound(int value) {
  if (pending_insert > 0) {
    sort(D, D + N);
    pending_insert = 0;
  }
  int pos = lower_pos(value);
  return (pos == N) ? make_pair(false, 0) : make_pair(true, D[pos]);
}

InternalBucket::InternalBucket() {
  pending_insert = -1;
}

InternalBucket::InternalBucket(int promotedValue, Bucket *left, Bucket *right) {
  pending_insert = -1;
  D[0] = promotedValue;
  C[0] = left;
  C[1] = right;
  N = 1;
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

  CTree() {
    root = new LeafBucket();
  }

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

  pair<bool, int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    int i = 0;
    Bucket *b = root;
    pair<InternalBucket*, int> parents[10];
    // if (!root->is_leaf()) parents[i++] = make_pair((InternalBucket*) root, ((InternalBucket*) root)->child_pos(value));
      // fprintf(stderr, "leaf_split iii = %d\n", i);
    while (!b->is_leaf() || ((LeafBucket*) b)->next_bucket()) {
      // fprintf(stderr, "leaf_split ii = %d\n", i);
      for (; !b->is_leaf(); i++) {
        assert(i < 10);
        parents[i] = make_pair((InternalBucket*) b, ((InternalBucket*) b)->lower_pos(value));
        // fprintf(stderr, "child pos %d <= %d\n", parents[i].second, parents[i].first->size());
        if (parents[i].second < parents[i].first->size() && parents[i].first->data(parents[i].second) == value)
          return make_pair(true, value);
        b = parents[i].first->child(parents[i].second);
      }

      // fprintf(stderr, "leaf_split i = %d\n", i);
      ((InternalBucket*) b)->leaf_split([&] (int promotedValue, LeafBucket *nb) { // Split if the bucket has chain.
        while (true) {
          if (i == 0) {
            // fprintf(stderr, "NEW ROOT\n");
           b = root = new InternalBucket(promotedValue, b, nb); break; }

          i--;

          b = parents[i].first;
          if (!b->is_full()) {
            // fprintf(stderr, "PARENT INSERT\n");
            ((InternalBucket*) b)->internal_insert(promotedValue, nb); break; }

          InternalBucket *inb = ((InternalBucket*) b)->internal_split();
          int promotedValueInternal = ((InternalBucket*) b)->internal_promote_last();
          if (promotedValue >= promotedValueInternal) {
            inb->internal_insert(promotedValue, nb);
          } else {
            ((InternalBucket*) b)->internal_insert(promotedValue, nb);
          }
          promotedValue = promotedValueInternal;
          nb = inb;
        }
      });

      // debug();
    }

    auto ret = ((LeafBucket*) b)->leaf_lower_bound(value);
    while (!ret.first && i > 0) {
      auto &p = parents[--i];
      if (p.second < p.first->size()) {
        // fprintf(stderr, "here %d\n", p.first->data(p.second));
        return make_pair(true, p.first->data(p.second));
      }
    }
    // fprintf(stderr, "awww %d %d\n", ret.first, ret.second);
      // debug();
    return ret;
  }
};

}

#endif
