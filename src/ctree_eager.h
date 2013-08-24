#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cassert>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace chrono;

int nLeaves, nCap, nInternals;

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

namespace ctree {

#define BSIZE 60 // Must be divisible by two.
// #define BSIZE 4 // Must be divisible by two.

class Bucket {
 protected:
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N;

 public:
  int size() { return N; }
  bool is_full() { return size() == BSIZE; }
  bool is_leaf() { return pending_insert >= 0; }
  void optimize();
  void debug(int depth);
  bool check();
};

class LeafBucket : public Bucket {
 protected:
  int D[BSIZE];
  // Bucket *chain;  // If set, then t

 public:
  LeafBucket();
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }
  void leaf_insert(int v);
  LeafBucket* leaf_split();
  void leaf_optimize();
  int promote_last();
  pair<bool,int> leaf_lower_bound(int value);
};

class InternalBucket : public LeafBucket {
 protected:
  Bucket *C[BSIZE + 1];

 public:
  InternalBucket();
  InternalBucket(int promotedValue, Bucket *left, Bucket *right);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  Bucket*& upper_child(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
  pair<bool,int> internal_lower_bound(int value);
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

void Bucket::debug(int depth) {
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  fprintf(stderr, "N = %d (%d%s), ", ((LeafBucket*) this)->size(), pending_insert, is_leaf() ? ", leaf" : "");
  for (int i = 0; i < ((LeafBucket*) this)->size(); i++) {
    fprintf(stderr, "[%d] %d, ", i, ((LeafBucket*) this)->data(i));
  }
  fprintf(stderr, "\n");
  if (!is_leaf()) {
    for (int i = 0; i <= ((LeafBucket*) this)->size(); i++) {
      ((InternalBucket*) this)->child(i)->debug(depth + 1);
    }
  }
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
  nLeaves++;
  nCap += BSIZE;
}

void LeafBucket::leaf_insert(int value) {
  assert(is_leaf());
  assert(N >= 0 && N < BSIZE);
  pending_insert++;
  D[N++] = value;
  // assert(check());
}

LeafBucket* LeafBucket::leaf_split() {
  LeafBucket *nb = new LeafBucket();
  if (pending_insert > 0) {
    nth_element(D, D + N/2, D + N);
    nb->pending_insert = 1;
  }
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
  }
  N /= 2;
  nb->N = N;
  // assert(check());
  return nb;
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

pair<bool,int> LeafBucket::leaf_lower_bound(int value) {
  if (pending_insert > 0) {
    sort(D, D + N);
    pending_insert = 0;
  }

  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  // std::lower_bound(D, D + N, value) - D;

  return (pos == N) ? make_pair(false, 0) : make_pair(true, D[pos]);
}

InternalBucket::InternalBucket() {
  pending_insert = -1;
  nInternals++;
}

InternalBucket::InternalBucket(int promotedValue, Bucket *left, Bucket *right) {
  nInternals++;
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

Bucket*& InternalBucket::upper_child(int value) {
  return child(upper_bound(D, D + N, value) - D);
}

pair<bool,int> InternalBucket::internal_lower_bound(int value) {
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  auto ret = C[pos]->is_leaf()
    ? ((LeafBucket*) C[pos])->leaf_lower_bound(value)
    : ((InternalBucket*) C[pos])->internal_lower_bound(value);
  if (ret.first) {
    return ret;
  } else if (pos < N) {
    return make_pair(true, D[pos]);
  }
  return ret;
}

class CTree {
  Bucket *root;

 public:
  int depth;
  const char* version = "eager";

  double t1, t2, t3;

  CTree() {
    root = new LeafBucket();
    depth = 0;
  }

  void debug() {
    root->debug(0);
    fprintf(stderr, "\n");
  }

  void insert(int value) {
    // fprintf(stderr, "ins %d\n", value);
    auto it = insert(root, NULL, value, 0);
    assert(!it.second);
    // root->debug(0);
  }

  void optimize() {
    // fprintf(stderr, "opt\n");
    root->optimize();
  }

  pair<int, Bucket*> leaf_insert(LeafBucket *&b, InternalBucket *p, int value, int depth) {
    if (b->is_full()) {
      LeafBucket *nb = b->leaf_split();
      int promotedValue = b->promote_last();
      if (value >= promotedValue) {
        nb->leaf_insert(value);
      } else {
        b->leaf_insert(value);
      }
      if (p) {
        if (p->is_full()) {
          return make_pair(promotedValue, nb);
        } else {
          p->internal_insert(promotedValue, nb);
        }
      } else {
        b = new InternalBucket(promotedValue, b, nb); // New root.
      }
    } else {
      b->leaf_insert(value);
    }
    return make_pair(0, (Bucket*) NULL);
  }

  pair<int, Bucket*> internal_insert(InternalBucket *&b, InternalBucket *p, int value, int depth) {
    auto it = insert(b->upper_child(value), b, value, depth + 1);
    if (it.second) {
      assert(b->is_full());
      InternalBucket *nb = b->internal_split();
      int promotedValue = b->internal_promote_last();
      if (it.first >= promotedValue) {
        nb->internal_insert(it.first, it.second);
      } else {
        b->internal_insert(it.first, it.second);
      }
      if (p) {
        if (p->is_full()) {
          return make_pair(promotedValue, nb);
        } else {
          p->internal_insert(promotedValue, nb);
        }
      } else {
        b = new InternalBucket(promotedValue, b, nb); // New root.
      }
    }
    return make_pair(0, (Bucket*) NULL);
  }

  pair<int, Bucket*> insert(Bucket *&b, InternalBucket *p, int value, int depth) {
    this->depth = max(this->depth, depth);
    if (b->is_leaf()) {
      return leaf_insert((LeafBucket*&) b, p, value, depth);
    } else {
      return internal_insert((InternalBucket*&) b, p, value, depth);
    }
  }

  bool erase(int value) {
    return true;
  }

  pair<bool,int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    pair<bool, int> ret;
    // t3 += time_it([&] {
     ret = root->is_leaf()
      ? ((LeafBucket*) root)->leaf_lower_bound(value)
      : ((InternalBucket*) root)->internal_lower_bound(value);
    // });
    return ret;
  }
};

}

#endif
