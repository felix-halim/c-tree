#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cassert>
#include <algorithm>

using namespace std;

namespace ctree {

#define BSIZE 256 // Must be divisible by two.
// #define BSIZE 4 // Must be divisible by two.

class Bucket {
 public:
  int D[BSIZE];
  Bucket **C;      // Children (for internal nodes only).
  Bucket *chain;  // Chained buckets (for leaf nodes only).
  int P;          // Pending insert (for leaf nodes only).
  int N;

  Bucket() {}
  Bucket(bool leaf) {
    N = P = 0;
    chain = NULL;
    C = NULL;
    if (!leaf) C = new Bucket*[BSIZE + 1];
  }

  bool is_full() { return N == BSIZE; }
  bool is_leaf() { return C == NULL; }

  void insert_leaf(int v) {
    assert(is_leaf());
    assert(N < BSIZE);
    D[N++] = v;
    P++;
    // assert(check());
  }

  void insert_internal(int value, Bucket *b) {
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

  Bucket* split() {
    Bucket *nb = new Bucket(is_leaf());
    if (is_leaf()) {
      if (P) {
        nth_element(D, D + N/2, D + N);
        nb->P = 1;
      }
      for (int i = N/2, j = 0; i < N; i++)
        nb->D[j++] = D[i];
    } else {
      for (int i = N/2, j = 0; i < N; i++) {
        nb->D[j++] = D[i];
        nb->C[j] = C[i + 1];
      }
    }
    N /= 2;
    nb->N = N;
    // assert(check());
    return nb;
  }

  void optimize() {
    if (is_leaf()) {
      assert(P);
      sort(D, D + N);
      P = 0;
    } else {
      for (int i = 0; i <= N; i++) {
        C[i]->optimize();
      }
    }
  }

  void debug(int depth) {
    for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
    fprintf(stderr, "N = %d (%d%s), ", N, P, is_leaf() ? ", leaf" : "");
    for (int i = 0; i < N; i++) {
      fprintf(stderr, "[%d] %d, ", i, D[i]);
    }
    fprintf(stderr, "\n");

    if (!is_leaf()) {
      for (int i = 0; i <= N; i++) {
        C[i]->debug(depth + 1);
      }
    }
  }

  bool check() {
    assert(N >=0 && N <= BSIZE);
    assert(is_leaf() || C[N]);
    for (int i = 0; i < N; i++) {
      if (is_leaf()) {
      } else if (i > 0) {
        assert(C[i]);
        assert(D[i - 1] < D[i]);
      }
    }
    return true;
  }
};

class CTree {
  Bucket *root;

 public:
  int depth;

  CTree() {
    root = new Bucket(true);
    depth = 0;
  }

  void debug() {
    root->debug(0);
    fprintf(stderr, "\n");
  }

  void insert(int value) {
    // fprintf(stderr, "insert %d\n", value);
    // fprintf(stderr, "a");
    auto it = insert(root, NULL, value, 0);
    // fprintf(stderr, "b");
    assert(!it.second);
    // root->debug(0);
  }

  void optimize() {
    root->optimize();
  }

  pair<int, Bucket*> insert(Bucket *&b, Bucket *p, int value, int depth) {
    this->depth = max(this->depth, depth);
    if (b->is_leaf()) {
      if (b->is_full()) {
        Bucket *nb = b->split();
        nth_element(b->D, b->D + b->N-1, b->D + b->N);
        int promotedValue = b->D[--b->N];
        if (value >= promotedValue) {
          nb->insert_leaf(value);
        } else {
          b->insert_leaf(value);
        }
        if (p) {
          if (p->is_full()) {
            // fprintf(stderr, "FULL %d\n", promotedValue);
            return make_pair(promotedValue, nb);
          } else {
            p->insert_internal(promotedValue, nb);
          }
        } else {
          // fprintf(stderr, "xxxxx %d\n", promotedValue);
          Bucket *root = new Bucket(false); // New root.
          root->D[0] = promotedValue;
          root->C[0] = b;
          root->C[1] = nb;
          root->N = 1;
          b = root;
        }
      } else {
        b->insert_leaf(value);
      }
    } else {
      int pos = upper_bound(b->D, b->D + b->N, value) - b->D;
      Bucket *oldB = b->C[pos];
      auto it = insert(b->C[pos], b, value, depth + 1);
      assert(!p || oldB == b->C[pos]);
      if (it.second) {
        // fprintf(stderr, "RECV %d\n", it.first);
        // it.second->debug(10);
        assert(b->is_full());
        Bucket *nb = b->split();
        nb->C[0] = b->C[b->N];
        int promotedValue = b->D[--b->N];
        if (it.first >= promotedValue) {
          nb->insert_internal(it.first, it.second);
        } else {
          b->insert_internal(it.first, it.second);
        }
        if (p) {
          if (p->is_full()) {
            return make_pair(promotedValue, nb);
          } else {
            p->insert_internal(promotedValue, nb);
          }
        } else {
          // fprintf(stderr, "new root %d\n", promotedValue);
          Bucket *root = new Bucket(false); // New root.
          root->D[0] = promotedValue;
          root->C[0] = b;
          root->C[1] = nb;
          root->N = 1;
          b = root;
        }
      }
    }
    return make_pair(0, (Bucket*) NULL);
  }

  bool erase(int value) {
    return true;
  }

  pair<bool,int> lower_bound(int value) {
    auto ret = lower_bound(root, value);
    // root->debug();
    return ret;
  }

  pair<bool,int> lower_bound(Bucket *b, int value) {
    if (b->is_leaf()) {
      if (b->P) {
        sort(b->D, b->D + b->N);
        b->P = 0;
      }
    }
    int pos = std::lower_bound(b->D, b->D + b->N, value) - b->D;
    // fprintf(stderr, "pos = %d, p = %d, is_leaf = %d\n", pos, b->P, b->is_leaf());
    if (b->D[pos] == value) return make_pair(true, value);
    if (b->is_leaf()) {
      return (pos == b->N) ? make_pair(false, 0) : make_pair(true, b->D[pos]);
    } else {
      auto ret = lower_bound(b->C[pos], value);
      if (ret.first) {
        return ret;
      } else if (pos + 1 < b->N) {
        return make_pair(true, b->D[pos + 1]);
      }
      return ret;
    }
  }
};

}

#endif
