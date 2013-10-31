#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree.h"
#include "test.h"

using namespace std;
using namespace ctree;

CTree<int> c;

void init(int *arr, int N) {
  // c.load("ctree"); return;

  c.batch_insert(arr, N);
  // fprintf(stderr, "doi \n" );
  // for (int i = 0; i < N; i++) {
  //   c.insert(arr[i]);
  // }
  // assert(c.check());
  #ifdef EAGER
    c.optimize();
    // locked = 1;
  #endif
  // c.save("ctree");
  // // c.debug();
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

int lower_bound(int value) {
  auto it = c.lower_bound(value);
  int ret = it.first ? it.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

int select(int a, int b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  int ret1 = it1.first ? it1.second : 0;
  int ret2 = it2.first ? it2.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

int sum(int a, int b) {
  // auto it1 = c.lower_bound(a);
  // auto it2 = c.lower_bound(b);
  int ret = 0;
  // for (int val; it1 != it2 && it1.next(val); ret += val);
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

void results(Statistics &s) {
  // assert(c.check());
  #ifdef EAGER
    s.note = "Eager";
  #else
    s.note = "Lazy";
  #endif
  s.n_leaves = nLeaves;
  s.n_capacity = nCap;
  s.n_internals = nInternals;
  s.max_depth = c.max_depth();
  s.slack = c.slack();
  s.in_size = INTERNAL_BSIZE;
  s.ln_size = LEAF_BSIZE;
  c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
