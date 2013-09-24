#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree.h"
#include "test.h"

using namespace std;
using namespace ctree;

CTree<long long> c;

void init(long long *arr, int N) {
  // c.load("ctree"); return;

  c.batch_insert(arr, N);
  // fprintf(stderr, "doi \n" );
  // for (int i = 0; i < N; i++) {
  //   c.insert(arr[i]);
  // }
  // assert(c.check());
  #ifdef EAGER
    c.optimize();
    locked = 1;
  #endif
  // c.save("ctree");
  // // c.debug();
}

void insert(long long value) {
  c.insert(value);
}

void erase(long long value) {
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

long long query(long long value) {
  auto it = c.lower_bound(value);
  long long ret = it.first ? it.second : 0;
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
