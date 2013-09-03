#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree_exp_leafsize.h"
#include "test.h"

using namespace std;
using namespace ctree;

CTree c;

void init(int *arr, int N) {
  // c.load("ctree"); return;

  c.batch_insert(arr, N);
  // fprintf(stderr, "doi \n" );
  // for (int i = 0; i < N; i++) {
  //   c.insert(arr[i]);
  // }
  // assert(c.check());
  // c.optimize();
  // c.save("ctree");
  // locked = 1;
  // // c.debug();
  // fprintf(stderr, "depth = %d, slack = %d\n", c.max_depth(), c.slack());
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  bool ok = c.erase(value);
  assert(ok);
}

int query(int value) {
  auto it = c.lower_bound(value);
  int ret = it.first ? it.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

void results(Statistics &s) {
  assert(c.check());
  s.note = "Exp Size";
  s.n_leaves = nLeaves;
  s.n_capacity = nCap;
  s.n_internals = nInternals;
  s.max_depth = c.max_depth();
  s.slack = c.slack();
  s.in_size = INTERNAL_BSIZE;
  s.ln_size = LEAF_BSIZE;
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
