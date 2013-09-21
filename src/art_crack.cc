#include <cstdio>
#include <cassert>
#include <algorithm>
#include "art_crack.h"
#include "test.h"

using namespace std;
using namespace art_crack;

ArtCrack c;

void init(int *arr, int N) {
  c.batch_insert(arr, N);
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

int query(int value) {
  auto it = c.lower_bound(value);
  int ret = it.first ? it.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

void results(Statistics &s) {
  s.n_leaves = nLeaves;
  s.n_capacity = nCap;
  s.n_internals = nInternals;
  s.ln_size = LEAF_BSIZE;
  c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
