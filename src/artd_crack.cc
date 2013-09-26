#include <cstdio>
#include <cassert>
#include <algorithm>
#include "artd_crack.h"
#include "test.h"

using namespace std;
using namespace art_crack;

ArtCrack<long long> c;

inline uintptr_t getLeafValue(Node* node) {
  // The the value stored in the pseudo-leaf
  return reinterpret_cast<uintptr_t>(node) >> 1;
}

void init(long long *arr, int N) {
  c.batch_insert(arr, N);
  // art_debug = 1;
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
  // art_debug = 1;
  // if (value == 277929528LL) art_debug = 1;
  auto it = c.lower_bound(value);
  long long ret = it.first ? it.second : 0;
  // fprintf(stdout, "%lld (%lld)\n", ret, value); fflush(stdout);
  // if (value == 277929528LL) fprintf(stderr, "%d %lld\n", it.first, it.second), exit(1);
  // art_debug = 0;
  return ret;
}

void results(Statistics &s) {
  s.n_leaves = nLeaves;
  s.n_capacity = nCap;
  s.n_internals = nInternals;
  s.ln_size = LEAF_BSIZE;
  c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
