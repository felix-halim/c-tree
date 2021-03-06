#include <cstdio>
#include <cassert>
#include <algorithm>
#include "art_crack.h"
#include "test.h"

using namespace std;
using namespace art_crack;

ArtCrack<int> c;

inline uintptr_t getLeafValue(Node* node) {
  // The the value stored in the pseudo-leaf
  assert(isLeaf(node));
  return c.bucket_head_value(reinterpret_cast<uintptr_t>(node)>>1);
}

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
  // art_debug = 1;
  if (value == 1904630696647192210LL) art_debug = 1;
  auto it = c.lower_bound(value);
  int ret = it.first ? it.second : 0;
  // fprintf(stdout, "%lld (%lld)\n", ret, value);
  if (value == 1904630696647192210LL) fprintf(stderr, "%d %lld\n", it.first, it.second);
  art_debug = 0;
  return ret;
}

void results(Statistics &s) {
  s.n_leaves = nLeaves;
  s.n_capacity = nCap;
  s.n_internals = nInternals;
  s.ln_size = LEAF_BSIZE;
  c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
