#include <cstdio>
#include <cassert>
#include "stx/btree_multiset"
#include "test.h"

stx::btree_multiset<unsigned> b;

void init(unsigned *arr, unsigned N) {
  b.insert(arr, arr + N);
}

void insert(unsigned value) {
  b.insert(value);
}

void erase(unsigned value) {
  auto it = b.lower_bound(value);
  assert(it != b.end());
  b.erase(it);
}

unsigned lower_bound(unsigned value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

unsigned select(unsigned lo, unsigned hi) {
  auto it1 = b.lower_bound(lo);
  auto it2 = b.lower_bound(hi);
  unsigned ret1 = (it1 == b.end()) ? 0 : *it1;
  unsigned ret2 = (it2 == b.end()) ? 0 : *it2;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

void results(Statistics &s) {
  // assert(c.check());
  s.N = b.size();
  // s.n_leaves = c.num_of_buckets();
  // s.n_capacity = c.num_of_buckets() * c.bucket_size();
  // s.n_internals = 1;
  // s.max_depth = 2;
  // s.slack = c.slack();
  // s.in_size = c.root_size();
  // s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
