#include <cstdio>
#include <cassert>
#include "stx/btree_multiset"
#include "test.h"

stx::btree_multiset<int> b;

void init(int *arr, int N) {
  b.insert(arr, arr + N);
}

void insert(int value) {
  b.insert(value);
}

void erase(int value) {
  auto it = b.lower_bound(value);
  assert(it != b.end());
  b.erase(it);
}

int lower_bound(int value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

int select(int lo, int hi) {
  auto it1 = b.lower_bound(lo);
  auto it2 = b.lower_bound(hi);
  int ret1 = (it1 == b.end()) ? 0 : *it1;
  int ret2 = (it2 == b.end()) ? 0 : *it2;
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
