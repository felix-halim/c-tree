#include <cstdio>
#include <cassert>
#include "google/btree_set.h"
#include "test.h"

btree::btree_multiset<int> b;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    b.insert(arr[i]);
}

void insert(int value) {
  b.insert(value);
}

void erase(int value) {
  auto it = b.lower_bound(value);
  assert(it != b.end());
  b.erase(it);
}

int query(int value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

void results(Statistics &s) {
  // assert(c.check());
  s.note = "Google";
  // s.n_leaves = c.num_of_buckets();
  // s.n_capacity = c.num_of_buckets() * c.bucket_size();
  // s.n_internals = 1;
  // s.max_depth = 2;
  // s.slack = c.slack();
  // s.in_size = c.root_size();
  // s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
