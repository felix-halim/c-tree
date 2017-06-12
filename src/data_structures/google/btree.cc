#include <cstdio>
#include <cassert>
#include "google/btree_set.h"
#include "../../tester.h"

btree::btree_multiset<long long> c;

void initexp() {}
void destroyexp() {}

// op = 1: inserts the value.
void insert(long long value) {
  c.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  return c.erase(value);
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  auto it = c.lower_bound(a);
  return (it == c.end()) ? 0 : *it;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  long long sum = 0;
  for (auto it = c.lower_bound(a); it != c.end(); it++) {
    long long v = *it;
    if (v >= b) break;
    sum += v;
  }
  return sum;
}

  // assert(c.check());
  // s.n_leaves = c.num_of_buckets();
  // s.n_capacity = c.num_of_buckets() * c.bucket_size();
  // s.n_internals = 1;
  // s.max_depth = 2;
  // s.slack = c.slack();
  // s.in_size = c.root_size();
  // s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
