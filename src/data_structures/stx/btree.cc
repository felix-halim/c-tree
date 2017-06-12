#include <cstdio>
#include <cassert>
#include "stx/btree_multiset"
#include "../../tester.h"

stx::btree_multiset<long long> s;

void initexp() {}
void destroyexp() {}

// op = 1: inserts the value.
void insert(long long value) {
  s.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  auto it = s.lower_bound(value);
  if (it != s.end()) {
    if (*it == value) {
      s.erase(it);
      return true;
    }
  }
  return false;
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  auto it1 = s.lower_bound(a);
  auto it2 = s.lower_bound(b);
  long long ret1 = (it1 == s.end()) ? 0 : *it1;
  long long ret2 = (it2 == s.end()) ? 0 : *it2;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  auto it1 = s.lower_bound(a);
  // auto it2 = s.lower_bound(b);
  long long sum = 0;
  while (it1 != s.end()) {
    long long v = *(it1++);
    if (v >= b) break;
    // fprintf(stderr, "  %lld\n", *it1);
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
