#include <cstdio>
#include <cassert>
#include "comb.h"
#include "test.h"

// Comb<int, std::less<int>, true, 3200, 125, 50> c;
Comb<long long> c;

void init(long long *arr, int N) {
  for (int i = 0; i < N; i++)
    c.insert(arr[i]);
}

void insert(long long value) {
  c.insert(value);
}

void erase(long long value) {
  // fprintf(stderr, "erase %d\n", value);
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

long long query(long long value) {
  long long val = 0;
  long long ret = c.lower_bound(value).next(val) ? val : 0;
  fprintf(stdout, "%lld (%lld)\n", ret, value);
  return ret;
}

void results(Statistics &s) {
  // assert(c.check());
  s.note = "Lazy";
  s.n_leaves = c.num_of_buckets();
  s.n_capacity = c.num_of_buckets() * c.bucket_size();
  s.n_internals = 1;
  s.max_depth = 2;
  s.slack = c.slack();
  s.in_size = c.root_size();
  s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
