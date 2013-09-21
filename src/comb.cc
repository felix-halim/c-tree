#include <cstdio>
#include <cassert>
#include "comb.h"
#include "test.h"

// Comb<int, std::less<int>, true, 3200, 125, 50> c;
Comb<int> c;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    c.insert(arr[i]);
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  // fprintf(stderr, "erase %d\n", value);
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

int query(int value) {
  int val = 0;
  int ret = c.lower_bound(value).next(val) ? val : 0;
  fprintf(stdout, "%d (%d)\n", ret, value);
  return ret;
}

void results(Statistics &s) {
  assert(c.check());
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
