#include <cstdio>
#include <cassert>
#include "test.h"

unsigned *arr, N;
double t;

void init(unsigned *iarr, unsigned iN) {
  t = time_it([&] {
    N = iN;
    arr = new unsigned[N];
    for (int i = 0; i < (int) N; i++)
      arr[i] = iarr[i];
  });
  sort(arr, arr + N);
}

void insert(unsigned value) {}
void erase(unsigned value) {}

unsigned lower_bound(unsigned value) {
  auto it = lower_bound(arr, arr + N, value);
  return (it == arr + N) ? 0 : *it;
}

unsigned select(unsigned a, unsigned b) {
  unsigned i1 = lower_bound(arr, arr + N, a) - arr;
  unsigned i2 = lower_bound(arr, arr + N, b) - arr;
  return ((i1 == N) ? 0 : arr[i1]) + ((i2 == N) ? 0 : arr[i2]);
}

void results(Statistics &s) {
  // assert(c.check());
  // s.n_leaves = c.num_of_buckets();
  // s.n_capacity = c.num_of_buckets() * c.bucket_size();
  // s.n_internals = 1;
  // s.max_depth = 2;
  // s.slack = c.slack();
  // s.in_size = c.root_size();
  // s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
