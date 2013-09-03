#include <cstdio>
#include <cassert>
#include "test.h"

int *arr, N;
double t;

void init(int *iarr, int iN) {
  t = time_it([&] {
    N = iN;
    arr = new int[N];
    for (int i = 0; i < N; i++)
      arr[i] = iarr[i];
  });
  sort(arr, arr + N);
}

void insert(int value) {}
void erase(int value) {}

int query(int value) {
  auto it = lower_bound(arr, arr + N, value);
  return (it == arr + N) ? 0 : *it;
}

void results(Statistics &s) {
  // assert(c.check());
  char runtime[100];
  sprintf(runtime, "%.6lf", t);
  s.note = string("insert = ") + runtime;
  // s.n_leaves = c.num_of_buckets();
  // s.n_capacity = c.num_of_buckets() * c.bucket_size();
  // s.n_internals = 1;
  // s.max_depth = 2;
  // s.slack = c.slack();
  // s.in_size = c.root_size();
  // s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
