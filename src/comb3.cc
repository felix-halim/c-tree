#include <cstdio>
#include <cassert>
#include "comb3.h"
#include "test.h"

Comb<int, less<int>, 2048, 64> c;

void init(int *arr, int N) {
  c.load(arr, N);
  // for (int i = 0; i < N; i++) {
  //   insert(arr[i]);
  //   if (i % 1000) lower_bound(arr[i]);
  // }
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
    if (!ok) fprintf(stderr, "ERASE error %d\n", value);
    assert(ok);
  #endif
}

int lower_bound(int value) {
  auto it = c.lower_bound(value);
  int ret = it.bucket ? *it : 0;
  // fprintf(stdout, "%d (%d)\n", ret, value);
  return ret;
}

int select(int a, int b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  int ret1 = (it1 == c.end()) ? 0 : *it1;
  int ret2 = (it2 == c.end()) ? 0 : *it2;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

int count(int a, int b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  return c.count(it1, it2);
}

void results(Statistics &s) {
  // assert(c.check());
  s.note = "Lazy";
  s.n_leaves = c.num_of_buckets();
  s.n_capacity = c.capacity();
  s.n_internals = c.n_internals();
  s.max_depth = 2;
  s.slack = c.slack();
  s.in_size = INTERNAL_BSIZE;
  s.ln_size = c.size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
