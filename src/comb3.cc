#include <cstdio>
#include <cassert>
#include "comb3.h"
#include "test.h"

Comb<int, less<int>, 2048> c;

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
  s.n_leaves = 0;//c.num_of_buckets();
  s.n_capacity = 0;//c.num_of_buckets() * c.bucket_size();
  s.n_internals = 1;
  s.max_depth = 2;
  s.in_size = 0;//c.root_size();
  s.ln_size = 0;//c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
