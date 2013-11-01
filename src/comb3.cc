#include <cstdio>
#include <cassert>
#include "comb3.h"
#include "test.h"

Comb<int> c;

void init(int *arr, int N) {
  c.load(arr, N);
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
  int val = 0;
  int ret = c.lower_bound(value).next(val) ? val : 0;
  // fprintf(stdout, "%d (%d)\n", ret, value);
  return ret;
}

int select(int a, int b) {
  int aa, bb;
  int ret1 = c.lower_bound(a).next(aa) ? aa : 0;
  int ret2 = c.lower_bound(b).next(bb) ? bb : 0;
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
