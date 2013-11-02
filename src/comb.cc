#include <cstdio>
#include <cassert>
#include "comb.h"
#include "test.h"

#if defined(COMB800)
  Comb<int, std::less<int>, false, 800, 30, 12> c(100000000);
#elif defined(COMB1600)
  Comb<int, std::less<int>, false, 1600, 62, 25> c(100000000);
#elif defined(COMB3200)
  Comb<int, std::less<int>, false, 3200, 125, 50> c(100000000);
#elif defined(COMB6400)
  Comb<int> c(100000000);
#endif

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
  fprintf(stdout, "%d (%d)\n", ret, value);
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
  s.n_leaves = c.num_of_buckets();
  s.n_capacity = c.num_of_buckets() * c.bucket_size();
  s.n_internals = 1;
  s.max_depth = 2;
  s.slack = c.slack();
  s.in_size = c.root_size();
  s.ln_size = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
