#include <cstdio>
#include <cassert>
#include "comb.h"
#include "test.h"

#if defined(COMB800)
  Comb<unsigned, std::less<unsigned>, false, 800, 30, 12> c(100000000);
#elif defined(COMB1600)
  Comb<unsigned, std::less<unsigned>, false, 1600, 62, 25> c(100000000);
#elif defined(COMB3200)
  Comb<unsigned, std::less<unsigned>, false, 3200, 125, 50> c(100000000);
#elif defined(COMB6400)
  Comb<unsigned> c(100000000);
#endif

void init(unsigned *arr, unsigned N) {
  c.load(arr, N);
}

void insert(unsigned value) {
  c.insert(value);
}

void erase(unsigned value) {
  // fprintf(stderr, "erase %d\n", value);
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

unsigned lower_bound(unsigned value) {
  unsigned val = 0;
  unsigned ret = c.lower_bound(value).next(val) ? val : 0;
  // printf("%d (%d)\n", ret, value);
  return ret;
}

unsigned select(unsigned a, unsigned b) {
  unsigned aa, bb;
  unsigned ret1 = c.lower_bound(a).next(aa) ? aa : 0;
  unsigned ret2 = c.lower_bound(b).next(bb) ? bb : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

unsigned count(unsigned a, unsigned b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  return c.count(it1, it2);
}

void results(Statistics &s) {
  // assert(c.check());
  s.n_index = c.root_size();
  s.n_bytes = c.num_of_buckets() * c.bucket_size() * sizeof(unsigned);
  s.n_leaf = c.num_of_buckets();
  s.n_slack_int = c.slack();
  s.n_large = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
