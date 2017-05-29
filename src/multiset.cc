#include <cstdio>
#include <cassert>
#include "multiset.h"
#include "test.h"

using namespace trimmer;

trimmer::multiset<unsigned, std::less<unsigned>, Mallocator> c;

void init(unsigned *arr, unsigned N) {
  c.insert(arr, arr + N);
}

void insert(unsigned value) {
  c.insert(value);
}

void erase(unsigned value) {
  fprintf(stderr, "erase %d\n", value);
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    assert(ok);
  #endif
}

unsigned lower_bound(unsigned value) {
  auto it = c.lower_bound(value);
  // printf("%d (%d)\n", value, value);
  // return 0;
  return it == c.end() ? 0 : *it;
}

unsigned select(unsigned a, unsigned b) {
  return lower_bound(a) + lower_bound(b);
}

unsigned count(unsigned a, unsigned b) {
  return 0;
}

void results(Statistics &s) {
  // assert(c.check());
  // s.n_index = c.root_size();
  // s.n_bytes = c.num_of_buckets() * c.bucket_size() * sizeof(unsigned);
  // s.n_leaf = c.num_of_buckets();
  // s.n_slack_int = c.slack();
  // s.n_large = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
