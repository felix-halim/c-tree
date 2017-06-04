#include <cstdio>
#include <cassert>
#include "comb.h"
#include "tester.h"

#if defined(COMB800)
  Comb<unsigned, std::less<unsigned>, false, 800, 30, 12> c(100000000);
#elif defined(COMB1600)
  Comb<unsigned, std::less<unsigned>, false, 1600, 62, 25> c(100000000);
#elif defined(COMB3200)
  Comb<unsigned, std::less<unsigned>, false, 3200, 125, 50> c(100000000);
#else
  Comb<unsigned> c(100000000);
#endif

void initexp() {}
void destroyexp() {}

// op = 1: inserts the value.
void insert(long long value) {
  c.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  return c.erase(value);
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  return c.count(it1, it2);
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  long long sum = 0;
  while (it1 != it2) {
    sum += *it1.next();
  }
  return sum;
}


  // c.load(arr, N);

  // unsigned val = 0;
  // unsigned ret = c.lower_bound(value).next(val) ? val : 0;

  // assert(c.check());
  // s.n_index = c.root_size();
  // s.n_bytes = c.num_of_buckets() * c.bucket_size() * sizeof(unsigned);
  // s.n_leaf = c.num_of_buckets();
  // s.n_slack_int = c.slack();
  // s.n_large = c.bucket_size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
