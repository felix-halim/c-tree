#include <cstdio>
#include "multiset.h"
#include "tester.h"

trimmer::multiset<unsigned, std::less<unsigned>, Mallocator> s;

// op = 1: inserts the value.
void insert(long long value) {
  s.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  return s.erase(value);
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  // auto it1 = s.lower_bound(a);
  // auto it2 = s.lower_bound(b);
  return 0; // it2 - it1;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  auto it1 = s.lower_bound(a);
  return *it1;
  // auto it2 = s.lower_bound(b);
  // long long sum = 0;
  // while (!(it1 == it2)) {
  //   sum += *(it1++);
  // }
  // return sum;
}

// s.insert(arr, arr + N);

// assert(s.check());
// s.n_index = s.root_size();
// s.n_bytes = s.num_of_buckets() * s.bucket_size() * sizeof(unsigned);
// s.n_leaf = s.num_of_buckets();
// s.n_slack_int = s.slack();
// s.n_large = s.bucket_size();
// s.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
