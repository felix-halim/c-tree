#include <cstdio>
#include <cassert>
#include "stx/btree_multiset"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

stx::btree_multiset<int> b;

void init(int *arr, int N) {
  b.insert(arr, arr + N);
}

void insert(int value) {
  b.insert(value);
}

void erase(int value) {
  auto it = b.lower_bound(value);
  assert(it != b.end());
  b.erase(it);
}

int query(int value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

void results(double insert_time, double query_time, int checksum) {
  printf("%.6lf,%.6lf,%d\n", insert_time, query_time, checksum);
}
