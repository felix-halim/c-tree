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
  bool ok = b.erase(value);
  assert(ok);
}

int query(int value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

void results(double insert_time, double query_time, int checksum) {
  printf("btree_stx_insert_time: %9.6lf, btree_stx_query_time: %9.6lf, btree_stx_csum: %d, ",
    insert_time, query_time, checksum);
}
