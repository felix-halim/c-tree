#include <cstdio>
#include <cassert>
#include "stx/btree_multiset"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

stx::btree_multiset<int> b;

static void init(int *arr, int N) {
  b.insert(arr, arr + N);
}

static void insert(int value) {
  b.insert(value);
}

static void erase(int value) {
  bool ok = b.erase(value);
  assert(ok);
}

static int query(int value) {
  auto it = b.lower_bound(value);
  return (it == b.end()) ? 0 : *it;
}

static void results(double insert_time, double query_time) {
  printf("btree_stx_insert_time: %9.6lf, btree_stx_query_time: %9.6lf, ",
    insert_time, query_time);
}
