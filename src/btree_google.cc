#include <cstdio>
#include <cassert>
#include "google/btree_set.h"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

btree::btree_multiset<int> b;

static void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    b.insert(arr[i]);
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

static void results(double insert_time, double query_time, int checksum) {
  printf("btree_google_insert_time: %9.6lf, btree_google_query_time: %9.6lf, btree_google_csum: %d, ",
    insert_time, query_time, checksum);
}
