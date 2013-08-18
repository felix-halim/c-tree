#include <cstdio>
#include <cassert>
#include "comb2.h"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

using namespace comb;

// Comb<int, std::less<int>, true, 3200, 125, 50> c;
Comb<int, std::less<int>, false, 4098, 256, 128> c;

static void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    c.insert(arr[i]);
}

static void insert(int value) {
  c.insert(value);
}

static void erase(int value) {
  bool ok = c.erase(value);
  assert(ok);
}

static int query(int value) {
  int val = 0;
  return c.lower_bound(value).next(val) ? val : 0;
}

static void results(double insert_time, double query_time, int checksum) {
  printf("comb2_insert_time: %9.6lf, comb2_query_time: %9.6lf, comb2_csum: %d, ",
    insert_time, query_time, checksum);
  // fprintf(stderr, "crack_time = %8.4lf, ", crack_time);
  // fprintf(stderr, "sort_time = %8.4lf, ", sort_time);
  // fprintf(stderr, "lower_time = %8.4lf\n", lower_time);
}
