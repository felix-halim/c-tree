#include <cstdio>
#include <cassert>
#include "comb.h"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

// Comb<int, std::less<int>, true, 3200, 125, 50> c;
Comb<int> c;

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

static void results(double insert_time, double query_time) {
  printf("comb_insert_time: %9.6lf, comb_query_time: %9.6lf, ",
    insert_time, query_time);
}
