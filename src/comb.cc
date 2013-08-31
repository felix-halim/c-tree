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

void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    c.insert(arr[i]);
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  // fprintf(stderr, "erase %d\n", value);
  bool ok = c.erase(value);
  assert(ok);
}

int query(int value) {
  int val = 0;
  int ret = c.lower_bound(value).next(val) ? val : 0;
  fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

void results(double insert_time, double query_time, int checksum) {
  printf("%.6lf,%.6lf,%d\n", insert_time, query_time, checksum);
}
