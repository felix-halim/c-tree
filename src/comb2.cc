#include <cstdio>
#include <cassert>
#include "comb2.h"

#ifdef NOUP
  // #include "test_noup.h"
#else
  // #include "test_lfhv.h"
#endif

using namespace comb;

// Comb<int, std::less<int>, true, 3200, 125, 50> c;
Comb<int, std::less<int>, false, 4098, 256, 128> c;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++)
    c.insert(arr[i]);
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  bool ok = c.erase(value);
  assert(ok);
}

int query(int value) {
  int val = 0;
  return c.lower_bound(value).next(val) ? val : 0;
}

int select(int a, int b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  int ret1 = it1.first ? it1.second : 0;
  int ret2 = it2.first ? it2.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret1 + ret2;
}

void results(double insert_time, double query_time, int checksum) {
  printf("%.6lf,%.6lf,%d\n", insert_time, query_time, checksum);
}
