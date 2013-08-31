#include <cstdio>
#include <cassert>
#include "test_noup.h"

int *arr, N;
double t;

void init(int *iarr, int iN) {
  t = time_it([&] {
    N = iN;
    arr = new int[N];
    for (int i = 0; i < N; i++)
      arr[i] = iarr[i];
  });
  sort(arr, arr + N);
}

int query(int value) {
  auto it = lower_bound(arr, arr + N, value);
  return (it == arr + N) ? 0 : *it;
}

void results(double insert_time, double query_time, int checksum) {
  printf("%.6lf,%.6lf,%d,%.6lf\n", insert_time, query_time, checksum, t);
}
