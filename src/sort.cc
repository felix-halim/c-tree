#include <cstdio>
#include <cassert>
#include "test_noup.h"

static int *arr, N;

static void init(int *iarr, int iN) {
  N = iN;
  arr = new int[N];
  for (int i = 0; i < N; i++)
    arr[i] = iarr[i];
  sort(arr, arr + N);
}

static int query(int value) {
  auto it = lower_bound(arr, arr + N, value);
  return (it == arr + N) ? 0 : *it;
}

static void results(double insert_time, double query_time) {
  printf("sort_insert_time: %9.6lf, sort_query_time: %9.6lf, ",
    insert_time, query_time);
}
