#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree.h"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

using namespace std;
using namespace ctree;

CTree c;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++) {
    c.insert(arr[i]);
  }
  // c.debug();
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  bool ok = c.erase(value);
  assert(ok);
}

int query(int value) {
  auto it = c.lower_bound(value);
  int ret = it.first ? it.second : 0;
  // fprintf(stderr, "%d (%d)\n", ret, value);
  return ret;
}

void results(double insert_time, double query_time, int checksum) {
  printf("ctree_insert_time: %9.6lf, ctree_query_time: %9.6lf, ctree_csum: %d, ",
    insert_time, query_time, checksum);
  fprintf(stderr, "t1 = %.6lf\n", c.t1);
  fprintf(stderr, "t2 = %.6lf\n", c.t2);
  fprintf(stderr, "t3 = %.6lf\n", c.t3);
  fprintf(stderr, "nLeaves = %d, nDes = %d\n", nLeaves, nDes);
  fprintf(stderr, "nCap = %d\n", nCap);
  fprintf(stderr, "nInternals = %d\n", nInternals);
}
