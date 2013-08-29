#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree_exp_leafsize.h"

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
  assert(c.check());
  // while (c.optimize());
  // c.compact();
  // while (c.optimize());
  // // c.debug();
  // fprintf(stderr, "depth = %d, slack = %d\n", c.max_depth(), c.slack());
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
  printf("%.6lf,%.6lf,%d,", insert_time, query_time, checksum);
  printf("\"%s\",%d,%d,%d,%.6lf,%.6lf,%.6lf\n", c.version, nLeaves, nCap, nInternals, c.t1, c.t2, c.t3);
  // fprintf(stderr, "depth = %d, slack = %d\n", c.max_depth(), c.slack());
  assert(c.check());
}
