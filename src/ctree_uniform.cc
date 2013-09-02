#include <cstdio>
#include <cassert>
#include <algorithm>
#include "ctree_uniform.h"

#ifdef NOUP
  #include "test_noup.h"
#else
  #include "test_lfhv.h"
#endif

using namespace std;
using namespace ctree;

CTree c;

void init(int *arr, int N) {
  // c.load("ctree"); return;

  c.batch_insert(arr, N);
  // fprintf(stderr, "doi \n" );
  // for (int i = 0; i < N; i++) {
  //   c.insert(arr[i]);
  // }
  // assert(c.check());
  // c.optimize();
  // c.save("ctree");
  // locked = 1;
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
  assert(c.check());
  printf("%.6lf,%.6lf,%d,", insert_time, query_time, checksum);
  printf("\"%s\",%d,%d,%d,", version, nLeaves, nCap, nInternals);
  printf("%d,%d,%d,%d,", c.max_depth(), c.slack(), INTERNAL_BSIZE, LEAF_BSIZE);
  c.print_allocators();
  puts("");
  // c.print_stats();
}
