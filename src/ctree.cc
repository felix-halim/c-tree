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

static void init(int *arr, int N) {
  for (int i = 0; i < N; i++) {
    c.insert(arr[i]);
    if (!(i & (i + 1))) fprintf(stderr, ".");
  }
}

static void insert(int value) {
  c.insert(value);
}

static void erase(int value) {
  bool ok = c.erase(value);
  assert(ok);
}

static int query(int value) {
  auto it = c.lower_bound(value);
  return it.first ? it.second : 0;
}

static void results(double insert_time, double query_time, int checksum) {
  printf("ctree_insert_time: %9.6lf, ctree_query_time: %9.6lf, ctree_csum: %d, ",
    insert_time, query_time, checksum);
  fprintf(stderr, "made = %d\n", c.depth);
}
