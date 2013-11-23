#include <cstdio>
#include <cassert>

#define LARGE_SIZE 4096
#define SMALL_SIZE 64

#define LARGE_TOUCH 10
#define SMALL_TOUCH 10

#include "comb_art.h"
#include "test.h"

Comb c;

void init(int *arr, int N) {
  // art_debug = 1;
  c.load(arr, N);
  // for (int i = 0; i < N; i++) {
  //   insert(arr[i]);
  //   if (i % 1000) lower_bound(arr[i]);
  // }
}

void insert(int value) {
  c.insert(value);
}

void erase(int value) {
  // fprintf(stderr, "erase %d\n", value);
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    if (!ok) fprintf(stderr, "ERASE error %d\n", value);
    assert(ok);
  #endif
}

int lower_bound(int value) {
  int ret = c.lower_bound(value);
  // fprintf(stdout, "%d (%d)\n", ret, value); fflush(stdout);
  return ret;
}

int select(int a, int b) {
  return 0;
}

int count(int a, int b) {
  return 0;
}

void results(Statistics &s) {
  // assert(c.check());
  s.n_index = n_index;
  s.n_bytes = 0;
  s.n_slack = 0;
  s.n_internal = n_small;
  s.n_leaf = n_large;
  s.bt_int_sz = 0;
  s.bt_leaf_sz = 0;
  s.large_touch = LARGE_TOUCH;
  s.small_touch = SMALL_TOUCH;
}
