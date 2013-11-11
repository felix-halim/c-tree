#include <cstdio>
#include <cassert>

#define LARGE_SIZE 4096
#define SMALL_SIZE 64

#define LARGE_TOUCH 0
#define SMALL_TOUCH 2

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
  s.note = "Lazy";
  s.n_leaves = n_large;//c.num_of_buckets();
  s.n_capacity = n_small;//c.capacity();
  s.n_internals = n_index;//c.n_internals();
  s.max_depth = c.n_chains();
  s.slack = 0;//c.slack();
  s.in_size = 0;//INTERNAL_BSIZE;
  s.ln_size = 0;//c.size();
  // c.alloc_sizes(s.ia_free, s.ia_size, s.la_free, s.la_size);
}
