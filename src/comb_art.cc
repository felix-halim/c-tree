#include <cstdio>
#include <cassert>

#define LARGE_SIZE 4096
#define SMALL_SIZE 64

#ifndef LARGE_TOUCH
  #define LARGE_TOUCH 1000000000
#endif

#ifndef SMALL_TOUCH
  #define SMALL_TOUCH 1000000000
#endif

#include "comb_art.h"
#include "tester.h"

inline uintptr_t getLeafValue(Node* node) {
  // The the value stored in the pseudo-leaf
  uintptr_t v = getData((uintptr_t) node);
  return isPointer(v) ? data((Bucket<unsigned>*) v) : getData(v);
}

Comb<long long> c;

void initexp() {}
void destroyexp() {}

// op = 1: inserts the value.
void insert(long long value) {
  c.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  #ifdef NDEBUG
    c.erase(value);
  #else
    bool ok = c.erase(value);
    if (!ok) fprintf(stderr, "ERASE error %d\n", value);
    assert(ok);
  #endif
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  auto it1 = c.lower_bound(a);
  auto it2 = c.lower_bound(b);
  return 0;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  return c.lower_bound(a);
  // auto it1 = c.lower_bound(a);
  // Crack the lower bound so that it is in the correct position.
  // c.lower_bound(b);
  // long long sum = 0;
  // while (it1.has_next()) {
  //   long long v = *it1.next();
  //   if (v >= b) break;
  //   // fprintf(stderr, "got %lld\n", v); 
  //   sum += v;
  // }
  // return sum;
}

  // assert(c.check());
    // s.N = N;
    // s.n_index = n_index;
    // s.n_bytes = n_bytes;
    // s.n_slack_int = n_slack_art;
    // s.n_slack_leaf = n_slack_leaves;
    // s.n_internal = n_internal;
    // s.n_leaf = n_leaf;
    // s.n_small = n_small;
    // s.n_large = n_large;
    // s.n_chained = n_chain;
    // s.art_n4 = art_n4;
    // s.art_n16 = art_n16;
    // s.art_n48 = art_n48;
    // s.art_n256 = art_n256;
  // s.large_touch = LARGE_TOUCH;
  // s.small_touch = SMALL_TOUCH;
