#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <climits>

#include <algorithm>
#include <set>

#include "../common.h"
#include "../random.h"

// #define CHECK_UNIQUE

void random_input(int n) {
  Random random(140384);

  #ifdef CHECK_UNIQUE
    std::set<long long> s;
    fprintf(stderr, "checking for input uniqueness\n");
  #endif

  for (int i = 0; i < n; i++) {
    long long num = llabs(random.nextLong());
    printf("%d %lld\n", INSERT, num);

    #ifdef CHECK_UNIQUE
      assert(!s.count(num));
      s.insert(num);
    #endif
  }
}

void noup_range_query(int q, double selectivity) {
  long long max_value = (1ULL << 63) - 1;
  assert(max_value == LLONG_MAX);

  long long gap = std::max(1LL, (long long) (max_value * selectivity));

  Random random(140384);
  for (int i = 0; i < q; i++) {
    long long lo = llabs(random.nextLong());
    long long hi = lo + gap;
    if (hi < lo) hi = LLONG_MAX;
    printf("%d %lld %lld\n", SUM, lo, hi);
  }
}
