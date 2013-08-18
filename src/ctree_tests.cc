#include <cstdio>
#include <set>
#include <chrono>
#include <vector>
#include <functional>
#include <algorithm>
#include "ctree.h"
#include "random.h"

using namespace std;
using namespace ctree;
using namespace chrono;

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)
#define PASSED "\033[92mPASSED\033[0m"
#define FAILED "\033[91mFAILED\033[0m"
#define ASSERT_TRUE(v)  \
  if (!(v)) {           \
    printf(FAILED "\nLine %d: assert failed\n", __LINE__); \
    fflush(stdout);     \
    abort();            \
  }

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

vector<pair<string,function<void()>>> tests {

  { "comb small", [] {
    Random rng;
    constexpr int N = 100000000;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt();

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.debug();

    random_shuffle(arr, arr + N);
    int csum = 0;
    REP(i, N) {
      auto it = c.lower_bound(arr[i]);
      // fprintf(stderr, "find %d, res = %d, %d\n", arr[i], it.first, it.second);
      ASSERT_TRUE(it.first);
      ASSERT_TRUE(it.second == arr[i]);
      csum = csum * 13 + it.second;
    }
    fprintf(stderr, "csum = %d\n", csum);
  }},

  { "comb correctness", [] {
    return;
    int val;
    Random rng;
    CTree c;
    multiset<int> mset;
    REP(i, 1000000) {
      if (!(i & (i + 1))) printf("."), fflush(stdout);
      int num = rng.nextInt(100000);
      if (mset.size() > 100000) {
        auto it = mset.find(num);
        bool ok = c.erase(num);
        if (it != mset.end()) {
          mset.erase(it);
          ASSERT_TRUE(ok);
        } else {
          ASSERT_TRUE(!ok);
        }
      } else {
        if (rng.nextBoolean()) {
          mset.insert(num);
          c.insert(num);
        } else if (mset.count(num)) {
          auto it = c.lower_bound(num);
          ASSERT_TRUE(it.first);
          ASSERT_TRUE(it.second == num);
        } else {
          auto it = c.lower_bound(num);
          ASSERT_TRUE(!it.first || it.second != num);
        }
      }
    }
    printf(" ");
  }},

  { "comb", [] {
    Random rng;
    CTree c;
    double insert_time = time_it([&] {
      REP(i, 1000000)
        c.insert(rng.nextInt());
    });

    double query_time = time_it([&] {
      int csum = 0, val;
      REP(i, 10000) {
        auto it = c.lower_bound(10);
        csum = csum * 13 + it.second;
      }
      fprintf(stderr, "csum = %d\n", csum);
    });
    fprintf(stderr, "insert = %6.3lf, query = %6.3lf\n", insert_time, query_time);
  }},

};


int main() {
  for (auto it : tests) {
    printf("%s ... ", it.first.c_str());
    it.second();
    puts(PASSED);
  }
}
