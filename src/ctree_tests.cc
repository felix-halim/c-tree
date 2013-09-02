#include <cstdio>
#include <set>
#include <chrono>
#include <vector>
#include <functional>
#include <algorithm>
#include "ctree_uniform.h"
#include "random.h"

using namespace std;
using namespace ctree;
using namespace chrono;

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)
#define PASSED "\033[92mPASSED\033[0m"
#define FAILED "\033[91mFAILED\033[0m"
#define ASSERT_TRUE(v)  \
  if (!(v)) {           \
    printf(FAILED "\nLine %d: assert failed\n", __LINE__); \
    fflush(stdout);     \
    abort();            \
  }

vector<pair<string,function<void()>>> tests {

  { "ctree erase small", [] {
    Random rng;
    int N = 70, E = 30;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt(100);

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.debug();

    random_shuffle(arr, arr + N);

    REP(i, E) {
      bool ok = c.erase(arr[--N]);
      assert(ok);
    }
    // c.debug();

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

  { "ctree erase medium", [] {
    Random rng;
    int N = 1000, E = 800;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt(1000);

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.debug();

    random_shuffle(arr, arr + N);

    REP(i, E) {
      bool ok = c.erase(arr[--N]);
      assert(ok);
    }

    // c.debug();

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

  { "ctree small", [] {
    Random rng;
    constexpr int N = 70;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt(100);

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.optimize();
    c.debug();

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

  { "ctree medium", [] {
    Random rng;
    constexpr int N = 100;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt(300);

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.optimize();
    c.debug();

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

  { "ctree large", [] {
    Random rng;
    constexpr int N = 200;
    int *arr = new int[N];
    REP(i, N) arr[i] = rng.nextInt(300);

    CTree c;
    REP(i, N) c.insert(arr[i]);
    // c.optimize();
    c.debug();

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

  { "ctree erase correctness", [] {
    Random rng;
    CTree c;
    multiset<int> mset;
    REP(i, 10000) {
      if (!(i & (i + 1))) printf("."), fflush(stdout);
      int num = rng.nextInt(100);
      if (mset.size() > 100) {
        auto it = mset.find(num);
        bool ok = c.erase(num);
        // c.debug();
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
          // c.debug();
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

  { "ctree", [] {
    Random rng;
    CTree c;
    double insert_time = time_it([&] {
      REP(i, 1000000)
        c.insert(rng.nextInt());
    });
    c.optimize();
    // c.debug();

    double query_time = time_it([&] {
      int csum = 0;
      REP(i, 10000) {
        auto it = c.lower_bound(10);
        if (it.first) csum = csum * 13 + it.second;
      }
      fprintf(stderr, "csum = %d\n", csum);
    });
    fprintf(stderr, "insert = %6.3lf, query = %6.3lf\n", insert_time, query_time);
  }},

};


int main() {
  for (auto it : tests) {
    // if (strstr(it.first.c_str(), "erase")) continue;
    printf("%s ... ", it.first.c_str());
    it.second();
    puts(PASSED);
  }
}
