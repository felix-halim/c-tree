#include <cstdio>
#include <set>
#include <chrono>
#include <vector>
#include <functional>
#include <algorithm>
#include "comb2.h"
#include "random.h"

using namespace std;
using namespace comb;
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

  { "individual bit sortedness", [] {
    REP(i, 64) {
      Bucket<int> b;
      b.init();
      REP(j, 63) b.add_cracker_index(j, j);
      b.piece_set_sorted(i, true);
      REP(j, 64) {
        if (i == j) {
          ASSERT_TRUE(b.piece_is_sorted(j));
        } else {
          ASSERT_TRUE(!b.piece_is_sorted(j));
        }
      }
    }
  }},

  { "random bits sortedness", [] {
    Random r(140384);
    REP(t, 1000) {
      set<int> on;
      int nbits = r.nextInt(64);
      REP(i, nbits) on.insert(r.nextInt(64));

      Bucket<int> b;
      b.init();
      REP(j, 63) b.add_cracker_index(j, j);
      for (int i : on)
        b.piece_set_sorted(i, true);
      REP(i, 64) {
        if (on.count(i)) {
          ASSERT_TRUE(b.piece_is_sorted(i));
        } else {
          ASSERT_TRUE(!b.piece_is_sorted(i));
        }
      }
    }
  }},

  { "invalidate bits sortedness", [] {
    Random r(140384);
    REP(t, 1000) {
      set<int> on;
      int nbits = r.nextInt(64);
      REP(i, nbits) on.insert(r.nextInt(64));

      Bucket<int> b;
      b.init();
      REP(j, 63) b.add_cracker_index(j, j);
      for (int i : on)
        b.piece_set_sorted(i, true);

      int j = r.nextInt(64);
      b.piece_set_unsorted_onwards(j); // Invalidate bits j onwards.

      REP(i, 64) {
        if (on.count(i) && i < j) {
          ASSERT_TRUE(b.piece_is_sorted(i));
        } else {
          ASSERT_TRUE(!b.piece_is_sorted(i));
        }
      }
    }
  }},


  { "add bits and sortedness", [] {
    set<int> on { 2, 10, 17, 18, 24, 48, 50, 62 };

    Bucket<int> b;
    b.init();
    REP(j, 62) b.add_cracker_index(j, j * 5);
    for (int i : on)
      b.piece_set_sorted(i, true);

    REP(i, 63) {
      if (on.count(i)) {
        ASSERT_TRUE(b.piece_is_sorted(i));
      } else {
        ASSERT_TRUE(!b.piece_is_sorted(i));
      }
    }

    b.add_cracker_index(18, 18 * 5 - 1);

    set<int> new_on;
    for (int i : on) new_on.insert(i < 18 ? i : (i + 1));

    REP(i, 64) {
      if (new_on.count(i)) {
        ASSERT_TRUE(b.piece_is_sorted(i));
      } else {
        ASSERT_TRUE(!b.piece_is_sorted(i));
      }
    }
  }},

  { "remove bits and sortedness", [] {
    set<int> on { 2, 10, 17, 18, 19, 48, 50, 62 };

    Bucket<int> b;
    b.init();
    REP(j, 62) b.add_cracker_index(j, j * 5);
    for (int i : on)
      b.piece_set_sorted(i, true);

    REP(i, 63) {
      if (on.count(i)) {
        ASSERT_TRUE(b.piece_is_sorted(i));
      } else {
        ASSERT_TRUE(!b.piece_is_sorted(i));
      }
    }

    b.remove_cracker_index(18);

    on.erase(18);
    set<int> new_on;
    for (int i : on) new_on.insert(i < 18 ? i : (i - 1));

    REP(i, 62) {
      if (new_on.count(i)) {
        ASSERT_TRUE(b.piece_is_sorted(i));
      } else {
        ASSERT_TRUE(!b.piece_is_sorted(i));
      }
    }
  }},

  { "flush pending insert", [] {

  }},

  { "rough middle partition", [] {
    // Random rng(140384);
    // std::less<int> cmp;
    // int arr[1024];
    // REP(t, 2000) {
    //   int N = rng.nextInt(1024);
    //   int D = rng.nextInt(100000);
    //   REP(i, N) arr[i] = rng.nextInt(D);
    //   rough_middle_partition(arr, arr + N, cmp, rng);
    //   int val = arr[N / 2];
    //   REP(i, N) {
    //     if (i < N / 2) {
    //       ASSERT_TRUE(arr[i] <= val);
    //     } else {
    //       ASSERT_TRUE(arr[i] >= val);
    //     }
    //   }
    // }
  }},

  { "fusion", [] {
  }},

  { "crack", [] {
  }},

  { "erase", [] {
  }},

  { "comb small", [] {
    Random rng;
    constexpr int N = 1000;
    int arr[N];
    REP(i, N) arr[i] = rng.nextInt(100);

    Comb<int, std::less<int>, false, 256, 32, 8> c;
    REP(i, N) c.insert(arr[i]);

    random_shuffle(arr, arr + N);
    int val, csum = 0;
    REP(i, N) {
      auto it = c.lower_bound(arr[i]);
      ASSERT_TRUE(it.next(val));
      csum = csum * 13 + val;
    }
    fprintf(stderr, "csum = %d\n", csum);
  }},

  { "comb correctness", [] {
    int val;
    Random rng;
    Comb<int> c;
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
          ASSERT_TRUE(it.next(val));
          ASSERT_TRUE(val == num);
        } else {
          auto it = c.lower_bound(num);
          ASSERT_TRUE(!it.next(val) || val != num);
        }
      }
    }
    printf(" ");
  }},

  { "comb", [] {
    Random rng;
    Comb<int> c;
    double insert_time = time_it([&] {
      REP(i, 1000000)
        c.insert(rng.nextInt());
    });

    double query_time = time_it([&] {
      int csum = 0, val;
      REP(i, 10000) {
        auto it = c.lower_bound(10);
        if (it.next(val)) csum = csum * 13 + val;
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
