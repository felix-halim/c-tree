#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <algorithm>
#include <chrono>
#include <vector>

#include "common.h"

using namespace std;
using namespace std::chrono;

// op = 1: inserts the value.
void insert(long long value);

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value);

// op = 3: count values in range [a, b).
long long count(long long a, long long b);

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b);

static long long read_ll() {
  long long a;
  return scanf("%lld", &a) != EOF ? a : 0;
}

static vector<pair<long long, vector<long long>>> read_batch_by_op(
    long long &prev_op, int &prev_size, const int n) {
  vector<pair<long long, vector<long long>>> ret;
  vector<long long> arr;
  for (long long i = 0; i < n; i++) {
    long long op = read_ll();
    if (prev_op != op) {
      if (prev_op > 0) {
        if (arr.size() > 0) {
          ret.push_back(make_pair(prev_op, arr));
          arr.clear();
        }
        prev_size = 1;
      }
      prev_op = op;
    }
    if (op == 0) break;
    if (arr.size() * 2 >= prev_size) {
      assert(prev_op == op);
      ret.push_back(make_pair(prev_op, arr));
      prev_size *= 2;
      arr.clear();
    }
    switch (op) {
      case INSERT: arr.push_back(read_ll()); break;
      case ERASE: arr.push_back(read_ll()); break;
      case COUNT: arr.push_back(read_ll()); arr.push_back(read_ll()); break;
      case SUM: arr.push_back(read_ll()); arr.push_back(read_ll()); break;
      default: assert(false);
    }
  }
  if (prev_op > 0 && arr.size() > 0) {
    ret.push_back(make_pair(prev_op, arr));
  }
  return ret;
}

template<typename F>
static double time_it(F f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

static double run(long long op, const vector<long long> &arr, long long &chk) {
  switch (op) {
    case INSERT: return time_it([&] { for (long long a : arr) insert(a); });
    case ERASE: return time_it([&] {
        for (long long a : arr) {
          chk = chk * 3 + (erase(a) ? 1 : 0);
        }
      });
    case COUNT: return time_it([&] {
        auto i = arr.begin();
        while (i != arr.end()) {
          chk = chk * 13 + count(*i, *(i + 1));
          i += 2;
        }
      });
    case SUM: return time_it([&] {
        auto i = arr.begin();
        while (i != arr.end()) {
          chk = chk * 13 + sum(*i, *(i + 1));
          i += 2;
        }
      });
  }
  assert(false);
  return 0;
}

int main(int argc, char *argv[]) {
  int nth = 0, batch_size = 1;
  long long op = -1, chk = 0;
  double tot_runtime = 0;
  while (true) {
    auto batches = read_batch_by_op(op, batch_size, 10000);
    if (batches.size() == 0) break;
    for (auto batch : batches) {
      double runtime = run(batch.first, batch.second, chk);
      nth += batch.second.size() / (batch.first <= 2 ? 1 : 2);
      tot_runtime += runtime;
      printf("%d %.6lf %lld\n", nth, tot_runtime, chk);
    }
  }
}
