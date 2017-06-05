#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>

#include <algorithm>
#include <vector>

#include "common.h"
#include "time_it.h"

using namespace std;

void initexp();
void destroyexp();

// op = 1: inserts the value.
void insert(long long value);

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value);

// op = 3: count values in range [a, b).
long long count(long long a, long long b);

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b);

static long long read_ll() {
  static char buf[1 << 20];
  static int i = 1 << 20;
  static int len = 1 << 20;
  if (i > (1 << 20) - 100) {
    for (int j = i; j < len; j++) {
      buf[j - i] = buf[j];
    }
    len -= i;
    i = 0;
    size_t read = fread(buf + len, sizeof(char), sizeof(buf) - len, stdin);
    len += read;
  }

  while (i < len && !isdigit(buf[i])) i++;
  // the_end = i >= len;
  long long a = 0;
  while (i < len && isdigit(buf[i])) a = a * 10 + (buf[i++] - '0');
  return a;
}

static vector<pair<long long, vector<long long>>> read_batch_by_op(
    long long &prev_op, unsigned &prev_size, const int n) {
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
  unsigned nth = 0, batch_size = 1;
  long long op = -1, chk = 0;
  double tot_runtime = 0;
  initexp();
  while (true) {
    auto batches = read_batch_by_op(op, batch_size, 100000);
    if (batches.size() == 0) break;
    for (auto batch : batches) {
      double runtime = run(batch.first, batch.second, chk);
      nth += batch.second.size() / (batch.first <= 2 ? 1 : 2);
      tot_runtime += runtime;
      printf("%d %.6lf %lld\n", nth, tot_runtime, chk);
    }
  }
  destroyexp();
}
