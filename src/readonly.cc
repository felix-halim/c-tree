#include <cstdio>
#include <cassert>

#include <memory>
#include <vector>

#include "tester.h"

using namespace std;

class Bucket {
 public:
  long long arr[1 << 13];
  int idx;
  int full;

  bool add(long long val) {
    if (full) return false;
    arr[idx++] = val;
    full = idx == (1 << 13);
    return true;
  }
};

vector<std::unique_ptr<Bucket>> buckets;
std::unique_ptr<Bucket> last;

vector<long long> arr(100000000);
unsigned idx = 0;
long long total = 0;

void initexp() {
  arr.clear();
  last = make_unique<Bucket>();
}

void destroyexp() {}

// op = 1: inserts the value.
void insert(long long value) {
  if (!last->add(value)) {
    buckets.push_back(std::move(last));
    last = make_unique<Bucket>();
    last->add(value);
  }
  // arr.push_back(value);
  // arr[idx++] = value;
  // assert(idx <= arr.size());
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  return true;
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  return -1;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  return arr.size() + buckets.size();
}
