#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

#define LEAF_BSIZE 12

#include "artd_crack.h"

using namespace art_crack;

inline uintptr_t getLeafValue(Node* node) {
  // The the value stored in the pseudo-leaf
  return reinterpret_cast<uintptr_t>(node) >> 1;
}

ArtCrack<long long> c;

int main() {
  Random r(140384);
  art_debug = 1;
  vector<long long> arr;
  for (int i = 1; i < 10000; i++) {
    long long v = r.nextInt();
    if (find(arr.begin(), arr.end(), v) != arr.end()) { i--; continue; }
    arr.push_back(v);
    c.insert(v);
  }
  for (int i = 1; i < 100000; i++) {
    long long v = r.nextInt();
    auto it = c.lower_bound(v);
    sort(arr.begin(), arr.end());
    int N = lower_bound(arr.begin(), arr.end(), v) - arr.begin();
    if (N == arr.size()) {
      assert(!it.first);
    } else {
      fprintf(stderr, "got %lld, ex %lld\n", it.second, arr[N]);
      assert(it.first);
      assert(it.second == arr[N]);
    }

    int j = r.nextInt(arr.size());
    bool ok = c.erase(arr[j]);
    assert(ok);
    arr[j] = r.nextInt();
    c.insert(arr[j]);
  }
}
