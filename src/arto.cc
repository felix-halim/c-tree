#include <cstdio>
#include <cassert>
#include <algorithm>
#include "arto.h"
#include "test.h"

inline uintptr_t getLeafValue(Node* node) {
   // The the value stored in the pseudo-leaf
   return reinterpret_cast<uintptr_t>(node)>>1;
}

using namespace std;

Node* tree = NULL;

void init(long long *arr, int N) {
  for (int i = 0; i < N; i++) {
    insert(arr[i]);
  }
  always_flush = 1;
  // art_debug = 1;
}

void insert(long long value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  insert(tree,&tree,key,0,value64,8);
}

void erase(long long value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  // assert(lookup(tree,key,8,0,8));
  erase(tree,&tree,key,8,0,8);
  // assert(!lookup(tree,key,8,0,8));
}

long long query(long long value) {
  // art_debug = value == 754275843;
  // art_debug = 1;
  ART_DEBUG("\nquery %lld\n", value);
  // fprintf(stderr, "query %lld\n", value);
  uint64_t value64 = value;
  uint8_t key[8];
  loadKey(value64, key);
   
  // if (value == 61554031375105761) art_debug = 1;
  Node* leaf=lower_bound(tree,&tree,key,8,0,8);
  long long ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stdout, "%lld (%lld)\n", ret, value);
  }
    // art_debug = 0;
  return ret;
}

void results(Statistics &s) {
}
