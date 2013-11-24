#include <cstdio>
#include <cassert>
#include <algorithm>
#include "art.h"
#include "test.h"

inline uintptr_t getLeafValue(Node* node) {
   // The the value stored in the pseudo-leaf
   return reinterpret_cast<uintptr_t>(node)>>1;
}

using namespace std;

Node* tree = NULL;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++) {
    insert(arr[i]);
  }
  // fprintf(stderr, "hash = %llu\n", (unsigned long long) hash_tree(tree));
  // art_debug = 1;
}

void insert(uintptr_t value64) {
  // static int nth = 0; nth++;
  uint8_t key[8];
  // uint64_t value64 = value;
  // value64 = (value64 << 30) | nth;
  loadKey(value64, key);
  insert(tree,&tree,key,0,value64,8);
}

void erase(int value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  // assert(lookup(tree,key,8,0,8));
  erase(tree,&tree,key,8,0,8);
  // assert(!lookup(tree,key,8,0,8));
}

int lower_bound(int value) {
  uint64_t value64 = value;
  value64 = (value64 << 30);
  uint8_t key[8];
  loadKey(value64, key);
   
  Node* leaf=lower_bound(tree,key,8,0,8);
  int ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stdout, "%lld (%lld)\n", ret, value);
  }
  return ret;
}

int select(int a, int b) {
  return lower_bound(a) + lower_bound(b);
}

void results(Statistics &s) {
}
