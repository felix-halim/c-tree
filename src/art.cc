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
  // n4 = 0;
  // n16 = 0;
  // n48 = 0;
  // n256 = 0;
  // // art_debug = 1;
  // nsplit = 0;
  // nadv = 0;
}

void insert(int value) {
  uint8_t key[8];
  uint64_t value64 = value;
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

int query(int value) {
  uint64_t value64 = value;
  uint8_t key[8];
  loadKey(value64, key);
   
  // if (value == 61554031375105761) art_debug = 1;
  Node* leaf=lower_bound(tree,key,8,0,8);
  int ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stdout, "%lld (%lld)\n", ret, value);
  }
    // art_debug = 0;
  return ret;
}

void results(Statistics &s) {
  s.in_size = nsplit;
  s.ln_size = nadv;
  s.ia_free = n4;
  s.ia_size = n16;
  s.la_free = n48;
  s.la_size = n256;
}
