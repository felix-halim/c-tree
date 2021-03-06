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
int N;

void init(unsigned *arr, unsigned n) {
  N = n;
  #ifdef EAGER
    bulk_insert(tree, arr, N); // Lazy insert, chain buckets.
  #else
    pending_bulk_insert(tree, arr, N); // Lazy insert, chain buckets.
  #endif
  // art_debug = 1;
}

void insert(unsigned value64) {
  uint8_t key[8];
  loadKey(value64, key);
  // insert(&tree,key,0,value64,8);
  insert(tree,&tree,key,0,value64,8);
  N++;
}

void erase(unsigned value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  // assert(lookup(&tree,key,8,0,8));
  // erase(tree,&tree,key,8,0,8);
  erase(tree,&tree,key,8,0,8);
  // assert(!lookup(&tree,key,8,0,8));
  N--;
}

unsigned lower_bound(unsigned value) {
  uint64_t value64 = value;
  uint8_t key[8];
  loadKey(value64, key);

  // #ifndef EAGER
  // static int nq = 0;
  // if (nq++ == 10000000) {
  //   fprintf(stderr, "flush all inserts ... ");
  //   rec_flush_pending(tree,8,0,8);
  //   fprintf(stderr, "done\n");
  // }
  // #endif

  Node* leaf=lower_bound(tree,key,8,0,8);
  unsigned ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stdout, "%d (%d)\n", ret, value); fflush(stdout);
  }
  return ret;
}

unsigned select(unsigned a, unsigned b) {
  return lower_bound(a) + lower_bound(b);
}

void results(Statistics &s) {
  s.N = N;
}
