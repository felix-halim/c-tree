#include <cstdio>
#include <cassert>
#include <algorithm>
#include "art.h"
#include "test.h"

using namespace std;

Node* tree=NULL;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++) {
    insert(arr[i]);
  }
}

void insert(int value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  insert(tree,&tree,key,0,value64,8);
}

void erase(int value) {
  // art_delete(&t, (char*) &value, 4);
}

int query(int value) {
  uint64_t value64 = value;
  uint8_t key[8];
  loadKey(value64, key);
  Node* leaf=lower_bound(tree,key,8,0,8);
  int ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stderr, "%d (%d)\n", ret, value);
  }
  return ret;
}

void results(Statistics &s) {
}
