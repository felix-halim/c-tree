#include <cstdio>
#include <cassert>
#include <algorithm>
#include "art.h"
#include "test.h"

using namespace std;

Node* tree=NULL;

void init(int *arr, int N) {
  for (int i = 0; i < N; i++) {
    uint8_t key[8];
    loadKey(arr[i], key);
    insert(tree,&tree,key,0,arr[i],8);
  }
}

void insert(int value) {
  // art_insert(&t, (char*) &value, 4, (void*)&value);
}

void erase(int value) {
  // art_delete(&t, (char*) &value, 4);
}

int query(int value) {
  uint8_t key[8];loadKey(value, key);
  Node* leaf=lookup(tree,key,8,0,8);
  if (isLeaf(leaf)) return getLeafValue(leaf);
  return 0;
}

void results(Statistics &s) {
}
