#include "assert.h"
#include "art.h"

Node* tree = NULL;

void lookup(uint64_t value64) {
  uint8_t key[8];
  loadKey(value64, key);
  Node *leaf = ::lookup(tree,key,8,0,8);
  assert(leaf);
  assert(isLeaf(leaf));
  assert(getLeafValue(leaf) == value64);
}

uint64_t lbp(uint64_t value64) {
  fprintf(stderr, "lbp %llu\n", value64);
  uint8_t key[8];
  loadKey(value64, key);
  Node *leaf = ::lower_bound_prev(tree,key,8,0,8);
  fprintf(stderr, "lbp2 %p\n", leaf);
  if (!leaf) return -1;
  assert(leaf);
  assert(isLeaf(leaf));
  fprintf(stderr, "%u\n", getLeafValue(leaf));
  return getLeafValue(leaf);
}

uint64_t lb(uint64_t value64) {
  fprintf(stderr, "lb %llu\n", value64);
  uint8_t key[8];
  loadKey(value64, key);
  Node *leaf = ::lower_bound(tree,key,8,0,8);
  if (!leaf) return -1;
  assert(leaf);
  assert(isLeaf(leaf));
  return getLeafValue(leaf);
}

uint64_t get_root(uint64_t value) {
  uint64_t ret = lbp(value);
  if (ret != -1) return ret;
  return lb(value);
}

void insert(uint64_t value64) {
  fprintf(stderr, "insert %llu\n", value64);
  uint8_t key[8];
  loadKey(value64, key);
  ::insert(tree,&tree,key,0,value64,8);
}


void erase(uint64_t value64) {
  uint8_t key[8];
  loadKey(value64, key);
  ::erase(tree,&tree,key,8,0,8);
  fprintf(stderr, "erased %llu, root = %p\n", value64, tree);
  assert(!::lookup(tree,key,8,0,8));
}

int main() {
  fprintf(stderr, "sizeof = %d\n", sizeof(uintptr_t));

  // insert(34234);
  // insert(34235);
  // lookup(34234);
  // erase(34234);
  // erase(34235);
  // fprintf(stderr, "here\n");


  // insert(100);
  // assert(100 == get_root(100));
  // assert(100 == get_root(101));
  // assert(100 == get_root(11));
  // insert(200);
  // assert(100 == get_root(100));
  // assert(100 == get_root(101));
  // assert(100 == get_root(11));
  // assert(100 == get_root(199));
  // assert(200 == get_root(200));
  // assert(200 == get_root(201));
  // assert(200 == get_root(211));
  // assert(200 == get_root(221));
  // // assert(200 == get_root(311));
  // insert(300);
  // assert(100 == get_root(100));
  // assert(100 == get_root(101));
  // assert(100 == get_root(11));
  // assert(100 == get_root(199));
  // assert(200 == get_root(200));
  // assert(200 == get_root(201));
  // assert(200 == get_root(211));
  // assert(200 == get_root(221));
  // assert(200 == get_root(299));
  // assert(300 == get_root(300));
  // assert(300 == get_root(301));


// insert root at 1443471374, b = 1
// remove root at 1443471374, b = 1
// insert root at 555337088, b = 1
// insert root at 904785014, b = 78126
// remove root at 555337088, b = 1


  uint64_t value = (1443471374LL << 29) | 1;
  fprintf(stderr, "dvalue = %llu, root = %p\n", value, tree);
  fprintf(stderr, "value = %llu\n", value);
  insert(value);
  fprintf(stderr, "dvalue = %llu, root = %p\n", value, tree);
  erase(value);
  fprintf(stderr, "value = %llu\n", value);

  value = (555337088LL << 29) | 1;
  insert(value);
  value = (904785014LL << 29) | 78126;
  insert(value);
  erase(value);
}
