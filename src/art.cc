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

void init(unsigned *arr, unsigned N) {
  for (int i = 0; i < (int) N; i++) {
    insert(arr[i]);
  }
  // fprintf(stderr, "hash = %llu\n", (unsigned long long) hash_tree(tree));
  // art_debug = 1;
}

void insert(unsigned value64) {
  uint8_t key[8];
  loadKey(value64, key);
  if (!lookup(tree,key,8,0,8)) {
    insert(tree,&tree,key,0,value64,8);
  } else {
    fprintf(stderr, "D"); // Duplicate.
  }
}

void erase(unsigned value) {
  uint8_t key[8];
  uint64_t value64 = value;
  loadKey(value64, key);
  // assert(lookup(tree,key,8,0,8));
  erase(tree,&tree,key,8,0,8);
  // assert(!lookup(tree,key,8,0,8));
}

unsigned lower_bound(unsigned value64) {
  uint8_t key[8];
  loadKey(value64, key);
   
  Node* leaf=lower_bound(tree,key,8,0,8);
  unsigned ret = 0;
  if (isLeaf(leaf)) {
    ret = getLeafValue(leaf);
    // fprintf(stdout, "%lld (%lld)\n", ret, value);
  }
  return ret;
}

unsigned select(unsigned a, unsigned b) {
  return lower_bound(a) + lower_bound(b);
}

void results(Statistics &s) {
  s.n_bytes = 0;
  s.n_slack_int = 0;
  s.n_internal = 0,
  s.n_leaf = 0;
  s.art_n4 = 0;
  s.art_n16 = 0;
  s.art_n48 = 0;
  s.art_n256 = 0;

  art_visit(tree, [&](Node *n) {
    if (isLeaf(n)) {
      s.n_leaf++;
      uintptr_t v = getData((uintptr_t) n);
      assert(v);
      assert(!isPointer(v));
      s.n_bytes += sizeof(uintptr_t);
    } else {
      s.n_internal++;
      switch (n->type) {
        case NodeType4: {
           Node4* node = static_cast<Node4*>(n);
           s.n_slack_int += 4 - node->count;
           s.art_n4++;
           s.n_bytes += sizeof(NodeType4);
           break;
        }
        case NodeType16: {
           Node16* node=static_cast<Node16*>(n);
           s.n_slack_int += 16 - node->count;
           s.art_n16++;
           s.n_bytes += sizeof(NodeType16);
           break;
        }
        case NodeType48: {
           Node48* node=static_cast<Node48*>(n);
           s.n_slack_int += 48 - node->count;
           s.art_n48++;
           s.n_bytes += sizeof(NodeType48);
           break;
        }
        case NodeType256: {
           Node256* node=static_cast<Node256*>(n);
           s.n_slack_int += 256 - node->count;
           s.art_n256++;
           s.n_bytes += sizeof(NodeType256);
           break;
        }
      }
    }
  });
}
