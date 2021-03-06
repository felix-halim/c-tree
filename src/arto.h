/*
  Cracked Adaptive Radix Tree
  Felix Halim, 2013
  felix.halim@gmail.com
 */

#include <stdlib.h>    // malloc, free
#include <string.h>    // memset, memcpy
#include <stdint.h>    // integer types
#include <emmintrin.h> // x86 SSE intrinsics
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>  // gettime
#include <algorithm>   // std::random_shuffle

// Constants for the node types
static const int8_t NodeType4=0;
static const int8_t NodeType16=1;
static const int8_t NodeType48=2;
static const int8_t NodeType256=3;

// The maximum prefix length for compressed paths stored in the
// header, if the path is longer it is loaded from the database on
// demand
static const unsigned maxPrefixLength=9;

int art_debug = 0;

#ifdef DNDEBUG
   #define ART_DEBUG(...)
#else
   #define ART_DEBUG(...) { if (art_debug) fprintf(stderr, __VA_ARGS__); }
#endif

#define BSIZE 1024

struct Bucket {
   int D[BSIZE], N;
   Bucket *next;

   Bucket() { N = 0; next = 0; };
};

// Shared header of all inner nodes
struct Node {
   // length of the compressed path (prefix)
   uint32_t prefixLength;
   // number of non-null children
   int16_t count;
   // node type
   int8_t type;
   // compressed path (prefix)
   uint8_t prefix[maxPrefixLength];

   Bucket *next, *tail;

   Node(int8_t type) : prefixLength(0),count(0),type(type),next(0),tail(0) {}
};

// Node with up to 4 children
struct Node4 : Node {
   uint8_t key[4];
   Node* child[4];

   Node4() : Node(NodeType4) {
      memset(key,0,sizeof(key));
      memset(child,0,sizeof(child));
   }
};

// Node with up to 16 children
struct Node16 : Node {
   uint8_t key[16];
   Node* child[16];

   Node16() : Node(NodeType16) {
      memset(key,0,sizeof(key));
      memset(child,0,sizeof(child));
   }
};

static const uint8_t emptyMarker=48;

// Node with up to 48 children
struct Node48 : Node {
   uint8_t childIndex[256];
   Node* child[48];

   Node48() : Node(NodeType48) {
      memset(childIndex,emptyMarker,sizeof(childIndex));
      memset(child,0,sizeof(child));
   }
};

// Node with up to 256 children
struct Node256 : Node {
   Node* child[256];

   Node256() : Node(NodeType256) {
      memset(child,0,sizeof(child));
   }
};

inline Node* makeLeaf(uintptr_t tid) {
   // Create a pseudo-leaf
   return reinterpret_cast<Node*>((tid<<1)|1);
}

inline uintptr_t getLeafValue(Node* node);

inline bool isLeaf(Node* node) {
   // Is the node a leaf?
   return reinterpret_cast<uintptr_t>(node)&1;
}

uint8_t flipSign(uint8_t keyByte) {
   // Flip the sign bit, enables signed SSE comparison of unsigned values, used by Node16
   return keyByte^128;
}

void loadKey(uintptr_t tid,uint8_t key[]) {
   // Store the key of the tuple into the key vector
   // Implementation is database specific
   reinterpret_cast<uint64_t*>(key)[0]=__builtin_bswap64(tid);
}

// This address is used to communicate that search failed
Node* nullNode=NULL;

static inline unsigned ctz(uint16_t x) {
   // Count trailing zeros, only defined for x>0
#ifdef __GNUC__
   return __builtin_ctz(x);
#else
   // Adapted from Hacker's Delight
   unsigned n=1;
   if ((x&0xFF)==0) {n+=8; x=x>>8;}
   if ((x&0x0F)==0) {n+=4; x=x>>4;}
   if ((x&0x03)==0) {n+=2; x=x>>2;}
   return n-(x&1);
#endif
}

Node** findChild(Node* n,uint8_t keyByte) {
   // Find the next child for the keyByte
   // ART_DEBUG("FC = %p, %d\n", n, n->type);
   switch (n->type) {
      case NodeType4: {
         // ART_DEBUG("FC4 = %u\n", keyByte);
         Node4* node=static_cast<Node4*>(n);
         for (int i=0;i<node->count;i++)
            if (node->key[i]==keyByte)
               return &node->child[i];
         return &nullNode;
      }
      case NodeType16: {
         // ART_DEBUG("FC16 = %u\n", keyByte);
         Node16* node=static_cast<Node16*>(n);
         __m128i cmp=_mm_cmpeq_epi8(_mm_set1_epi8(flipSign(keyByte)),_mm_loadu_si128(reinterpret_cast<__m128i*>(node->key)));
         unsigned bitfield=_mm_movemask_epi8(cmp)&((1<<node->count)-1);
         if (bitfield)
            return &node->child[ctz(bitfield)]; else
            return &nullNode;
      }
      case NodeType48: {
         // ART_DEBUG("FC48 = %u\n", keyByte);
         Node48* node=static_cast<Node48*>(n);
         if (node->childIndex[keyByte]!=emptyMarker)
            return &node->child[node->childIndex[keyByte]]; else
            return &nullNode;
      }
      case NodeType256: {
         // ART_DEBUG("FC256 = %u\n", keyByte);
         Node256* node=static_cast<Node256*>(n);
         return &(node->child[keyByte]);
      }
   }
   throw; // Unreachable
}

Node* minimum(Node* node) {
   // Find the leaf with smallest key
   if (!node)
      return NULL;

   if (isLeaf(node))
      return node;

   switch (node->type) {
      case NodeType4: {
         Node4* n=static_cast<Node4*>(node);
         return minimum(n->child[0]);
      }
      case NodeType16: {
         Node16* n=static_cast<Node16*>(node);
         return minimum(n->child[0]);
      }
      case NodeType48: {
         Node48* n=static_cast<Node48*>(node);
         unsigned pos=0;
         while (n->childIndex[pos]==emptyMarker)
            pos++;
         return minimum(n->child[n->childIndex[pos]]);
      }
      case NodeType256: {
         Node256* n=static_cast<Node256*>(node);
         unsigned pos=0;
         while (!n->child[pos])
            pos++;
         return minimum(n->child[pos]);
      }
   }
   throw; // Unreachable
}

Node* maximum(Node* node) {
   // Find the leaf with largest key
   if (!node)
      return NULL;

   if (isLeaf(node))
      return node;

   switch (node->type) {
      case NodeType4: {
         Node4* n=static_cast<Node4*>(node);
         return maximum(n->child[n->count-1]);
      }
      case NodeType16: {
         Node16* n=static_cast<Node16*>(node);
         return maximum(n->child[n->count-1]);
      }
      case NodeType48: {
         Node48* n=static_cast<Node48*>(node);
         unsigned pos=255;
         while (n->childIndex[pos]==emptyMarker)
            pos--;
         return maximum(n->child[n->childIndex[pos]]);
      }
      case NodeType256: {
         Node256* n=static_cast<Node256*>(node);
         unsigned pos=255;
         while (!n->child[pos])
            pos--;
         return maximum(n->child[pos]);
      }
   }
   throw; // Unreachable
}

bool leafMatches(Node* leaf,uint8_t key[],unsigned keyLength,unsigned depth,unsigned maxKeyLength) {
   // Check if the key of the leaf is equal to the searched key
   if (depth!=keyLength) {
      uint8_t leafKey[maxKeyLength];
      loadKey(getLeafValue(leaf),leafKey);
      for (unsigned i=depth;i<keyLength;i++)
         if (leafKey[i]!=key[i]) {
            ART_DEBUG("LEAF NOT MATCH\n");
            return false;
         }
   }
   return true;
}

unsigned prefixMismatch(Node* node,uint8_t key[],unsigned depth,unsigned maxKeyLength) {
   // Compare the key with the prefix of the node, return the number matching bytes
   unsigned pos;
   if (node->prefixLength>maxPrefixLength) {
      for (pos=0;pos<maxPrefixLength;pos++)
         if (key[depth+pos]!=node->prefix[pos])
            return pos;
      uint8_t minKey[maxKeyLength];
      loadKey(getLeafValue(minimum(node)),minKey);
      for (;pos<node->prefixLength;pos++)
         if (key[depth+pos]!=minKey[depth+pos])
            return pos;
   } else {
      for (pos=0;pos<node->prefixLength;pos++)
         if (key[depth+pos]!=node->prefix[pos])
            return pos;
   }
   return pos;
}

Node* lookup(Node* node,uint8_t key[],unsigned keyLength,unsigned depth,unsigned maxKeyLength) {
   // Find the node with a matching key, optimistic version
   if (!node) return NULL;

   bool skippedPrefix=false; // Did we optimistically skip some prefix without checking it?

   while (node!=NULL) {
      if (isLeaf(node)) {
         if (!skippedPrefix&&depth==keyLength) // No check required
            return node;

         if (depth!=keyLength) {
            // Check leaf
            uint8_t leafKey[maxKeyLength];
            loadKey(getLeafValue(node),leafKey);
            for (unsigned i=(skippedPrefix?0:depth);i<keyLength;i++)
               if (leafKey[i]!=key[i])
                  return NULL;
         }
         return node;
      }

      if (node->prefixLength) {
         if (node->prefixLength<maxPrefixLength) {
            for (unsigned pos=0;pos<node->prefixLength;pos++)
               if (key[depth+pos]!=node->prefix[pos])
                  return NULL;
         } else
            skippedPrefix=true;
         depth+=node->prefixLength;
      }

      node=*findChild(node,key[depth]);
      depth++;
   }

   return NULL;
}

void flush_inserts(Node *&node, int depth, int maxKeyLength);

Node* lower_bound(Node *node, Node **nodeRef, uint8_t key[], unsigned keyLength, unsigned depth, unsigned maxKeyLength, bool skippedPrefix=false, bool bigger = false) {
   ART_DEBUG("lower_boundx depth = %u, %p, bigger = %d\n", depth, node, bigger);
   if (!node) return NULL;

   if (isLeaf(node)) {
      if (!skippedPrefix && depth == keyLength) // No check required
         return node;

      if (depth && depth != keyLength && !bigger) {
         uint8_t leafKey[maxKeyLength];
         loadKey(getLeafValue(node),leafKey);
         #ifndef NDEBUG
            for (unsigned i = 0; i < keyLength; i++) {
               ART_DEBUG("i = %d, %u %u\n", i, leafKey[i], key[i]);
            }
         #endif
         for (unsigned i=(skippedPrefix?0:depth);i<keyLength;i++) {
            ART_DEBUG("i => %u, %u %u\n", i, leafKey[i], key[i]);
            if (leafKey[i] < key[i]) {
               ART_DEBUG("ARRRR!!! \n");
               return NULL;
            } else if (leafKey[i] > key[i]) {
               ART_DEBUG("SHORTCUT\n");
               break;
            }
         }
      }
      ART_DEBUG("YES\n");
      return node;
   }

   flush_inserts(*nodeRef, depth, maxKeyLength);
   node = *nodeRef;

   // ART_DEBUG("prefixLength = %u\n", node->prefixLength);
   if (node->prefixLength) {
      if (node->prefixLength < maxPrefixLength && !bigger) {
         for (unsigned pos=0; pos < node->prefixLength; pos++)
            if (key[depth+pos] > node->prefix[pos]) {
               ART_DEBUG("ARRRRGGGGGHhHHHHAaAJSKHJAKSJDHKAJSHDKASJDHAKSJHDAKJSDHAKSJDHAKSJDHASKDJHAS\n");
               return NULL;
            } else if (key[depth+pos] < node->prefix[pos]) {
               bigger = 1;
               break;
            }
      } else
         skippedPrefix=true;
      depth+=node->prefixLength;
   }

   Node *n = node;
   int keyByte = key[depth++];

   switch (n->type) {
      case NodeType4: {
            ART_DEBUG("LoweBound4, count = %d\n", node->count);
            Node4* node = static_cast<Node4*>(n);
            for (int i = 0; i < node->count; i++) {
               Node *c = node->child[i];
               ART_DEBUG("i = %d, key4 = %d >= %d, %lu\n", i, node->key[i], keyByte, isLeaf(c) ? getLeafValue(c) : 0);
               if (node->key[i] >= keyByte || bigger) {
                  // ART_DEBUG("got it\n");
                  Node *ret = lower_bound(c, &node->child[i], key, keyLength, depth, maxKeyLength, skippedPrefix, bigger || node->key[i] > keyByte);
                  if (ret) return ret;
               }
            }
         }
         break;

      case NodeType16: {
            ART_DEBUG("LoweBound16 %d\n", node->count);
            int pos = 0;
            Node16* node = static_cast<Node16*>(n);
            if (!bigger) {
               __m128i cmp=_mm_cmplt_epi8(_mm_set1_epi8(flipSign(keyByte)),_mm_loadu_si128(reinterpret_cast<__m128i*>(node->key)));
               uint16_t bitfield=_mm_movemask_epi8(cmp)&(0xFFFF>>(16-node->count));
               pos = bitfield ? ctz(bitfield) : node->count;
               while (pos > 0 && flipSign(node->key[pos - 1]) >= keyByte) pos--;
            }
            ART_DEBUG("pos = %d\n", pos);
            while (pos < node->count) {
               Node *ret = lower_bound(node->child[pos], &node->child[pos], key, keyLength, depth, maxKeyLength, skippedPrefix, bigger || flipSign(node->key[pos]) > keyByte);
               if (ret) return ret;
               pos++;
            }
         }
         break;

      case NodeType48: {
            ART_DEBUG("LoweBound48\n");
            Node48* node=static_cast<Node48*>(n);
            if (bigger) keyByte = 0;
            while (keyByte < 256) {
               if (node->childIndex[keyByte] != emptyMarker) {
                  int i = node->childIndex[keyByte];
                  Node *ret = lower_bound(node->child[i], &node->child[i], key, keyLength, depth, maxKeyLength, skippedPrefix, bigger);
                  if (ret) return ret;
               }
               keyByte++;
               bigger = 1;
            }
         }
         break;

      case NodeType256: {
            ART_DEBUG("LoweBound256 %d\n", node->count);
            Node256* node=static_cast<Node256*>(n);
            if (bigger) keyByte = 0;
            while (keyByte < 256) {
               ART_DEBUG("LoweBound256 keyByte = %d, bigger = %d\n", keyByte, bigger);
               if (node->child[keyByte]) {
                  Node *ret = lower_bound(node->child[keyByte], &node->child[keyByte], key, keyLength, depth, maxKeyLength, skippedPrefix, bigger);
                  if (ret) return ret;
               }
               keyByte++;
               bigger = 1;
            }
         }
         break;

      default: assert(0);
   }

   ART_DEBUG("NOOOOOOOOOOO\n");
   return NULL;
}

Node* lower_bound_prev(Node* node, uint8_t key[], unsigned keyLength, unsigned depth, unsigned maxKeyLength, bool skippedPrefix=false, bool is_less = false) {
   // ART_DEBUG("lower_prev1 %u, leaf = %d\n", depth, isLeaf(node));
   if (!node) return NULL;

   if (isLeaf(node)) {
      if (!skippedPrefix && depth == keyLength) // No check required
         return node;

      if (depth && depth != keyLength && !is_less) {
         uint8_t leafKey[maxKeyLength];
         loadKey(getLeafValue(node),leafKey);
         for (unsigned j = 0; j < keyLength; j++) {
            ART_DEBUG("%u: %u %u\n", j, leafKey[j], key[j]);
         }
         ART_DEBUG("THE VALUE = %lu\n", getLeafValue(node) >> 30);
         for (unsigned i=(skippedPrefix?0:depth);i<keyLength;i++) {
            if (leafKey[i] > key[i]) {
               ART_DEBUG("GAGAL %u %u\n", leafKey[i], key[i]);
               return NULL;
            } else if (leafKey[i] < key[i]) {
               ART_DEBUG("WOHOOOO");
               break;
            }
         }
      }
      ART_DEBUG("VALUE = %lu, is_less = %d, %lu\n", getLeafValue(node) >> 30, is_less, sizeof(uintptr_t));
      return node;
   }

   // ART_DEBUG("lower_prev1 plen = %u\n", node->prefixLength);
   if (node->prefixLength) {
      if (node->prefixLength < maxPrefixLength && !is_less) {
         for (unsigned pos=0; pos < node->prefixLength; pos++)
            if (key[depth+pos] < node->prefix[pos]) {
               ART_DEBUG("GAGAL2 %u %u, pos = %u\n", key[depth+pos], node->prefix[pos], pos);
               return NULL;
            } else if (key[depth+pos] > node->prefix[pos]) {
               is_less = 1;
               break;
            }
      } else {
         skippedPrefix=true;
      }
      depth+=node->prefixLength;
   }

   Node *n = node;
   int keyByte = key[depth++];

   switch (n->type) {
      case NodeType4: {
            ART_DEBUG("Node4P %d\n", node->count);
            Node4* node = static_cast<Node4*>(n);
            for (int i = node->count - 1; i >= 0; i--) {
               if (isLeaf(node->child[i])) {
                  ART_DEBUG("child %d, value = %lu\n", i, getLeafValue(node->child[i]) >> 30);
               }
               if (node->key[i] <= keyByte || is_less) {
                  Node *ret = lower_bound_prev(node->child[i], key, keyLength, depth, maxKeyLength, skippedPrefix, is_less || node->key[i] < keyByte);
                  if (ret) return ret;
               }
            }
         }
         break;

      case NodeType16: {
            ART_DEBUG("Node16P count = %d, is_less = %d\n", node->count, is_less);
            Node16* node = static_cast<Node16*>(n);
            int pos = node->count - 1;
            if (!is_less) {
               // ART_DEBUG("Node16 is_less %d\n", node->count);
               __m128i cmp=_mm_cmplt_epi8(_mm_set1_epi8(flipSign(keyByte)),_mm_loadu_si128(reinterpret_cast<__m128i*>(node->key)));
               uint16_t bitfield=_mm_movemask_epi8(cmp)&(0xFFFF>>(16-node->count));
               pos = bitfield ? ctz(bitfield) : node->count;
               while (pos >= node->count || flipSign(node->key[pos]) > keyByte) pos--;
               while (pos + 1 < node->count && flipSign(node->key[pos + 1]) <= keyByte) pos++;
            }
            ART_DEBUG("pos = %d\n", pos);
            assert(pos < node->count);
            while (pos >= 0) {
               // ART_DEBUG("Node16 pos %d, %p, %u %u\n", pos, node->child[pos], flipSign(node->key[pos]), keyByte);
               if (is_less || flipSign(node->key[pos]) <= keyByte) {
                  Node *ret = lower_bound_prev(node->child[pos], key, keyLength, depth, maxKeyLength, skippedPrefix, is_less || flipSign(node->key[pos]) < keyByte);
                  if (ret) return ret;
               } else {
                  // ART_DEBUG("FAIL\n");
               }
               pos--;
            }
         }
         break;

      case NodeType48: {
            ART_DEBUG("Node48P %d is_less = %d, keyByte = %d\n", node->count, is_less, keyByte);
            Node48* node=static_cast<Node48*>(n);
            if (is_less) keyByte = 255;
            while (keyByte >= 0) {
               if (node->childIndex[keyByte] != emptyMarker) {
                  Node *c = node->child[node->childIndex[keyByte]];
                  ART_DEBUG("N48C %d = %lld\n", keyByte, isLeaf(c) ? getLeafValue(c) >> 30 : -1LL);
                  Node *ret = lower_bound_prev(c, key, keyLength, depth, maxKeyLength, skippedPrefix, is_less);
                  if (ret) return ret;
               }
               keyByte--;
               is_less = 1;
            }
         }
         break;

      case NodeType256: {
            ART_DEBUG("Node256P %d\n", node->count);
            Node256* node=static_cast<Node256*>(n);
            if (is_less) keyByte = 255;
            while (keyByte >= 0) {
               if (node->child[keyByte]) {
                  Node *ret = lower_bound_prev(node->child[keyByte], key, keyLength, depth, maxKeyLength, skippedPrefix, is_less);
                  if (ret) return ret;
               }
               keyByte--;
               is_less = 1;
            }
         }
         break;

      default: assert(0);
   }
   return NULL;
}

Node* lookupPessimistic(Node* node,uint8_t key[],unsigned keyLength,unsigned depth,unsigned maxKeyLength) {
   // Find the node with a matching key, alternative pessimistic version

   while (node!=NULL) {
      if (isLeaf(node)) {
         if (leafMatches(node,key,keyLength,depth,maxKeyLength))
            return node;
         return NULL;
      }

      if (prefixMismatch(node,key,depth,maxKeyLength)!=node->prefixLength)
         return NULL; else
         depth+=node->prefixLength;

      node=*findChild(node,key[depth]);
      depth++;
   }

   return NULL;
}

// Forward references
void insertNode4(Node4* node,Node** nodeRef,uint8_t keyByte,Node* child);
void insertNode16(Node16 *&node,uint8_t keyByte,Node* child);
void insertNode48(Node48 *&node,uint8_t keyByte,Node* child);
void insertNode256(Node256 *&node,uint8_t keyByte,Node* child);

unsigned min(unsigned a,unsigned b) {
   // Helper function
   return (a<b)?a:b;
}

void copyPrefix(Node* src,Node* dst) {
   // Helper function that copies the prefix from the source to the destination node
   dst->prefixLength=src->prefixLength;
   memcpy(dst->prefix,src->prefix,min(src->prefixLength,maxPrefixLength));
}

void insert(Node *node, Node **nodeRef,uint8_t key[],unsigned depth,uintptr_t value,unsigned maxKeyLength, bool eager);
void flush_inserts(Node* node,Node** nodeRef,uint8_t key[],unsigned depth,uintptr_t value,unsigned maxKeyLength, bool eager) {
   // ART_DEBUG("continue insert %d\n", depth);

   // Handle prefix of inner node
   if (node->prefixLength) {
      unsigned mismatchPos=prefixMismatch(node,key,depth,maxKeyLength);
      if (mismatchPos!=node->prefixLength) {
         // Prefix differs, create new node
         ART_DEBUG("Prefix Differs\n");
         Node4* newNode=new Node4();
         *nodeRef=newNode;
         newNode->next = node->next;
         newNode->tail = node->tail;
         node->next = node->tail = 0;
         newNode->prefixLength=mismatchPos;
         memcpy(newNode->prefix,node->prefix,min(mismatchPos,maxPrefixLength));
         // Break up prefix
         if (node->prefixLength<maxPrefixLength) {
            ART_DEBUG("Break Prefix\n");
            insertNode4(newNode,nodeRef,node->prefix[mismatchPos],node);
            node->prefixLength-=(mismatchPos+1);
            memmove(node->prefix,node->prefix+mismatchPos+1,min(node->prefixLength,maxPrefixLength));
         } else {
            ART_DEBUG("Break Prefix2\n");
            node->prefixLength-=(mismatchPos+1);
            uint8_t minKey[maxKeyLength];
            loadKey(getLeafValue(minimum(node)),minKey);
            insertNode4(newNode,nodeRef,minKey[depth+mismatchPos],node);
            memmove(node->prefix,minKey+depth+mismatchPos+1,min(node->prefixLength,maxPrefixLength));
         }
         insertNode4(newNode,nodeRef,key[depth+mismatchPos],makeLeaf(value));
         return;
      }
      depth+=node->prefixLength;
   }

   assert(node == *nodeRef);
   // Recurse
   Node** child=findChild(node,key[depth]);
   if (*child) {
      ART_DEBUG("Recurse %d from %p to %p; ", depth, node, *child);
      insert(*child,child,key,depth+1,value,maxKeyLength,eager);
      return;
   }

   // Insert leaf into inner node
   Node* newNode=makeLeaf(value);
   // ART_DEBUG("Insert Inner\n");
   switch (node->type) {
      case NodeType4: insertNode4(static_cast<Node4*>(node),nodeRef,key[depth],newNode); break;
      case NodeType16: insertNode16((Node16*&) *nodeRef,key[depth],newNode); break;
      case NodeType48: insertNode48((Node48*&) *nodeRef,key[depth],newNode); break;
      case NodeType256: insertNode256((Node256*&) *nodeRef,key[depth],newNode); break;
   }
   ART_DEBUG("done inner \n");
}

void flush_inserts(Node *&node, int depth, int maxKeyLength) {
   if (!node->next) return;

   ART_DEBUG("flushing_inserts %d, next = %p\n", depth, node->next);
   while (node->next) {
      uint8_t key[8];
      ART_DEBUG("flushing bucket = %p\n", node->next);
      bool eager = node->next == node->tail;
      for (int i = 0; i < node->next->N; i++) {
         loadKey(node->next->D[i], key);
         // ART_DEBUG("LOADED KEY %d / %d\n", i, node->next->N);
         flush_inserts(node, &node, key, depth, node->next->D[i], maxKeyLength, eager);
         assert(node->next);
         // ART_DEBUG("f %d / %d\n", i, node->next->N);
      }
      Bucket *next = node->next->next;
      ART_DEBUG("deleting %p, next = %p, host = %p\n", node->next, next, node);
      delete node->next;
      // ART_DEBUG("deleting done, host = %p\n", *nodeRef);
      node->next = next;
   }
   node->tail = 0;
   ART_DEBUG("flushing_inserts done\n");
}

void insert(Node *node, Node **nodeRef,uint8_t key[],unsigned depth,uintptr_t value,unsigned maxKeyLength, bool eager) {
   // Insert the leaf value into the tree
   // ART_DEBUG("Insert depth = %d, %u, node = %p\n", depth, depth < maxKeyLength ? key[depth] : 0, node);

   if (node==NULL) {
      ART_DEBUG("null leaf %d %lu\n", depth, value);
      *nodeRef=makeLeaf(value);
      return;
   }

   if (isLeaf(node)) {
      ART_DEBUG("insert leaf %d %lu\n", depth, value);
      // Replace leaf with Node4 and store both leaves in it
      uint8_t existingKey[maxKeyLength];
      loadKey(getLeafValue(node),existingKey);
      unsigned newPrefixLength=0;
      while (existingKey[depth+newPrefixLength]==key[depth+newPrefixLength])
         newPrefixLength++;

      Node4* newNode=new Node4();
      newNode->prefixLength=newPrefixLength;
      memcpy(newNode->prefix,key+depth,min(newPrefixLength,maxPrefixLength));
      *nodeRef=newNode;

      insertNode4(newNode,nodeRef,existingKey[depth+newPrefixLength],node);
      insertNode4(newNode,nodeRef,key[depth+newPrefixLength],makeLeaf(value));
      return;
   }

   if (eager) {
      flush_inserts(node, nodeRef, key, depth, value, maxKeyLength, eager);
   } else {
      // Buffer inserts
      if (!node->next) {
         node->next = node->tail = new Bucket();
         ART_DEBUG("NEW NEXT\n");
         assert(!node->tail->next);
      }
      assert(node->tail);
      if (node->tail->N == BSIZE) {
         assert(!node->tail->next);
         ART_DEBUG("HOST %p -> ... -> %p -> ", node, node->tail);
         node->tail->next = new Bucket();
         node->tail = node->tail->next;
         ART_DEBUG("%p\n", node->tail);
      }
      assert(node->tail->N < BSIZE);
      node->tail->D[node->tail->N++] = value;
      ART_DEBUG("append %lu to %p\n", value, node->tail);
   }
}

void insertNode4(Node4* node,Node** nodeRef,uint8_t keyByte,Node* child) {
   // Insert leaf into inner node
   if (node->count<4) {
      // Insert element
      ART_DEBUG("insert4 %d\n", node->count);
      int pos;
      for (pos=0;(pos<node->count)&&(node->key[pos]<keyByte);pos++);
      memmove(node->key+pos+1,node->key+pos,node->count-pos);
      memmove(node->child+pos+1,node->child+pos,(node->count-pos)*sizeof(uintptr_t));
      node->key[pos]=keyByte;
      node->child[pos]=child;
      node->count++;
      ART_DEBUG("insert4d %d\n", node->count);
   } else {
      // Grow to Node16
      ART_DEBUG("insert4to16 %d\n", node->count);
      Node16* newNode=new Node16();
      *nodeRef=newNode;
      newNode->count=4;
      copyPrefix(node,newNode);
      for (unsigned i=0;i<4;i++)
         newNode->key[i]=flipSign(node->key[i]);
      memcpy(newNode->child,node->child,node->count*sizeof(uintptr_t));
      newNode->next = node->next;
      newNode->tail = node->tail;
      delete node;
      return insertNode16((Node16*&) *nodeRef,keyByte,child);
   }
}

void insertNode16(Node16* &node,uint8_t keyByte,Node* child) {
   // Insert leaf into inner node
   if (node->count<16) {
      ART_DEBUG("insert16 %d\n", node->count);
      // Insert element
      uint8_t keyByteFlipped=flipSign(keyByte);
      __m128i cmp=_mm_cmplt_epi8(_mm_set1_epi8(keyByteFlipped),_mm_loadu_si128(reinterpret_cast<__m128i*>(node->key)));
      uint16_t bitfield=_mm_movemask_epi8(cmp)&(0xFFFF>>(16-node->count));
      unsigned pos=bitfield?ctz(bitfield):node->count;
      memmove(node->key+pos+1,node->key+pos,node->count-pos);
      memmove(node->child+pos+1,node->child+pos,(node->count-pos)*sizeof(uintptr_t));
      node->key[pos]=keyByteFlipped;
      node->child[pos]=child;
      node->count++;
   } else {
      // Grow to Node48
      ART_DEBUG("insert16to48 %d\n", node->count);
      Node48* newNode=new Node48();
      memcpy(newNode->child,node->child,node->count*sizeof(uintptr_t));
      for (int i=0;i<node->count;i++)
         newNode->childIndex[flipSign(node->key[i])]=i;
      copyPrefix(node,newNode);
      newNode->count=node->count;
      newNode->next = node->next;
      newNode->tail = node->tail;
      delete node;
      node = (Node16*) newNode;
      return insertNode48((Node48*&) node,keyByte,child);
   }
}

void insertNode48(Node48 *&node,uint8_t keyByte,Node* child) {
   // Insert leaf into inner node
   if (node->count<48) {
      // Insert element
      ART_DEBUG("insert48 %d\n", node->count);
      unsigned pos=node->count;
      if (node->child[pos])
         for (pos=0;node->child[pos]!=NULL;pos++);
      node->child[pos]=child;
      node->childIndex[keyByte]=pos;
      node->count++;
   } else {
      // Grow to Node256
      ART_DEBUG("insert48to256 %d\n", node->count);
      Node256* newNode=new Node256();
      for (unsigned i=0;i<256;i++)
         if (node->childIndex[i]!=48)
            newNode->child[i]=node->child[node->childIndex[i]];
      newNode->count=node->count;
      copyPrefix(node,newNode);
      newNode->next = node->next;
      newNode->tail = node->tail;
      delete node;
      node = (Node48*) newNode;
      return insertNode256((Node256*&) node,keyByte,child);
   }
}

void insertNode256(Node256* &node,uint8_t keyByte,Node* child) {
   // Insert leaf into inner node
   ART_DEBUG("insert256 %d\n", node->count);
   node->count++;
   node->child[keyByte]=child;
}

// Forward references
void eraseNode4(Node4* node,Node** nodeRef,Node** leafPlace);
void eraseNode16(Node16* node,Node** nodeRef,Node** leafPlace);
void eraseNode48(Node48* node,Node** nodeRef,uint8_t keyByte);
void eraseNode256(Node256* node,Node** nodeRef,uint8_t keyByte);

void erase(Node* node,Node** nodeRef,uint8_t key[],unsigned keyLength,unsigned depth,unsigned maxKeyLength) {
   // Delete a leaf from a tree

   ART_DEBUG("ERASE %p %d\n", node, depth);

   if (!node) {
      ART_DEBUG("NO NODE\n");
      return;
   }

   if (isLeaf(node)) {
      // Make sure we have the right leaf
      if (leafMatches(node,key,keyLength,depth,maxKeyLength)) {
         ART_DEBUG("DONE DELETE %p %p\n", node, *nodeRef);
         *nodeRef=NULL;
         ART_DEBUG("DONE DELETE %p %p\n", node, *nodeRef);
      } else {
         ART_DEBUG("DONE LEAF NOT MATCH\n");
      }
      return;
   }

   flush_inserts(*nodeRef, depth, maxKeyLength);
   node = *nodeRef;

   // Handle prefix
   // ART_DEBUG("PREFIX LEN = %d\n", node->prefixLength);
   if (node->prefixLength) {
      if (prefixMismatch(node,key,depth,maxKeyLength)!=node->prefixLength) {
         ART_DEBUG("PREFIX MISMATCH\n");
         return;
      }
      depth+=node->prefixLength;
   }

   // ART_DEBUG("PREFIX LEN2 = %d, depth = %u\n", node->prefixLength, depth);
   Node** child = findChild(node,key[depth]);
   // ART_DEBUG("PREFIX LEN3 = %p\n", child);
   // ART_DEBUG("CHILD = %p\n", *child);
   if (isLeaf(*child)&&leafMatches(*child,key,keyLength,depth,maxKeyLength)) {
      // Leaf found, delete it in inner node
      ART_DEBUG("LEAF FOUND\n");
      switch (node->type) {
         case NodeType4: eraseNode4(static_cast<Node4*>(node),nodeRef,child); break;
         case NodeType16: eraseNode16(static_cast<Node16*>(node),nodeRef,child); break;
         case NodeType48: eraseNode48(static_cast<Node48*>(node),nodeRef,key[depth]); break;
         case NodeType256: eraseNode256(static_cast<Node256*>(node),nodeRef,key[depth]); break;
      }
   } else {
      //Recurse
      erase(*child,child,key,keyLength,depth+1,maxKeyLength);
   }
}

void eraseNode4(Node4* node,Node** nodeRef,Node** leafPlace) {
   // Delete leaf from inner node
   unsigned pos=leafPlace-node->child;
   memmove(node->key+pos,node->key+pos+1,node->count-pos-1);
   memmove(node->child+pos,node->child+pos+1,(node->count-pos-1)*sizeof(uintptr_t));
   node->count--;

   ART_DEBUG("Erase4, count = %d\n", node->count);
   if (node->count==1) {
      // Get rid of one-way node
      Node* child=node->child[0];
      if (!isLeaf(child)) {
         // Concantenate prefixes
         unsigned l1=node->prefixLength;
         if (l1<maxPrefixLength) {
            node->prefix[l1]=node->key[0];
            l1++;
         }
         if (l1<maxPrefixLength) {
            unsigned l2=min(child->prefixLength,maxPrefixLength-l1);
            memcpy(node->prefix+l1,child->prefix,l2);
            l1+=l2;
         }
         // Store concantenated prefix
         memcpy(child->prefix,node->prefix,min(l1,maxPrefixLength));
         child->prefixLength+=node->prefixLength+1;
      }
      *nodeRef=child;
      delete node;
   }
}

void eraseNode16(Node16* node,Node** nodeRef,Node** leafPlace) {
   // Delete leaf from inner node
   unsigned pos=leafPlace-node->child;
   memmove(node->key+pos,node->key+pos+1,node->count-pos-1);
   memmove(node->child+pos,node->child+pos+1,(node->count-pos-1)*sizeof(uintptr_t));
   node->count--;
   ART_DEBUG("Erase16\n");

   if (node->count==3) {
      // Shrink to Node4
      Node4* newNode=new Node4();
      newNode->count=3;
      copyPrefix(node,newNode);
      for (unsigned i=0;i<3;i++)
         newNode->key[i]=flipSign(node->key[i]);
      memcpy(newNode->child,node->child,sizeof(uintptr_t)*3);
      *nodeRef=newNode;
      delete node;
   }
}

void eraseNode48(Node48* node,Node** nodeRef,uint8_t keyByte) {
   // Delete leaf from inner node
   node->child[node->childIndex[keyByte]]=NULL;
   node->childIndex[keyByte]=emptyMarker;
   node->count--;
   ART_DEBUG("Erase48\n");

   if (node->count==12) {
      // Shrink to Node16
      Node16 *newNode=new Node16();
      *nodeRef=newNode;
      copyPrefix(node,newNode);
      for (unsigned b=0;b<256;b++) {
         if (node->childIndex[b]!=emptyMarker) {
            newNode->key[newNode->count]=flipSign(b);
            newNode->child[newNode->count]=node->child[node->childIndex[b]];
            newNode->count++;
         }
      }
      delete node;
   }
}

void eraseNode256(Node256* node,Node** nodeRef,uint8_t keyByte) {
   // Delete leaf from inner node
   node->child[keyByte]=NULL;
   node->count--;
   ART_DEBUG("Erase256\n");

   if (node->count==37) {
      // Shrink to Node48
      Node48 *newNode=new Node48();
      *nodeRef=newNode;
      copyPrefix(node,newNode);
      for (unsigned b=0;b<256;b++) {
         if (node->child[b]) {
            newNode->childIndex[b]=newNode->count;
            newNode->child[newNode->count]=node->child[b];
            newNode->count++;
         }
      }
      delete node;
   }
}
