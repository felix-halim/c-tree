#include <cstdio>
#include <cassert>
#include "art.h"
#include "tester.h"

art_tree t;

void initexp() {
  fprintf(stderr, "init art\n");
  int res = art_tree_init(&t);
  assert(res == 0);
}

int count_cb(void *data, const unsigned char* key, uint32_t key_len, void *val) {
  (*((long long*)data))++;
  return 0;
}

void count_size() {
  long long out = 0;
  art_iter(&t, count_cb, &out);
  fprintf(stderr, "art size = %lld\n", out);
}

void destroyexp() {
  fprintf(stderr, "destroy art\n");
  // int res = art_tree_destroy(&t);
  // assert(res == 0);
}

static unsigned char* get_key(long long val) {
  static unsigned char key[8];
  key[7] = val & 0xFF; val >>= 8;
  key[6] = val & 0xFF; val >>= 8;
  key[5] = val & 0xFF; val >>= 8;
  key[4] = val & 0xFF; val >>= 8;
  key[3] = val & 0xFF; val >>= 8;
  key[2] = val & 0xFF; val >>= 8;
  key[1] = val & 0xFF; val >>= 8;
  key[0] = val & 0xFF;
  return key;
}

static long long decode_key(unsigned char *key) {
  long long k = 0;
  for (int i = 0; i < 8; i++) {
    k = (k << 8) | key[i];
  }
  return k;
}

// op = 1: inserts the value.
void insert(long long value) {
  void *old = art_insert(&t, get_key(value), 8, (void*) value);
  if (old) {
    fprintf(stderr, "D"); // Duplicate.
  }
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  void *res = art_delete(&t, get_key(value), 8);
  return res != NULL;
}


long long sentinel;

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  sentinel = b;
  long long out = 0;
  art_iter_lower(&t, (unsigned char*) &a, 8, count_cb, &out);
  return out;
}


int sum_cb(void *data, const unsigned char* key, uint32_t key_len, void *val) {
  long long z = (long long) val;
  long long &d = *((long long*) data);
  if (z >= sentinel) return 1;
  d += z;
  return 0;
}

int print_cb(void *data, const unsigned char* key, uint32_t key_len, void *val) {
  long long z = (long long) val;
  fprintf(stderr, "%lld\n", z);
  return 0;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  sentinel = b;
  long long out = 0;
  // art_iter(&t, print_cb, &out);
  art_iter_lower(&t, get_key(a), 8, sum_cb, &out);
  return out;
}
