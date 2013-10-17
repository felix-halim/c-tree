#include <cstdio>
#include <set>
#include <cassert>
#include <algorithm>

#include "test.h"
#include "crack.h"

using namespace std;

static int arr_cap;

void init(int *a, int N_) {
  ci.clear();
  msize = 0;
  N = N_;
  // marr = new int[N*2];
  arr_cap = max(450000000, N * 2);
  arr = new int[arr_cap];     // for updates expansion
  for (int i = 0; i < N; i++)
    arr[i] = a[i];  // copy all
}

void insert(int value) {
  if (pdel.count(value)){
    pdel.erase(pdel.lower_bound(value));    // don't insert if exists in pdel
  } else {
    pins.insert(value);
  }
}

void erase(int value) {
  if (pins.count(value)){
    pins.erase(pins.lower_bound(value));    // don't delete if exists in pins
  } else {
    pdel.insert(value);
  }
}

int view_query(int a, int b);
int count_query(int a, int b);

int lower_bound(int value) {
  return view_query(value, value + 1) ? value : 0;
}

int select(int a, int b) {
  return view_query(a, b) ? a : 0;
}

int count(int a, int b) {
  assert(N + 1000000 < arr_cap);
  return count_query(a, b);
}

void results(Statistics &s) {
}

Random rr;

void naive_random_crack(){
  value_type x = arr[rr.nextInt(N)];
  int L,R; find_piece(ci, N, x,L,R);
  add_crack(ci, N, x, partition(arr, x,L,R));
}

int view_query(int a, int b){
  #ifdef RANDOM_CRACK_PER_QUERY
    for (int i=0; i<RANDOM_CRACK_PER_QUERY; i++)
      naive_random_crack();
  #endif
  
  #ifdef RANDOM_CRACK_EVERY_NTH_QUERY
    static int nth = 0;
    if (++nth % RANDOM_CRACK_EVERY_NTH_QUERY == 0)
      naive_random_crack();
  #endif

  int cnt = crack(a,b);
  return cnt;
}

int count_query(int a, int b){
  view_query(a,b);

  int L=0, cnt=0;
  ci_iter it1, it2;
  if (ci.count(a)){
    it1 = ci.lower_bound(a);
    L = it1->second.pos;
  } else {
    int L1,R1;
    it1 = find_piece(ci, N, a, L1, R1);
    for (int i=L1; i<R1; i++)
      if (arr[i] >= a && arr[i] < b) cnt++;
    L = R1;
    it1++;
  }
  assert(it1 != ci.end());
  it1++;

  if (ci.count(b)){
    it2 = ci.lower_bound(b);
  } else {
    int L2, R2;
    it2 = find_piece(ci, N, b, L2, R2);
    for (int i=L2; i<R2; i++)
      if (arr[i] >= a && arr[i] < b) cnt++;
    if (it1 == it2) return cnt;
    if (it2 != ci.begin()) it2--;
  }

  while (it1 != ci.end()) {
    if (it1->first > it2->first) break;
    cnt += it1->second.prev_pos() - L;
    L = it1->second.pos;
    if (it1 == it2) break;
    it1++;
  }
  return cnt;
}
