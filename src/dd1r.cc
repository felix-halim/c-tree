#include <cstdio>
#include <set>
#include <cassert>
#include <algorithm>

#include "test.h"
#include "crack.h"

using namespace std;

#define MAX_NCRACK 1
#define CRACK_AT 128

void init(int *a, int N_) {
  ci.clear();
  msize = 0;
  N = N_;
  // marr = new int[N*2];
  arr = new int[N*2];     // for updates expansion
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

int lower_bound(int value) {
  return view_query(value, value + 1) ? value : 0;
}

int select(int a, int b) {
  return view_query(a, b) ? a : 0;
}

void results(Statistics &s) {
}

int ddr_find(value_type v){
  int L,R;
  return targeted_random_crack(ci,v,arr,N,L,R,MAX_NCRACK,CRACK_AT);
}

int view_query(int a, int b){
  merge_ripple(ci, arr, N, pins, pdel, a, b);  // merge qualified updates
  int i2 = ddr_find(b);  // unlimited cracks allowed plus one crack on v2
  int i1 = ddr_find(a);  // unlimited cracks allowed plus one crack on v1
  return i2 - i1;
}

int count_query(int a, int b){
  merge_ripple(ci, arr, N, pins, pdel, a, b);  // merge qualified updates
  int i2 = ddr_find(b);  // unlimited cracks allowed plus one crack on v2
  int i1 = ddr_find(a);  // unlimited cracks allowed plus one crack on v1
  int cnt = 0;
  for (int i=i1; i<i2; i++)
    if (arr[i]>=0) cnt++;
  return cnt;
}
