#include <cstdio>
#include <set>
#include <cassert>
#include <algorithm>
#include "test.h"

using namespace std;

#include "crack.h"

void init(unsigned *a, unsigned N_) {
  ci.clear();
  msize = 0;
  N = N_;
  marr = new unsigned[N*2];
  arr = new unsigned[N*2];     // for updates expansion
  for (int i = 0; i < N; i++)
    arr[i] = a[i];  // copy all
}

void insert(unsigned value) {
  if (pdel.count(value)){
    pdel.erase(pdel.lower_bound(value));    // don't insert if exists in pdel
  } else {
    pins.insert(value);
  }
}

void erase(unsigned value) {
  if (pins.count(value)){
    pins.erase(pins.lower_bound(value));    // don't delete if exists in pins
  } else {
    pdel.insert(value);
  }
}

unsigned view_query(unsigned a, unsigned b);

unsigned lower_bound(unsigned value) {
  return view_query(value, value + 1) ? value : 0;
}

unsigned select(unsigned a, unsigned b) {
  return view_query(a, b) ? a : 0;
}

void results(Statistics &s) {
}


unsigned view_query(unsigned a, unsigned b){
  merge_ripple(ci, arr, N, pins, pdel, a, b);  // merge qualified updates

  int L1,R1; find_piece(ci, N, a, L1,R1);
  int L2,R2; find_piece(ci, N, b, L2,R2);
  int i1 = R1, i2 = L2;
  
  msize = 0;

  if (L1==L2){    // a and b are in the same piece
    assert(R1==R2);
    if (L1 < R1){
      #ifdef MIN_PCSZ
        if (R1 - L1 > MIN_PCSZ){
          mdd1r_split_and_materialize<2>(L1,R1, a,b);
        } else {
          return crack(a,b);
        }
      #else
        mdd1r_split_and_materialize<2>(L1,R1, a,b);
      #endif
    }
  } else {    // b and b are in different piece, order doesn't matter
    if (L1 < R1){
      #ifdef MIN_PCSZ
        if (R1 - L1 > MIN_PCSZ){
          mdd1r_split_and_materialize<0>(L1,R1, a,b);  // do the same algo on the first piece
        } else {
          add_crack(ci, N, a, i1 = partition(arr, a, L1, R1));  // 2-split
        }
      #else
        mdd1r_split_and_materialize<0>(L1,R1, a,b);  // do the same algo on the first piece
      #endif
    }
    if (L2 < R2){
      #ifdef MIN_PCSZ
        if (R2 - L2 > MIN_PCSZ){
          mdd1r_split_and_materialize<1>(L2,R2, a,b);  // do the same algo on the second piece
        } else {
          add_crack(ci, N, b, i2 = partition(arr, b, L2, R2));  // 2-split
        }
      #else
        mdd1r_split_and_materialize<1>(L2,R2, a,b);  // do the same algo on the second piece
      #endif
    }
  }
  return msize + max(0, i2 - i1);
}

unsigned count_query(unsigned a, unsigned b){
  return view_query(a,b);
}
