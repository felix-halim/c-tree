#include <cstdio>
#include <set>
#include <cassert>
#include <algorithm>
#include "test.h"

using namespace std;


#include "crack.h"


void init(int *a, int N_) {
  ci.clear();
  msize = 0;
  N = N_;
  marr = new int[N*2];
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

int lower_bound(int value) {
  return view_query(value, value + 1) ? value : 0;
}

int select(int a, int b) {
  return view_query(a, b) ? a : 0;
}

void results(Statistics &s) {
}
