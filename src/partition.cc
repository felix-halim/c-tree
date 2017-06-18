/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <functional>
#include <algorithm>
#include <vector>

#include "random.h"
#include "time_it.h"

using namespace std;

static Random rng(140384);

static void ctree_sort(long long *arr, long long *tmp, int n, int depth) {
  if (n < 30) {
    sort(arr, arr + n);
    return;
  }

  int nlo = 0;

  {
    long long P = arr[rng.nextInt(n)];
    int lo[n], hi[n];
    int nhi = n - 1;
    for (int i = 0; i < n; i++) {
      lo[nlo] = hi[nhi] = i;
      tmp[nhi] = arr[i];
      int j = arr[i] < P;
      nhi -= 1 - j;
      nlo += j;
    }
    for (int i = n - 1, j = nlo - 1; hi[i] < nlo; i--, j--) {
      arr[hi[i]] = arr[lo[j]];
    }
  }
 
  ctree_sort(arr, tmp, nlo, depth + 1);
  ctree_sort(tmp + nlo, arr + nlo, (n - nlo), depth + 1);
  memcpy(arr + nlo, tmp + nlo, sizeof(long long) * (n - nlo));
}

// Best so far.
static void ctree_sort2(long long *arr, long long *tmp, int n, int depth) {
  if (n < 30) {
    sort(arr, arr + n);
    return;
  }

  long long P = arr[rng.nextInt(n)];
  int nlo = 0, nhi = n - 1;
  for (int i = 0; i < n; i++) {
    tmp[nlo] = tmp[nhi] = arr[i];
    int j = arr[i] < P;
    nhi -= 1 - j;
    nlo += j;
  }
 
  ctree_sort2(tmp, arr, nlo, depth + 1);
  ctree_sort2(tmp + nlo, arr + nlo, (n - nlo), depth + 1);

  memcpy(arr, tmp, sizeof(long long) * nlo);
  memcpy(arr + nlo, tmp + nlo, sizeof(long long) * (n - nlo));
}

static void ctree_sort3(long long *arr, int n, int depth) {
  if (n < 30) {
    sort(arr, arr + n);
    return;
  }

  long long P = arr[rng.nextInt(n)];
  long long lo[n], hi[n];
  int nlo = 0, nhi = 0;
  for (int i = 0; i < n; i++) {
    if (arr[i] < P) {
      lo[nlo++] = arr[i];
    } else {
      hi[nhi++] = arr[i];
    }
    // lo[nlo] = hi[nhi] = arr[i];
    // int j = arr[i] < P;
    // nhi += 1 - j;
    // nlo += j;
  }
 
  ctree_sort3(lo, nlo, depth + 1);
  ctree_sort3(hi, nhi, depth + 1);

  memcpy(arr, lo, sizeof(long long) * nlo);
  memcpy(arr + nlo, hi, sizeof(long long) * nhi);
}

static void ctree_sort4(long long *arr, int n, int depth) {
  if (n < 30) {
    sort(arr, arr + n);
    return;
  }

  long long P = arr[rng.nextInt(n)];
  long long tmp[n];
  int nlo = 0, nhi = n - 1;
  for (int i = 0; i < n; i++) {
    tmp[nlo] = tmp[nhi] = arr[i];
    int j = arr[i] < P;
    nhi -= 1 - j;
    nlo += j;
  }
 
  ctree_sort4(tmp, nlo, depth + 1);
  ctree_sort4(tmp + nlo, (n - nlo), depth + 1);

  memcpy(arr, tmp, sizeof(long long) * nlo);
  memcpy(arr + nlo, tmp + nlo, sizeof(long long) * (n - nlo));
}


// #define BSIZE 8192
#define BSIZE 4096

static long long arr[BSIZE];

int main() {
  Random rng(140384);

  double total_time = 0;
  long long chk = 0;
  for (int m = 0; m < 30000; m++) {
    for (int i = 0; i < BSIZE; i++) {
      arr[i] = rng.nextLong();
    }

    total_time += time_it([&]() {
      // ctree_sort(arr, tmp, BSIZE, 0);
      long long tmp[BSIZE]; ctree_sort2(arr, tmp, BSIZE, 0);
      // ctree_sort3(arr, BSIZE, 0);
      // ctree_sort4(arr, BSIZE, 0);
    });

    chk = chk * 13 + arr[0];
    for (int i = 1; i < BSIZE; i++) {
      // assert(arr[i - 1] < arr[i]);
      chk = chk * 13 + arr[i];
    }
  }

  fprintf(stderr, "tt = %.6lf, chk = %lld\n", total_time, chk);
  assert(chk == -3719803178456933530LL);
}
