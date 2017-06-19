/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <vector>

#include "random.h"
#include "ska_sort.h"
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

inline unsigned long long to_unsigned(long long l) {
  return static_cast<unsigned long long>(l) +
         static_cast<unsigned long long>(1ll << (sizeof(long long) * 8 - 1));
}

inline long long to_signed(unsigned long long l) {
  return l -
         static_cast<unsigned long long>(1ll << (sizeof(long long) * 8 - 1));
}

static void radix_sort(long long *arr, int n, int depth) {
  if (n < 256) {
    sort(arr, arr + n);
    return;
  }

  int occ[256]{0};

  int nshift = (7 - depth) * 8;
  int n16 = (n / 16) * 16;
  for (int i = 0; i < n16; i+=16) {
    for (int j = 0; j < 16; j++) {
      occ[(to_unsigned(arr[i+j]) >> nshift) & 255]++;
    }
  }
  for (int i = n16; i < n; i++) {
    occ[(to_unsigned(arr[i]) >> nshift) & 255]++;
  }

  int idx[256]{0};
  for (int i = 1; i < 256; i++) {
    idx[i] = idx[i - 1] + occ[i - 1];
  }
  for (int i = 1; i < 256; i++) {
    occ[i] += occ[i - 1];
  }

  for (int p = 255; p >= 0; p--) {
    while (idx[p] < occ[p]) {
      int rem = occ[p] - idx[p];
      // fprintf(stderr, "p = %d, rem = %d\n", p, rem);
      for (int i = 0, j = idx[p]; i < rem; i++, j++) {
        int k = (to_unsigned(arr[j]) >> nshift) & 255;
        swap(arr[idx[k]++], arr[j]);
      }
    }
  }

  for (int i = 0; i < 256; i++) {
    // fprintf(stderr, "%d => %d %d\n", i, idx[i], occ[i]);
    assert(idx[i] == occ[i]);
    int prev = i == 0 ? 0 : occ[i - 1];
    int sz = occ[i] - prev;
    if (sz > 0) { 
      radix_sort(arr + prev, sz, depth + 1);
    }
  }
  // fprintf(stderr, "%d %d %d\n", cnt, n, depth);
}

// #define BSIZE 8192
#define BSIZE 100000000

static long long arr[BSIZE];
static long long tmp[BSIZE];

int main() {
  Random rng(140384);

  double total_time = 0;
  long long chk = 0;
  for (int m = 0; m < 1; m++) {
    for (int i = 0; i < BSIZE; i++) {
      arr[i] = rng.nextLong();
    }

    total_time += time_it([&]() {
      // ctree_sort(arr, tmp, BSIZE, 0);
      // sort(arr, arr + BSIZE);

      // sort(arr, arr + BSIZE, [](long long a, long long b) {
      //   return to_unsigned(a) < to_unsigned(b);
      // });

      // ctree_sort2(arr, tmp, BSIZE, 0);

      // ska_sort(arr, arr + BSIZE, [&](const long long &e) { return e; });

      radix_sort(arr, BSIZE, 0);

      // ctree_sort3(arr, BSIZE, 0);
      // ctree_sort4(arr, BSIZE, 0);
    });
    fprintf(stderr, "done\n");

    chk = chk * 13 + arr[0];
    for (int i = 1; i < BSIZE; i++) {
      // assert(arr[i - 1] < arr[i]);
      chk = chk * 13 + arr[i];
    }
    fprintf(stderr, "done2\n");
  }

  fprintf(stderr, "tt = %.6lf, chk = %lld\n", total_time, chk);
  assert(chk == 894150326305700963LL);
  // assert(chk == -3719803178456933530LL);
}
