/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <functional>
#include <algorithm>
#include <vector>

#include "ska_sort.h"
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

inline unsigned long long to_unsigned(long long l)
{
    return static_cast<unsigned long long>(l) + static_cast<unsigned long long>(1ll << (sizeof(long long) * 8 - 1));
}

inline long long to_signed(unsigned long long l)
{
    return l - static_cast<unsigned long long>(1ll << (sizeof(long long) * 8 - 1));
}

static void radix_sort(long long *arr, long long *tmp, int n, int depth) {
  if (n < 256) {
    sort(arr, arr + n);
    return;
  }
  int occ[256] {0};

  int nshift = (7 - depth) * 8;
  int n4 = (n/4) * 4;
  for (int i = 0; i < n4; i+=4) {
    occ[(to_unsigned(arr[i]) >> nshift) & 255]++;
    occ[(to_unsigned(arr[i+1]) >> nshift) & 255]++;
    occ[(to_unsigned(arr[i+2]) >> nshift) & 255]++;
    occ[(to_unsigned(arr[i+3]) >> nshift) & 255]++;
  }
  for (int i = n4; i < n; i++) {
    occ[(to_unsigned(arr[i]) >> nshift) & 255]++;
  }

  int idx[256] {0};
  for (int i = 1; i < 256; i++) {
    idx[i] = idx[i - 1] + occ[i - 1];
  }

  for (int i = 0; i < n; i++) {
    tmp[idx[(to_unsigned(arr[i]) >> nshift) & 255]++] = arr[i];
  }

  for (int i = 0, lo = 0; i < 256; i++) {
    if (occ[i] > 1) {
      radix_sort(tmp + lo, arr + lo, occ[i], depth + 1);
    }
    lo += occ[i];
  }

  memcpy(arr, tmp, sizeof(long long) * n);
}

// #define BSIZE 8192
#define BSIZE 100000000

static long long arr[BSIZE];
static long long tmp[BSIZE];

struct fhitem {
  unsigned char arr[8];
  fhitem(){}
  fhitem(long long num) {
    unsigned long long t = to_unsigned(num);
    for (int i = 0; i < 8; i++) {
      arr[7 - i] = t & 255;
      t >>= 8;
    }
  }
  long long num() {
    unsigned long long t = arr[0];
    for (int i = 1; i < 8; i++) {
      t = (t << 8) | arr[i];
    }
    return to_signed(t);
  }
  bool operator()(const fhitem const &that) const {
    for (int i = 0; i < 8; i++) {
      if (arr[i] != that.arr[i]) {
        return arr[i] < that.arr[i];
      }
    }
    return false;
  }
};

bool fhitemcmp(const fhitem &a, const fhitem &b) {
  for (int i = 0; i < 8; i++) {
    if (a.arr[i] != b.arr[i]) {
      return a.arr[i] < b.arr[i];
    }
  }
  return false;
}

static void radix_fhitems_sort(fhitem *arr, fhitem *tmp, int n, int depth) {
  if (n < 256) {
    sort(arr, arr + n, fhitemcmp);
    return;
  }
  int occ[256] {0};

  int nshift = (7 - depth) * 8;
  for (int i = 0; i < n; i++) {
    occ[arr[i].arr[depth]]++;
  }

  int idx[256] {0};
  for (int i = 1; i < 256; i++) {
    idx[i] = idx[i - 1] + occ[i - 1];
  }

  for (int i = 0; i < n; i++) {
    tmp[idx[arr[i].arr[depth]]] = arr[i];
  }

  for (int i = 0, lo = 0; i < 256; i++) {
    if (occ[i] > 1) {
      radix_fhitems_sort(tmp + lo, arr + lo, occ[i], depth + 1);
    }
    lo += occ[i];
  }

  memcpy(arr, tmp, sizeof(fhitem) * n);
}


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

      // fhitem *fhitems = new fhitem[BSIZE];
      // fhitem *fhitems_tmp = new fhitem[BSIZE];
      // for (int i = 0; i < BSIZE; i++) {
      //   fhitems[i] = fhitem(arr[i]);
      // }
      // radix_fhitems_sort(fhitems, fhitems_tmp, BSIZE, 0);
      // for (int i = 0; i < BSIZE; i++) {
      //   arr[i] = fhitems[i].num();
      // }

      unsigned long long * uarr = new unsigned long long[BSIZE];
      for (int i = 0; i < BSIZE; i++) {
        uarr[i] = to_unsigned(arr[i]);
      }
      for (int i = 0; i < BSIZE; i++) {
        arr[i] = to_signed(uarr[i]);
      }
      // radix_sort(arr, tmp, BSIZE, 0);

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
