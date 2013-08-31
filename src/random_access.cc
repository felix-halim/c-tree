#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <chrono>
#include "random.h"

using namespace std;
using namespace chrono;

#define SIZE_MB 1500
#define N (SIZE_MB * 1024 * 1024 / 4)

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

int *arr;
Random rng(140384);

void random_access(int MAXN, int Q) {
  fprintf(stderr, "N = %d, Q = %d\n", N, Q);
  arr = new int[N];
  double t1 = time_it([&] {
    for (int i = 0; i < Q; i++) {
      int j = rng.nextInt(((min(N, MAXN) >> 0) << 0) - 64);
      for (int k = 0; k < 64; k++, j++) {
        arr[j] = rng.nextInt();
      }
    }
  });

  int csum = 0;
  for (int i = 0; i < N; i++) {
    csum += arr[i];
  }

  fprintf(stderr, "time = %.6lf, csum = %d\n", t1, csum);
}

void linear_alloc(int MAXN) {
  vector<int*> v;
  double t1 = time_it([&] {
    for (int i = 0; i < MAXN; i++) {
      int *a = new int[64];
      for (int k = 0; k < 64; k++) {
        a[k] = rng.nextInt();
      }
      v.push_back(a);
    }
  });

  int csum = 0;
  for (int i = 0; i < (int) v.size(); i++) {
    csum += v[i][rng.nextInt(64)];
  }

  fprintf(stderr, "time = %.6lf, csum = %d\n", t1, csum);
}

int main(int argc, char *argv[]) {
  // random_access(atoi(argv[1]), atoi(argv[2]));
  double t1 = time_it([&]{
    arr = new int[100000000];
    for (int i = 0; i < 100000000; i++) {
      arr[i] = rng.nextInt();
    }
  });
  int *arr2;
  double t2 = time_it([&]{
    arr2 = new int[100000000];
    memcpy(arr2, arr, sizeof(int) * 100000000);
  });
  fprintf(stderr, "t1 = %.6lf, t2 = %.6lf, %d\n", t1, t2, arr2[rng.nextInt(100000000)]);
  linear_alloc(atoi(argv[1]));
}
