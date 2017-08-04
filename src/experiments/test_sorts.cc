/*
To run, go to the root folder and execute:

make test_sorts
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>

#include "../random.h"
#include "../time_it.h"
#include "ctree.h"
// #include "parallel_sort.h"
#include "ska_sort.h"
#include "vergesort.h"

#ifdef DBG
#define MAXN (1 << 20)
#define BSIZE (1 << 10)
#else
#define MAXN (1 << 26)
#define BSIZE (1 << 12)
#endif

static long long arr[MAXN];

using namespace std;

#define KRED "\x1B[31m"
#define KGRN "\x1B[32m"
#define KYEL "\x1B[33m"
#define KBLU "\x1B[34m"
#define KMAG "\x1B[35m"
#define KCYN "\x1B[36m"
#define KWHT "\x1B[37m"
#define RESET "\x1B[0m"

// Run tests on sort_algo for different level of sortedness.
static void run(const char *name,
                std::function<void(long long[], int)> sort_algo) {
  fprintf(stderr, "%s (N = %d):\n", name, MAXN);
  Random random(140384);
  int p100 = 100000;
  for (int p = p100;; p /= 2) {
    int misplaced = 0;
    long long sorted_chk = 0;
    for (int i = 0; i < MAXN; i++) {
      arr[i] = i;
      sorted_chk = sorted_chk * 13 + i;
      if (random.nextInt(p100) < p) {
        swap(arr[i], arr[random.nextInt(i + 1)]);
        misplaced++;
      }
    }

    double sort_time = time_it([&]() { sort_algo(arr, MAXN); });

    long long chk = 0;
    for (int i = 0; i < MAXN; i++) {
      chk = chk * 13 + arr[i];
    }
    fprintf(stderr, "%7.3lf%% sorted: %10.6lf s, rswaps = %8d%s\n", p / 1000.0,
            sort_time, misplaced,
            chk == sorted_chk ? "" : (" " KRED "[CHK Failed]" RESET));
    if (p == 0) {
      break;
    }
  }
  fprintf(stderr, "\n");
}

static void std_sort(long long a[], int n) { std::sort(a, a + n); }
static void verge_sort(long long a[], int n) { vergesort::vergesort(a, a + n); }
static void skasort(long long a[], int n) { ska_sort(a, a + n); }
static void ctree_sort(long long a[], int n) { ctreesort<BSIZE>(a, n); }

int main() {
  run("ctree_sort", ctree_sort);
  run("std::sort", std_sort);
  // run("std::parallel_sort", parallel_sort);
  run("vergesort", verge_sort);
  run("ska_sort", skasort);
}
