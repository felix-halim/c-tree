#include <cstring>
#include <chrono>
#include <map>
#include <algorithm>
#include "random.h"

using namespace std;
using namespace std::chrono;

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

map<int,int> checksum7 {
  { 1, 1179274802 },
  { 10, 365664126 },
  { 100, 708307941 },
  { 1000, 276117111 },
  { 10000, 1032820950 },
  { 100000, -1810972451 },
  { 1000000, -1937611658 },
  { 10000000, 1647471467 },
  { 100000000, -527620924 },
  { 1000000000, 2144904532 },
};

map<int,int> checksum8 {
  { 1, -295581051 },
  { 10, 1566733998 },
  { 100, -1248672030 },
  { 1000, -977199412 },
  { 10000, 778706157 },
  { 100000, -1587835339 },
  { 1000000, 425986485 },
  { 10000000, -1248377164 },
  { 100000000, -1370805606 },
  { 1000000000, -1643017426 },
};

void init(int *arr, int N);  // Initializes the initial values of N integers.
void insert(int value);      // Inserts the value.
int query(int value);        // Query for lower bound, returns 0 if not found.
void erase(int value);       // Deletes the value. The value guaranteed to exists.
void results(double insert_time, double query_time, int checksum);

int main(int argc, char *argv[]) {
  char *prog = argv[0];
  while (true) {
    char *p = strstr(prog, "/");
    if (p) prog = p + 1; else break;
  }
  char *hostname = argv[1];
  int N = atoi(argv[2]);

  Random rng(140384);
  int *iarr = new int[N];
  for (int i = 0; i < N; i++) iarr[i] = rng.nextInt();

  double insert_time = time_it([&] { init(iarr, N); });

  int csum = 0;
  double query_time = 0;
  for (int Q = 1, cur = 0; ; Q *= 10) {
    int nQ = Q - cur; cur = Q;
    query_time += time_it([&] {
      for (int i = 0; i < nQ; i++) {
        csum = csum * 13 + query(rng.nextInt());
        if (i % 1000 == 0) {
          for (int j = 0; j < 1000; j++) {
            int k = ((rng.nextInt() % N) + N) % N;
            erase(iarr[k]);
            insert(iarr[k] = rng.nextInt());
          }
        }
      }
    });
    printf("%lu,\"%s\",\"%s\",%d,%d,", system_clock::to_time_t(system_clock::now()), hostname, prog, N, Q);
    results(insert_time, query_time, csum);
    if (N == 10000000 && checksum7.count(Q) && checksum7[Q] != csum)
      fprintf(stderr, "\033[1;31mFAILED\033[0m checksum %d != %d\n", checksum7[Q], csum);
    if (N == 100000000 && checksum8.count(Q) && checksum8[Q] != csum)
      fprintf(stderr, "\033[1;31mFAILED\033[0m checksum %d != %d\n", checksum8[Q], csum);
    if (Q == 1000000000) break;
  }
  return 0;
}
