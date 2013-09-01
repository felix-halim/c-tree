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
  { 10, 727569345 },
  { 100, -905374377 },
  { 1000, -354145257 },
  { 10000, 322147057 },
  { 100000, -1715352722 },
  { 1000000, 174351911 },
  { 10000000, -141256775 },
  { 100000000, -1550769227 },
};

map<int,int> checksum8 {
  { 1, -295581051 },
  { 10, 441507466 },
  { 100, 1214133662 },
  { 1000, 1588124755 },
  { 10000, 980068087 },
  { 100000, 1735881447 },
  { 1000000, -831586212 },
  { 10000000, -813580773 },
  { 100000000, 512824327 },
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
  for (int Q = 1, cur = 0; Q <= 1000000000; Q *= 10) {
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
  }
  return 0;
}
