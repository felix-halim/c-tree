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

map<int,int> checksum {
  { 1, -295581051 },
  { 10, -48549690 },
  { 100, -2100668734 },
  { 1000, 1432162356 },
  { 10000, 1722829097 },
  { 100000, 1281102302 },
  { 1000000, -664634429 },
  { 10000000, 1113486809 },
  { 100000000, -769439752 },
};

void init(int *arr, int N);  // Initializes the initial values of N integers.
void insert(int value);      // Inserts the value.
int query(int value);        // Query for lower bound, returns 0 if not found.
void erase(int value);       // Deletes the value. The value guaranteed to exists.
void results(double insert_time, double query_time, int checksum);

int main(int argc, char *argv[]) {
  char *hostname = argv[1];
  int N = atoi(argv[2]);
  int Q = atoi(argv[3]);
  char *prog = argv[0];
  while (true) {
    char *p = strstr(prog, "/");
    if (p) prog = p + 1; else break;
  }
  printf("%lu,\"%s\",\"%s\",%d,%d,", system_clock::to_time_t(system_clock::now()), hostname, prog, N, Q);

  Random r(140384);
  int *iarr = new int[N];
  for (int i = 0; i < N; i++) iarr[i] = r.nextInt();

  int csum = 0;
  double insert_time = time_it([&] { init(iarr, N); });
  double query_time = time_it([&] {
    for (int i = 0; i < Q; i++) {
      csum = csum * 13 + query(r.nextInt());
      if (i % 1000 == 0) {
        for (int j = 0; j < 1000; j++) {
          int k = ((r.nextInt() % N) + N) % N;
          erase(iarr[k]);
          insert(iarr[k] = r.nextInt());
        }
      }
    }
  });

  results(insert_time, query_time, csum);
  assert(N != 100000000 || !checksum.count(Q) || checksum[Q] == csum);
  return 0;
}
