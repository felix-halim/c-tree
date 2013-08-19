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
  { 100000, -1927822250 },
  { 1000000, 880684722 },
  { 10000000, -645327377 },
  { 100000000, 1814544708 },
  { 1000000000, 1636094096 },
};

void init(int *arr, int N);  // Initializes the initial values of N integers.
int query(int value);        // Query for lower bound, returns 0 if not found.
void results(double insert_time, double query_time, int checksum);

int main(int argc, char *argv[]) {
  int N = atoi(argv[1]);
  int Q = atoi(argv[2]);

  Random r(140384);
  int *iarr = new int[N];
  int *queries = new int[Q];
  for (int i = 0; i < N; i++) iarr[i] = r.nextInt();
  for (int i = 0; i < Q; i++) queries[i] = r.nextInt();

  int csum = 0;
  double insert_time = time_it([&] { init(iarr, N); });
  double query_time = time_it([&]{
    for (int i = 0; i < Q; i++)
      csum = csum * 13 + query(queries[i]);
  });

  results(insert_time, query_time, csum);
  assert(N != 100000000 || !checksum.count(Q) || checksum[Q] == csum);
  return 0;
}
