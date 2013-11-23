#include <chrono>
#include <cstring>

#include <algorithm>

using namespace std;
using namespace std::chrono;

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)

template<typename F>
static double time_it(F f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

void insert(int value);       // Inserts the value.
void erase(int value);        // Deletes the value. The value guaranteed to exists.

#include "query.h"
#include "update.h"
#include "checksum.h"

struct Statistics {
  // Managed by test framework.
  int N;
  long long Q;
  double selectivity;
  int verified;
  string algorithm;
  string query_workload;
  string update_workload;
  double insert_time;
  double query_time;
  double update_time;
  unsigned long long checksum;

  // Each algorithm optionally fill in the following stats at the end of run:
  int n_index;
  int n_bytes;
  int n_slack;
  int n_internal;
  int n_leaf;
  int bt_int_sz;
  int bt_leaf_sz;
  int large_touch;
  int small_touch;
  int art_n4;
  int art_n16;
  int art_n48;
  int art_n256;

  void print_header() {
    printf("timestamp,algorithm,query_workload,update_workload,N,Q,selectivity,verified,insert_time,update_time,query_time,checksum,"
      "n_index,n_bytes,n_slack,n_internal,n_leaf,bt_int_sz,bt_leaf_sz,large_touch,small_touch,art_n4,art_n16,art_n48,art_n256,query_type\n");
  }

  void print() {
    printf("%lu", system_clock::to_time_t(system_clock::now()));
    printf(",\"%s\"", algorithm.c_str());
    printf(",\"%s\"", query_workload.c_str());
    printf(",\"%s\"", update_workload.c_str());
    printf(",%d", N);
    printf(",%lld", Q);
    printf(",%lf", selectivity);
    printf(",%d", verified);
    printf(",%.6lf", insert_time);
    printf(",%.6lf", update_time);
    printf(",%.6lf", query_time);
    printf(",%llu", checksum);
    printf(",%d", n_index);
    printf(",%d", n_bytes);
    printf(",%d", n_slack);
    printf(",%d", n_internal);
    printf(",%d", n_leaf);
    printf(",%d", bt_int_sz);
    printf(",%d", bt_leaf_sz);
    printf(",%d", large_touch);
    printf(",%d", small_touch);
    printf(",%d", art_n4);
    printf(",%d", art_n16);
    printf(",%d", art_n48);
    printf(",%d", art_n256);

    #if defined(COUNT_QUERY)
      puts(",COUNT");
    #elif defined(SUM_QUERY)
      puts(",SUM");
    #elif defined(SELECT_QUERY)
      puts(",SELECT");
    #else
      puts(",LOWER_BOUND");
    #endif
  }
};

static int verify(int U, Statistics &s) {
  switch (U) {
    case 0 : return (s.N == 10000000) ? checksum_match(noup_checksum7, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match(noup_checksum8, s.Q, s.checksum) : 0);
    case 1 : return (s.N == 10000000) ? checksum_match(lfhv_checksum7, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match(lfhv_checksum8, s.Q, s.checksum) : 0);
    case 2 : return (s.N == 10000000) ? checksum_match2(skew_checksum7, s.selectivity, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match2(skew_checksum8, s.selectivity, s.Q, s.checksum) : 0);
    case 3 : return s.N == 100000000 ? checksum_match(queue_checksum8, s.Q, s.checksum) : 0;
    case 4 : return s.N == 100000000 ? checksum_match(trash_checksum8, s.Q, s.checksum) : 0;
    case 5 : return s.N == 100000000 ? checksum_match(delete_checksum8, s.Q, s.checksum) : 0;
    case 6 : return s.N == 100000000 ? checksum_match(append_checksum8, s.Q, s.checksum) : 0;
  }
  return false;
}

static mt19937 gen(140384);
static uniform_int_distribution<> dis;
static Statistics s;

static char* algorithm_name(char *prog) {
  while (true) {
    char *p = strstr(prog, "/");
    if (p) prog = p + 1; else break;
  }
  return prog;
}

static vector<long long> generate_samples(long long MAXQ) {
  vector<long long> samples;
  samples.push_back(0);
  samples.push_back(MAXQ);
  for (long long i = 1; i <= MAXQ; i *= 2) samples.push_back(i);
  for (long long i = 1; i <= MAXQ; i *= 10) {
    samples.push_back(i);
    for (long long j = 0; j < 100; j += 5) samples.push_back(i * j / 100);
  }
  sort(samples.begin(), samples.end());
  samples.erase(unique(samples.begin(), samples.end()), samples.end());
  return samples;
}

void init(int *arr, int N);   // Initializes the initial values of N integers.
int lower_bound(int value);   // Query for lower bound, returns 0 if not found.
int select(int a, int b);     // Select values from [a, b), without fetching the values.
int count(int a, int b);      // Count values in range [a, b).
int sum(int a, int b);        // Sum values in range [a, b).
void results(Statistics &s);  // Optionally fill in statistics.

int main(int argc, char *argv[]) {
  if (argc != 6) {
    fprintf(stderr,
      "usage: ./%s "        // 0
      "input_file "         // 1
      "num_of_queries "     // 2
      "selectivity "        // 3
      "query_workload "     // 4
      "update_workload\n",  // 5
      argv[0]);
    exit(1);
  }

  long long MAXQ;
  int W, U;
  s.algorithm = algorithm_name(argv[0]);
  sscanf(argv[2], "%lld", &MAXQ);
  sscanf(argv[3], "%lf", &s.selectivity); // Only used if COUNT_QUERY or SUM_QUERY or SELECT_QUERY is defined.
  sscanf(argv[4], "%d", &W);
  sscanf(argv[5], "%d", &U);

  Update update(argv[1], U);

  Workload query_w(W, s.selectivity, MAXQ);
  query_w.set_max(update.max_element() + 1);

  s.N = update.the_N();

  fprintf(stderr, "N = %d\n", s.N);

  s.query_workload = workload_name[W];
  s.update_workload = update_workload[U];

  s.insert_time = time_it([&] { init(update.get_arr(), s.N); });

  fprintf(stderr, "I = %.3lf, Q = %lld\n", s.insert_time, MAXQ);

  vector<long long> samples = generate_samples(MAXQ);
  long long i = 1, *next_sample = &samples[0];
  for (s.Q = 0; ; ) {
    double update_time = 0;
    s.query_time += time_it([&] {
      long long nQ = - *next_sample; nQ += *(++next_sample); s.Q += nQ;

      for (int a, b; nQ--; ) {
        bool ok = query_w.query(a,b); // get query endpoints based on the workload
        assert(ok);

        s.checksum = s.checksum * 13 + 
          #if defined(COUNT_QUERY)
            count(a, b);
          #elif defined(SUM_QUERY)
            sum(a, b);
          #elif defined(SELECT_QUERY)
            select(a, b);
          #else
            lower_bound(a);
          #endif

        update_time += update.execute(i++);
      }
    });
    s.query_time -= update_time;
    s.update_time += update_time;
    results(s);
    s.verified = verify(U, s);
    s.print();
    fflush(stdout);
    if (s.Q >= MAXQ) break;
  }

  return 0;
}
