#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <algorithm>

using namespace std;
using namespace std::chrono;

#include "util.h"
#include "query.h"
#include "update.h"

void init(int *arr, int N);   // Initializes the initial values of N integers.
int lower_bound(int value);   // Query for lower bound, returns 0 if not found.
int select(int a, int b);     // Select values from [a, b), without fetching the values.
int count(int a, int b);      // Count values in range [a, b).
int sum(int a, int b);        // Sum values in range [a, b).
void insert(int value);       // Inserts the value.
void erase(int value);        // Deletes the value. The value guaranteed to exists.
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
  sscanf(argv[2], "%lld", &MAXQ);
  double sel = atof(argv[3]);
  int W = atoi(argv[4]);
  int U = atoi(argv[5]);
  Workload query_w(MAXQ, sel, W);
  fprintf(stderr, "Loading ... ");
  Update update(argv[1], U, W, [&](int mx) { query_w.set_max(mx + 1); });
  fprintf(stderr, "done. ");
  Statistics s(parse_algorithm_name(argv[0]), query_w.name(), update.name());

  fprintf(stderr, "N = %d, ", update.get_n());

  double insert_time = time_it([&] { init(update.get_arr(), update.get_n()); });

  fprintf(stderr, "I = %.3lf, Q = %lld\n", insert_time, query_w.maxq());

  vector<long long> samples = generate_samples(query_w.maxq());
  double total_query_time = 0, total_update_time = 0;
  unsigned long long checksum = 0;
  for (long long i = 1, *next_sample = &samples[0], Q = 0; ; ) {
    double update_time = 0;
    double runtime = time_it([&] {
      long long nQ = - *next_sample; nQ += *(++next_sample); Q += nQ;

      for (int a, b; nQ--; i++) {
        bool ok = query_w.query(a,b); // get query endpoints based on the workload
        // assert(ok);

        checksum = checksum * 13 + 
          #if defined(COUNT_QUERY)
            count(a, b);
          #elif defined(SUM_QUERY)
            sum(a, b);
          #elif defined(SELECT_QUERY)
            select(a, b);
          #else
            lower_bound(a);
          #endif

        update_time += update.execute(i,
          [](int v){ insert(v); }, // insert.
          [](int v){ erase(v); } // erase.
        );
      }
    });
    total_query_time += runtime - update_time;
    total_update_time += update_time;

    results(s);

    s.print(update.get_n(), Q, query_w.selectivity(), verify(update.code(), update.get_n(), Q, checksum),
      insert_time, total_query_time, total_update_time, checksum);

    if (Q >= query_w.maxq()) break;
  }
}
