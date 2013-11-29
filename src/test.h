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

void init(unsigned *arr, unsigned N);   // Initializes the initial values of N integers.
unsigned lower_bound(unsigned value);   // Query for lower bound, returns 0 if not found.
unsigned select(unsigned a, unsigned b);     // Select values from [a, b), without fetching the values.
unsigned count(unsigned a, unsigned b);      // Count values in range [a, b).
unsigned sum(unsigned a, unsigned b);        // Sum values in range [a, b).
void insert(unsigned value);       // Inserts the value.
void erase(unsigned value);        // Deletes the value. The value guaranteed to exists.
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
  Update update(argv[1], U, W, [&](unsigned mx) { query_w.set_max(mx + 1); });
  fprintf(stderr, "done. ");
  Statistics s(parse_algorithm_name(argv[0]), query_w.name(), update.name());

  fprintf(stderr, "N = %u, ", update.get_n());

  sort(update.get_arr(), update.get_arr() + update.get_n());
  double insert_time = time_it([&] { init(update.get_arr(), update.get_n()); });

  fprintf(stderr, "I = %.3lf, Q = %lld\n", insert_time, query_w.maxq());

  vector<long long> samples = generate_samples(query_w.maxq());
  double total_query_time = 0, total_update_time = 0;
  unsigned long long checksum = 0;
  for (long long i = 1, *next_sample = &samples[0], Q = 0; ; ) {
    double load_time = 0; // Should not be counted. It is time to load data from disk.
    double update_time = 0;
    double runtime = time_it([&] {
      long long nQ = - *next_sample; nQ += *(++next_sample); Q += nQ;

      for (unsigned a, b; nQ--; i++) {
        bool ok = query_w.query(a,b); // get query endpoints based on the workload
        assert(ok || W == 0);

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

        update.execute(i,
          [](unsigned v){ insert(v); }, // insert.
          [](unsigned v){ erase(v); }, // erase.
          [&](double load_t) { load_time += load_t; },
          [&](double update_t) { update_time += update_t; }
        );
      }
    });
    runtime -= load_time;
    total_query_time += runtime - update_time;
    total_update_time += update_time;

    results(s);

    s.print(Q, query_w.selectivity(), verify(update.code(), update.get_n(), Q, checksum),
      insert_time, total_query_time, total_update_time, checksum);

    if (Q >= query_w.maxq() || runtime > 3600 * 5) break;
  }
}
