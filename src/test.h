#include <chrono>
#include <cstring>

#include <algorithm>

#include "query.h"
#include "update.h"
#include "checksum.h"

using namespace std;
using namespace std::chrono;

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

struct Statistics {
  // Managed by test framework.
  int N;
  int Q;
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
  int n_leaves;
  int n_capacity;
  int n_internals;
  int max_depth;
  int slack;
  int in_size;
  int ln_size;
  int ia_free;
  int ia_size;
  int la_free;
  int la_size;
  string note;

  void print_header() {
    printf("timestamp,algorithm,query_workload,update_workload,N,Q,selectivity,verified,insert_time,update_time,query_time,checksum,"
      "n_leaves,n_capacity,n_internals,max_depth,slack,in_size,ln_size,ia_free,ia_size,la_free,la_size\n");
  }

  void print() {
    printf("%lu,", system_clock::to_time_t(system_clock::now()));
    printf("\"%s\",", algorithm.c_str());
    printf("\"%s\",", query_workload.c_str());
    printf("\"%s\",", update_workload.c_str());
    printf("%d,", N);
    printf("%d,", Q);
    printf("%lf,", selectivity);
    printf("%d,", verified);
    printf("%.6lf,", insert_time);
    printf("%.6lf,", update_time);
    printf("%.6lf,", query_time);
    printf("%llu,", checksum);
    printf("%d,", n_leaves);
    printf("%d,", n_capacity);
    printf("%d,", n_internals);
    printf("%d,", max_depth);
    printf("%d,", slack);
    printf("%d,", in_size);
    printf("%d,", ln_size);
    printf("%d,", ia_free);
    printf("%d,", ia_size);
    printf("%d,", la_free);
    printf("%d,", la_size);
    printf("\"%s\"", note.c_str());
    puts("");
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

void init(int *arr, int N);   // Initializes the initial values of N integers.
void insert(int value);       // Inserts the value.
int lower_bound(int value);   // Query for lower bound, returns 0 if not found.
int select(int a, int b);     // Select values from [a, b), without fetching the values.
void erase(int value);        // Deletes the value. The value guaranteed to exists.
void results(Statistics &s);  // Optionally fill in statistics.

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

  int MAXQ, W, U;
  s.algorithm = algorithm_name(argv[0]);
  sscanf(argv[2], "%d", &MAXQ);
  sscanf(argv[3], "%lf", &s.selectivity);
  sscanf(argv[4], "%d", &W);
  sscanf(argv[5], "%d", &U);
  Update update(argv[1], U);

  Workload query_w(W, s.selectivity);
  s.query_workload = workload_name[W];
  s.update_workload = update_workload[U];

  if (U != 6) {
    update.load();
    s.N = update.size() / 2;
    query_w.set_max(update.max_element() + 1);
    if (U == 3) update.prepare_queue(s.N);
  } else {
    update.load(100000);
    s.N = update.size();
    query_w.set_max(update.max_element() + 1);
  }

  fprintf(stderr, "N = %d\n", s.N);

  s.insert_time = time_it([&] { init(update.get_arr(), s.N); });

  fprintf(stderr, "I = %.6lf\n", s.insert_time);

  if (U == 5) {
    update.prepare_deletion(s.N);
    MAXQ = min(s.N, MAXQ);
  }

  for (s.Q = 1; ; s.Q *= 10) {
    double update_time = 0;
    double load_time = 0;
    s.query_time += time_it([&] {
      int nQ = s.Q - s.Q / 10; // nQ = how many queries needed.
      for (int i = 1, a, b; i <= nQ; i++) {
        // if (U == 3) {
        //   a = update.get_next_smallest();
        // } else {
          bool ok = query_w.query(a,b); // get query endpoints based on the workload
          if (!ok){ s.Q = i; MAXQ = -1; break; }
        // }

        s.checksum = s.checksum * 13 + lower_bound(a);

        switch (U) {
          // NOUP.
          case 0: break;

          // LFHV.
          case 1: if (i % 1000 == 0) update_time += time_it([&] {
                    REP(j, 1000) {
                      update.update(s.N, a, b);
                      erase(a);
                      insert(b);
                    }
                  });
                  break;

          // HFHV.
          case 2: if (i % 10 == 0) update_time += time_it([&] {
                    REP(j, 1000) {
                      update.update(s.N, a, b);
                      erase(a);
                      insert(b);
                    }
                  });
                  break;

          // QUEUE.
          case 3: if (i % 10 == 0) update_time += time_it([&] {
                    REP(j, 10) {
                      update.update_queue(s.N, a, b);
                      erase(a);
                      insert(b);
                    }
                  });
                  break;

          // TRASH.
          case 4: if (i == 10000) update_time += time_it([&] {
                    int *arr = update.get_arr();
                    REP(j, 1000000) {
                      insert(arr[s.N + j]);
                      // fprintf(stderr, "%d \n", arr[s.N + j]);
                    }
                  });
                  break;

          // DELETE.
          case 5: if (i % 1000 == 0) update_time += time_it([&] {
                    REP(j, 1000) {
                      erase(update.update_delete(s.N));
                    }
                  });
                  break;

          // APPEND.
          case 6: if (i % 1 == 0) update_time += time_it([&] {
                    if (MAXQ != -1) {
                      update.clear();
                      bool loaded = false;
                      load_time += time_it([&] {
                        loaded = update.load(100000);
                        query_w.set_max(update.max_element());
                        s.N += update.size();
                      });
                      if (loaded) {
                        int *arr = update.get_arr();
                        REP(j, update.size()) insert(arr[j]);
                      } else {
                        MAXQ = -1;
                      }
                    }
                  });
                  break;
        }
      }
    });
    update_time -= load_time;
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
