#include <chrono>
#include <cstring>

#include <algorithm>

#include "checksum.h"
#include "workload.h"

using namespace std;
using namespace std::chrono;

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
  }
  return false;
}

void init(int *arr, int N);   // Initializes the initial values of N integers.
void insert(int value);       // Inserts the value.
int query(int value);         // Query for lower bound, returns 0 if not found.
void erase(int value);        // Deletes the value. The value guaranteed to exists.
void results(Statistics &s);  // Optionally fill in statistics.

static mt19937 gen(140384);
static uniform_int_distribution<> dis;
static int U, next_smallest = 1010000000;
static Statistics s;

static const char *update_workload[] = {
  "NOUP",   // 0. Read only queries.
  "LFHV",   // 1. Update 1000 tuples every 1000 queries.
  "HFLV",   // 2. Update 10 tuples every 10 queries.
  "QUEUE",  // 3. Remove the largest value, insert new smaller value than any existing value.
  "TRASH",  // 4. Insert values in the middle of the domain.
  "DELETE", // 5. Delete 1000 tuples every 1000 queries.
  "APPEND", // 6. Insert 10M tuples every 1000 queries. 
};

static vector<int> read_dataset(char *fn) {
  FILE *in = fopen(fn, "rb");
  if (!in) { fprintf(stderr,"Error opening file %s\n", fn); exit(1); }

  vector<int> arr;
  int tmp[1024];
  while (!feof(in)) {
    int N = fread(tmp, sizeof(int), 1024, in);
    for (int i = 0; i < N; i++) arr.push_back(tmp[i]);
  }

  if (ferror(in)) { fprintf(stderr,"Error reading %s!\n", fn); exit(1); }
  fclose(in);

  return arr;
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

  char *prog = argv[0];
  while (true) {
    char *p = strstr(prog, "/");
    if (p) prog = p + 1; else break;
  }
  s.algorithm = prog;

  vector<int> arr = read_dataset(argv[1]);

  int MAXQ;
  sscanf(argv[2], "%d", &MAXQ);
  sscanf(argv[3], "%lf", &s.selectivity);

  int W; sscanf(argv[4], "%d", &W);
  int mx = *max_element(arr.begin(), arr.end()) + 1;
  Workload query_workload(mx, W, s.selectivity);
  s.query_workload = workload_name[W];

  sscanf(argv[5], "%d", &U);
  s.update_workload = update_workload[U];

  s.N = arr.size();
  if (U != 6) s.N /= 2;

  if (U == 3) {                     // QUEUE update workload.
    for (int i = 0; i < s.N; i++)
      if (arr[i] <= next_smallest)
        arr[i] += next_smallest;
    sort(arr.begin(), arr.begin() + s.N);
  } else if (U == 6) {              // APPEND update workload.
    assert(s.N >= 500000000);
    s.N = 10000000;
  }

  // for (int i = 0; i < s.N * 2; i++) {
  //   if (i % 100 == 0)
  //   fprintf(stderr, "%d -> %d\n", i, arr[i]);
  // }

  fprintf(stderr, "N = %d\n", s.N);

  s.insert_time = time_it([&] { init(&arr[0], s.N); });
  s.checksum = 0;
  s.query_time = 0;

  fprintf(stderr, "I = %.6lf\n", s.insert_time);

  if (U == 5) {
    shuffle(arr.begin(), arr.end() + s.N, gen);
  }

  uniform_int_distribution<> disN(0, s.N - 1);
  int qidx = s.N - 1;
  for (s.Q = 1; ; s.Q *= 10) {
    double update_time = 0;
    s.query_time += time_it([&] {
      int nQ = s.Q - s.Q / 10; // nQ = how many queries needed.
      for (int i = 1, a, b; i <= nQ; i++) {
        if (U == 3) {
          a = next_smallest;
        } else {
          bool ok = query_workload.query(a,b); // get query endpoints based on the workload
          if (!ok){ s.Q = i; break; }
        }

        s.checksum = s.checksum * 13 + query(a);

        switch (U) {
          // NOUP.
          case 0: break;

          // LFHV.
          case 1: if (i % 1000 == 0) update_time += time_it([&] {
                    for (int j = 0; j < 1000; j++) {
                      int k = disN(gen);
                      erase(arr[k]);
                      int l = s.N + disN(gen);
                      swap(arr[k], arr[l]);
                      insert(arr[k]);
                    }
                  });
                  break;

          // HFLV.
          case 2: if (i % 10 == 0) update_time += time_it([&] {
                    for (int j = 0; j < 10; j++) {
                      int k = disN(gen);
                      erase(arr[k]);
                      int l = s.N + disN(gen);
                      swap(arr[k], arr[l]);
                      insert(arr[k]);
                    }
                  });
                  break;

          // QUEUE.
          case 3: if (i % 10 == 0) update_time += time_it([&] {
                    // if ((i + 1) % 1000000 == 0)
                    //   fprintf(stderr, "arr[%d] = %d, next = %d\n", qidx,arr[qidx],next_smallest);
                    for (int j = 0; j < 10; j++) {
                      erase(arr[qidx]);
                      insert(next_smallest);
                      arr[qidx--] = next_smallest--;
                      if (qidx < 0) qidx = s.N - 1;
                      assert(next_smallest >= 0);
                    }
                  });
                  break;

          // TRASH.
          case 4: if (i == 10000) update_time += time_it([&] {
                    // Insert 1M in the middle domain.
                    for (int j = 0; j < 1000000; j++) {
                      insert(arr[s.N + j]);
                    }
                  });
                  break;

          // DELETE.
          case 5: if (i % 1000 == 0) update_time += time_it([&] {
                    for (int j = 0; j < 1000; j++) {
                      erase(arr[--s.N]);
                    }
                  });
                  break;

          // APPEND.
          case 6: if (i % 1000 == 0) update_time += time_it([&] {
                    if (s.N < 580000000) {
                      // insert 10M tuples
                      for (int j = 0; j < 10000000; j++) {
                        insert(arr[s.N++]);
                      }
                    }
                  });
                  break;
        }
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
