#include <chrono>
#include <map>
#include <cstring>
#include <algorithm>

using namespace std;
using namespace std::chrono;

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

map<int, unsigned long long> noup_checksum {
  { 1, 815231912ULL },
  { 10, 8837347545989574115ULL },
  { 100, 13550455485803658225ULL },
  { 1000, 8650331494777728818ULL },
  { 10000, 10629385045258940980ULL },
  { 100000, 7278354990761420101ULL },
  { 1000000, 3156700812791593580ULL },
  { 10000000, 754841950488114731ULL },
  { 100000000, 9956823453743028538ULL },
  { 1000000000, 2574210756483803256ULL },
};

map<int, unsigned long long> lfhv_checksum7 {
  { 1, 424586796ULL },
  { 10, 5111533890458962554ULL },
  { 100, 14384458847572363242ULL },
  { 1000, 15438438144713565350ULL },
  { 10000, 15267306927534078068ULL },
  { 100000, 3535642975512884138ULL },
  { 1000000, 8127813596529597114ULL },
  { 10000000, 2300420985112789568ULL },
  { 100000000, 3075781441988763947ULL },
  { 1000000000, 15168882957383223700ULL },
};

map<int, unsigned long long> lfhv_checksum8 {
  { 1, 815231912ULL },
  { 10, 9913565663100041150ULL },
  { 100, 15069744952245369793ULL },
  { 1000, 17912551117949039703ULL },
  { 10000, 396786536907475525ULL },
  { 100000, 1637257589636096178ULL },
  { 1000000, 339991071628281231ULL },
  { 10000000, 10939780323244126737ULL },
  { 100000000, 2805000507769535585ULL },
  { 1000000000, 4303073853386275283ULL },
};

map<int, map<int, unsigned long long>> skew_checksum7 {
  {1,
    {
    { 1, 424586796ULL },
    { 10, 4877741279594536584ULL },
    { 100, 1651525396633200272ULL },
    { 1000, 9361992673690044064ULL },
    { 10000, 9112746061879768640ULL },
    { 100000, 17238393772355088000ULL },
    { 1000000, 16521600776326779136ULL },
    { 10000000, 9870782629473069568ULL },
    { 100000000, 18060091962719123456ULL },
    { 1000000000, 18257314662972076032ULL },
    }
  }
};

map<int, map<int, unsigned long long>> skew_checksum8 {
  {1,
    {
    { 1, 815231912ULL },
    { 10, 9365553491223454448ULL },
    { 100, 8165526824854516448ULL },
    { 1000, 5393392765373772480ULL },
    { 10000, 15925377051034456960ULL },
    { 100000, 13777427989925339904ULL },
    { 1000000, 4474280651906801152ULL },
    { 10000000, 2911392226073029632ULL },
    { 100000000, 18060091962719123456ULL },
    { 1000000000, 18257314662972076032ULL },
    }
  }
};

struct Statistics {
  // Managed by test framework.
  int N, Q, selectivity, verified;
  string host;
  string algorithm;
  string workload;
  double insert_time;
  double query_time;
  unsigned long long checksum;

  // Each algorithm should fill in the following stats at the end of run:
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

  void print() {
    printf("%lu,", system_clock::to_time_t(system_clock::now()));
    printf("\"%s\",", host.c_str());
    printf("\"%s\",", algorithm.c_str());
    printf("\"%s\",", workload.c_str());
    printf("%d,", N);
    printf("%d,", Q);
    printf("%d,", selectivity);
    printf("%d,", verified);
    printf("%.6lf,", insert_time);
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

template<typename M>
bool checksum_match(M &m, int k, unsigned long long v) {
  if (m.count(k)) {
    if (m[k] == v) return true;
    fprintf(stderr, "\033[1;31mFAILED\033[0m checksum %llu != %llu\n", m[k], v);
  }
  return false;
}

template<typename M>
bool checksum_match2(M &m, int k1, int k2, unsigned long long v) {
  if (m.count(k1) && m[k1].count(k2)) {
    if (m[k1][k2] == v) return true;
    fprintf(stderr, "\033[1;31mFAILED\033[0m checksum %llu != %llu\n", m[k1][k2], v);
  }
  return false;
}

int verify(int U, Statistics &s) {
  switch (U) {
    case 0 : return s.N == 100000000 && checksum_match(noup_checksum, s.Q, s.checksum);
    case 1 : return (s.N == 10000000) ? checksum_match(lfhv_checksum7, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match(lfhv_checksum8, s.Q, s.checksum) : 0);
    case 2 : return (s.N == 10000000) ? checksum_match2(skew_checksum7, s.selectivity, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match2(skew_checksum8, s.selectivity, s.Q, s.checksum) : 0);
  }
  return false;
}

void init(int *arr, int N);  // Initializes the initial values of N integers.
void insert(int value);      // Inserts the value.
int query(int value);        // Query for lower bound, returns 0 if not found.
void erase(int value);       // Deletes the value. The value guaranteed to exists.
void results(Statistics &s);

Statistics s;

int main(int argc, char *argv[]) {
  char *prog = argv[0];
  while (true) {
    char *p = strstr(prog, "/");
    if (p) prog = p + 1; else break;
  }
  s.algorithm = prog;
  s.host = argv[1];
  s.N = atoi(argv[2]);
  int U = atoi(argv[3]);
  int MAX_Q = atoi(argv[4]);
  s.selectivity = s.N;
  if (argc > 5) s.selectivity = atoi(argv[5]);

  mt19937 gen(140384);
  uniform_int_distribution<> dis;
  uniform_int_distribution<> disN(0, s.N - 1);

  int *iarr = new int[s.N];
  for (int i = 0; i < s.N; i++) iarr[i] = dis(gen);

  s.insert_time = time_it([&] { init(iarr, s.N); });

  s.checksum = 0;
  s.query_time = 0;

  switch (U) {
    case 0 : s.workload = "NOUP"; break;
    case 1 : s.workload = "LFHV"; break;
    case 2 : {
      s.workload = "SKEW";
      auto mn = dis.min();
      auto mx = dis.max();
      // fprintf(stderr, "old query domain %d %d\n", mn, mx);
      long long sz = (long long) mx - mn;
      if (s.selectivity < sz) {
        mn = dis(gen) % (sz - s.selectivity);
        mx = mn + s.selectivity;
        dis = uniform_int_distribution<>(mn, mx);
        // fprintf(stderr, "new query domain %d %d\n", mn, mx);
      }
      break;
    }
    default : assert(0); abort(); break;
  }

  for (s.Q = 1; ; s.Q *= 10) {
    s.query_time += time_it([&] {
      int nQ = s.Q - s.Q / 10; // nQ = how many queries needed.
      for (int i = 0; i < nQ; i++) {
        s.checksum = s.checksum * 13 + query(dis(gen));
        if (U == 1 && i % 1000 == 0) {
          for (int j = 0; j < 1000; j++) {
            int k = disN(gen);
            erase(iarr[k]);
            insert(iarr[k] = dis(gen));
          }
        }
      }
    });

    results(s);
    s.verified = verify(U, s);

    s.print();
    fflush(stdout);
    if (s.Q >= MAX_Q) break;
  }
  return 0;
}
