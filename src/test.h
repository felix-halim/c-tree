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

map<int, unsigned long long> noup_checksum7 {
  { 1, 911793201573294536ULL },
  { 10, 6420604009734469797ULL },
  { 100, 7085843312534587758ULL },
  { 1000, 9956582071699372241ULL },
  { 10000, 8518346662228332032ULL },
  { 100000, 150630346973349967ULL },
  { 1000000, 14815815474562926180ULL },
  { 10000000, 17827943117006321667ULL },
  { 100000000, 17687247983157514810ULL },
  { 1000000000, 9738629927182098908ULL },
};

map<int, unsigned long long> noup_checksum8 {
  { 1, 1750697200359773610ULL },
  { 10, 12438330354614217067ULL },
  { 100, 14741302341397222675ULL },
  { 1000, 9425395587432092148ULL },
  { 10000, 11549427510549084983ULL },
  { 100000, 12545465948791404173ULL },
  { 1000000, 15637276492067202780ULL },
  { 10000000, 17170414591296787493ULL },
  { 100000000, 16844025207510256064ULL },
  { 1000000000, 13738747866145218018ULL },
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
  { 1, 1750697200359773610ULL },
  { 10, 12438330354614217067ULL },
  { 100, 14741302341397222675ULL },
  { 1000, 9425395587432092148ULL },
  { 10000, 16012959058611103843ULL },
  { 100000, 8179383034432303941ULL },
  { 1000000, 14629639428699753311ULL },
  { 10000000, 900779849998059449ULL },
  { 100000000, 11699558442089270172ULL },
  { 1000000000, 8697580139090060052ULL },
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
  },
  {10,
    {
    { 1, 815231912ULL },
    { 10, 9365553494349107804ULL },
    { 100, 8171037487378421780ULL },
    { 1000, 17169146298459277120ULL },
    { 10000, 2014747129663233980ULL },
    { 100000, 2598384795719555248ULL },
    { 1000000, 9819475034527542770ULL },
    { 10000000, 8692120839264508914ULL },
    { 100000000, 15620951258232046912ULL },
    { 1000000000, 700300658117644334ULL },
    }
  },
  {100,
    {
    { 1, 815231958ULL },
    { 10, 9365554019691036020ULL },
    { 100, 10922412009247266319ULL },
    { 1000, 11641928627107815566ULL },
    { 10000, 1199052141464273911ULL },
    { 100000, 2555351741377885890ULL },
    { 1000000, 3292500748377450424ULL },
    { 10000000, 17963198292074214506ULL },
    { 100000000, 16008672679042555171ULL },
    { 1000000000, 13228431060445324144ULL },
    }
  },
  {1000,
    {
    { 1, 815231958ULL },
    { 10, 9365554049909121970ULL },
    { 100, 16181518462155005893ULL },
    { 1000, 17166959429559299580ULL },
    { 10000, 3688150156723210678ULL },
    { 100000, 16889316009393866162ULL },
    { 1000000, 10967513538785242393ULL },
    { 10000000, 11349239286052037906ULL },
    { 100000000, 4531546032060817778ULL },
    { 1000000000, 7507862803611855209ULL },
    }
  },
  {10000,
    {
    { 1, 815232401ULL },
    { 10, 9365559338588609635ULL },
    { 100, 14994519453414616773ULL },
    { 1000, 14308756366165237384ULL },
    { 10000, 12023421841868224274ULL },
    { 100000, 11028320610099621192ULL },
    { 1000000, 12651980267904445848ULL },
    { 10000000, 7499569611410507432ULL },
    { 100000000, 5204245453749000211ULL },
    { 1000000000, 17973740366925032745ULL },
    }
  },
  {100000,
    {
    { 1, 815236782ULL },
    { 10, 9365611741620622532ULL },
    { 100, 16719980412430574505ULL },
    { 1000, 12502569532130943137ULL },
    { 10000, 4316617403387906955ULL },
    { 100000, 14996592735238995969ULL },
    { 1000000, 15869782845042943771ULL },
    { 10000000, 7681052631306348830ULL },
    { 100000000, 6985867537339164294ULL },
    { 1000000000, 14600947350663234873ULL },
    }
  },
  {1000000,
    {
    { 1, 815280591ULL },
    { 10, 9366135791235742197ULL },
    { 100, 14046875718249898173ULL },
    { 1000, 14234791087568095958ULL },
    { 10000, 7976510585028221979ULL },
    { 100000, 13078926292340494041ULL },
    { 1000000, 3969019173229452056ULL },
    { 10000000, 8188525352061995664ULL },
    { 100000000, 14273555022628740132ULL },
    { 1000000000, 10102588147120430653ULL },
    }
  },
  {10000000,
    {
    { 1, 815719012ULL },
    { 10, 9371380407504620367ULL },
    { 100, 11950101575393052230ULL },
    { 1000, 14332957116254203224ULL },
    { 10000, 15051951751392043501ULL },
    { 100000, 10311221489377729508ULL },
    { 1000000, 11879219536363063655ULL },
    { 10000000, 10116243532420672290ULL },
    { 100000000, 813609005315465710ULL },
    { 1000000000, 6052198173317751238ULL },
    }
  },
  {100000000,
    {
    { 1, 820207241ULL },
    { 10, 9425070719586968877ULL },
    { 100, 3433323067073822481ULL },
    { 1000, 9450543209421737376ULL },
    { 10000, 85233614717075882ULL },
    { 100000, 6582638904451919206ULL },
    { 1000000, 569479827167296108ULL },
    { 10000000, 16772647736709083161ULL },
    { 100000000, 2644058508593816476ULL },
    { 1000000000, 5086263021862614202ULL },
    }
  },
};

struct Statistics {
  // Managed by test framework.
  int N, Q, selectivity, verified;
  string host;
  string algorithm;
  string workload;
  double insert_time;
  double query_time;
  double update_time;
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
    case 0 : return (s.N == 10000000) ? checksum_match(noup_checksum7, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match(noup_checksum8, s.Q, s.checksum) : 0);
    case 1 : return (s.N == 10000000) ? checksum_match(lfhv_checksum7, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match(lfhv_checksum8, s.Q, s.checksum) : 0);
    case 2 : return (s.N == 10000000) ? checksum_match2(skew_checksum7, s.selectivity, s.Q, s.checksum) :
                    (s.N == 100000000 ? checksum_match2(skew_checksum8, s.selectivity, s.Q, s.checksum) : 0);
  }
  return false;
}

void init(long long *arr, int N);  // Initializes the initial values of N integers.
void insert(long long value);      // Inserts the value.
long long query(long long value);  // Query for lower bound, returns 0 if not found.
void erase(long long value);       // Deletes the value. The value guaranteed to exists.
void results(Statistics &s);

mt19937 gen(140384);
uniform_int_distribution<> dis;
long long dis_suffix = 0;
Statistics s;

inline static long long next_ll() {
  return dis(gen);
  // return (((long long) dis(gen)) << 31) | (dis_suffix++);
}

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

  long long *iarr = new long long[s.N * 2];
  int gap = 2147483647 / s.N / 2;
  iarr[0] = 1;
  long long mx = iarr[0];
  for (int i = 1; i < s.N * 2; i++) {
    iarr[i] = iarr[i - 1] + (next_ll() % gap) + 1;
    mx = max(mx, iarr[i]);
  }
  dis = uniform_int_distribution<>(iarr[0], mx + 1);
  random_shuffle(iarr, iarr + s.N * 2);

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

  uniform_int_distribution<> disN(0, s.N - 1);
  for (s.Q = 1; ; s.Q *= 10) {
    double update_time = 0;
    s.query_time += time_it([&] {
      int nQ = s.Q - s.Q / 10; // nQ = how many queries needed.
      for (int i = 0; i < nQ; i++) {
        s.checksum = s.checksum * 13 + query(next_ll());
        if (U == 1 && i && i % 1000 == 0) {
          update_time += time_it([&] {
            for (int j = 0; j < 1000; j++) {
              int k = disN(gen);
              // if (iarr[k] == -1) { j--; continue; }
              // erase(iarr[k]);
              // iarr[k] = -1;
              int l = s.N + disN(gen);
              // swap(iarr[k], iarr[l]);
              // insert(iarr[k]);
              if (iarr[l] == -1) { j--; continue; }
              insert(iarr[l]);
              iarr[l] = -1;
            }
          });
        }
      }
    });
    s.query_time -= update_time;
    s.update_time += update_time;

    results(s);
    s.verified = verify(U, s);

    s.print();
    fflush(stdout);
    if (s.Q >= MAX_Q) break;
  }
  return 0;
}
