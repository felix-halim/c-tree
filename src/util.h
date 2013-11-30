#include <map>

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)

template<typename F>
static double time_it(F f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

static char* parse_algorithm_name(char *prog) {
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
    for (long long j = 0; i >= 1000000000 && j < 100; j += 5) samples.push_back(i * j / 100);
  }
  sort(samples.begin(), samples.end());
  samples.erase(unique(samples.begin(), samples.end()), samples.end());
  return samples;
}

struct Statistics {
  string algorithm;
  string query_workload;
  string update_workload;

  Statistics(string a, string q, string u): algorithm(a), query_workload(q), update_workload(u) {
    n_index = 0;
    n_bytes = 0;
    n_slack_int = 0;
    n_slack_leaf = 0;
    n_internal = 0;
    n_leaf = 0;
    n_small = 0;
    n_large = 0;
    n_chained = 0;
    bt_int_sz = 0;
    bt_leaf_sz = 0;
    large_touch = 0;
    small_touch = 0;
    art_n4 = 0;
    art_n16 = 0;
    art_n48 = 0;
    art_n256 = 0;
  }

  // Each algorithm optionally fill in the following stats at the end of run:
  int N;
  int n_index;
  long long n_bytes;
  int n_slack_int;
  int n_slack_leaf;
  int n_internal;
  int n_leaf;
  int n_small;
  int n_large;
  int n_chained;
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
      "n_index,n_bytes,n_slack_int,n_slack_leaf,n_internal,n_leaf,n_small,n_large,n_chained,"
      "bt_int_sz,bt_leaf_sz,large_touch,small_touch,art_n4,art_n16,art_n48,art_n256,query_type\n");
  }

  void print(long long Q, double selectivity, int verified, double insert_time, double query_time, double update_time, unsigned long long checksum) {
    printf("%lu", system_clock::to_time_t(system_clock::now()));
    printf(",\"%s\"", algorithm.c_str());
    printf(",\"%s\"", query_workload.c_str());
    printf(",\"%s\"", update_workload.c_str());
    printf(",%d", N);
    printf(",%lld", Q);
    printf(",%lf", selectivity);
    printf(",\"%c\"", verified ? 'Y' : 'N');
    printf(",%.6lf", insert_time);
    printf(",%.6lf", update_time);
    printf(",%.6lf", query_time);
    printf(",%llu", checksum);
    printf(",%d", n_index);
    printf(",%lld", n_bytes);
    printf(",%d", n_slack_int);
    printf(",%d", n_slack_leaf);
    printf(",%d", n_internal);
    printf(",%d", n_leaf);
    printf(",%d", n_small);
    printf(",%d", n_large);
    printf(",%d", n_chained);
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

    fflush(stdout);
  }
};


std::map<int, unsigned long long> noup_checksum7 {
  { 1, 1443471823ULL },
  { 10, 15408736745479646270ULL },
  { 100, 18218735101919177500ULL },
  { 1000, 1492377536669543510ULL },
  { 10000, 6390096138010707052ULL },
  { 100000, 4310109292676179566ULL },
  { 1000000, 2905447459389835922ULL },
  { 10000000, 17617009477471388694ULL },
  { 100000000, 14227558382274216239ULL },
  { 1000000000, 4081601480252848911ULL },
};

std::map<int, unsigned long long> noup_checksum8 {
  { 1, 1443471377ULL },
  { 10, 15408731924548745167ULL },
  { 100, 15371389476403399536ULL },
  { 1000, 4508019715468061903ULL },
  { 10000, 2034706957392129817ULL },
  { 100000, 10400446070826407748ULL },
  { 1000000, 15530041155553518765ULL },
  { 10000000, 5555997414426845907ULL },
  { 100000000, 9188129546620505387ULL },
  { 1000000000, 9518872207859602874ULL },
};

std::map<int, unsigned long long> lfhv_checksum7 {
  { 1, 1443471823ULL },
  { 10, 15408736745479646270ULL },
  { 100, 18218735101919177500ULL },
  { 1000, 1492377536669543510ULL },
  { 10000, 13983177334309733995ULL },
  { 100000, 3039814775556676946ULL },
  { 1000000, 18126500374871495601ULL },
  { 10000000, 4264665950854348760ULL },
  { 100000000, 2042326541013504073ULL },
  { 1000000000, 7068591044984524264ULL },
};

std::map<int, unsigned long long> lfhv_checksum8 {
  { 1, 1443471377ULL },
  { 10, 15408731924548745167ULL },
  { 100, 15371389476403399536ULL },
  { 1000, 4508019715468061903ULL },
  { 10000, 2034706957392129817ULL },
  { 100000, 3273507535726318554ULL },
  { 1000000, 14322267694170196479ULL },
  { 10000000, 14719227275801663666ULL },
  { 100000000, 3798163431919307505ULL },
  { 1000000000, 17550971557407402583ULL },
};

std::map<int, unsigned long long> queue_checksum8 {
  { 1, 1443471381ULL },
  { 10, 16262347208527018447ULL },
  { 100, 8802644747984416ULL },
  { 1000, 1068309002083697ULL },
  { 10000, 13156602606304856011ULL },
  { 100000, 13749468600757118025ULL },
  { 1000000, 5838631702442076503ULL },
  { 10000000, 2296152760828406456ULL },
  { 100000000, 998426981192674219ULL },
  { 1000000000, 6070397024749655732ULL },
};

std::map<int, unsigned long long> trash_checksum8 {
  { 1, 1443471377ULL },
  { 10, 15408731924548745167ULL },
  { 100, 15371389476403399536ULL },
  { 1000, 4508019715468061903ULL },
  { 10000, 2034706957392129817ULL },
  { 100000, 13146551787786516041ULL },
  { 1000000, 8414686567344604143ULL },
  { 10000000, 12892445546246666539ULL },
  { 100000000, 6037015176804974254ULL },
  { 1000000000, 6859403095055459567ULL },
};

std::map<int, unsigned long long> delete_checksum8 {
  { 1, 1443471377ULL },
  { 10, 15408731924548745167ULL },
  { 100, 15371389476403399536ULL },
  { 1000, 4508019715468061903ULL },
  { 10000, 2034706957392129817ULL },
  { 100000, 18091595545939233641ULL },
  { 1000000, 18088951063046439492ULL },
  { 10000000, 1165036610406096989ULL },
  { 100000000, 7302612522330782360ULL },
};

std::map<int, unsigned long long> append_checksum8 {
  { 1, 154442717ULL },
  { 10, 1775014986347233215ULL },
  { 100, 7343837289042995878ULL },
  { 1000, 9361056522399319518ULL },
  { 10000, 6123724530754144970ULL },
  { 100000, 776250836994303551ULL },
  { 1000000, 18295530010266777255ULL },
  { 10000000, 6219452757308234767ULL },
  { 100000000, 18286807892613052455ULL },
  { 1000000000, 16340975687709899762ULL },
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

static int verify(int U, int N, int Q, unsigned long long checksum) {
  switch (U) {
    case 0 : return (N == 10000000) ? checksum_match(noup_checksum7, Q, checksum) :
                    (N == 100000000 ? checksum_match(noup_checksum8, Q, checksum) : 0);
    case 1 : return (N == 10000000) ? checksum_match(lfhv_checksum7, Q, checksum) :
                    (N == 100000000 ? checksum_match(lfhv_checksum8, Q, checksum) : 0);
    case 3 : return N == 100000000 ? checksum_match(queue_checksum8, Q, checksum) : 0;
    case 4 : return N == 100000000 ? checksum_match(trash_checksum8, Q, checksum) : 0;
    case 5 : return N == 100000000 ? checksum_match(delete_checksum8, Q, checksum) : 0;
    case 6 : return N == 100000000 ? checksum_match(append_checksum8, Q, checksum) : 0;
  }
  return false;
}
