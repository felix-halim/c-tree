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
    for (long long j = 0; j < 100; j += 5) samples.push_back(i * j / 100);
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
  int n_index;
  int n_bytes;
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

  void print(int N, long long Q, double selectivity, int verified, double insert_time, double query_time, double update_time, unsigned long long checksum) {
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



std::map<int, std::map<int, unsigned long long>> skew_checksum7 {
  {1,
    {
    { 1, 1010000016ULL },
    { 10, 11603089914351322464ULL },
    { 100, 4707422041215887680ULL },
    { 1000, 4813038930404913512ULL },
    { 10000, 18165742350779668520ULL },
    { 100000, 6245446159144308328ULL },
    { 1000000, 16422320693352767464ULL },
    { 10000000, 12268924908534872808ULL },
    { 100000000, 3876364641539258600ULL },
    { 1000000000, 2367851660476687592ULL },
    }
  }
};

std::map<int, std::map<int, unsigned long long>> skew_checksum8 {
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
    // case 2 : return (N == 10000000) ? checksum_match2(skew_checksum7, selectivity, Q, checksum) :
    //                 (N == 100000000 ? checksum_match2(skew_checksum8, selectivity, Q, checksum) : 0);
    case 3 : return N == 100000000 ? checksum_match(queue_checksum8, Q, checksum) : 0;
    case 4 : return N == 100000000 ? checksum_match(trash_checksum8, Q, checksum) : 0;
    case 5 : return N == 100000000 ? checksum_match(delete_checksum8, Q, checksum) : 0;
    case 6 : return N == 100000000 ? checksum_match(append_checksum8, Q, checksum) : 0;
  }
  return false;
}
