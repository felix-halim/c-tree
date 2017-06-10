// TODO: break this down into separate files.
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <random>

using namespace std;

const char* update_workload[10] = {
  "NOUP",    // 0. Read only queries.
  "LFHV",    // 1. Update 1000 tuples every 1000 queries.
  "HFLV",    // 2. Update 10 tuples every 10 queries.
  "QUEUE",   // 3. Remove the largest value, insert new smaller value than any existing value.
  "TRASH",   // 4. Insert values in the middle of the domain.
  "DELETE",  // 5. Delete 1000 tuples every 1000 queries.
  "APPENDB", // 6. Insert 100K tuples every query. 
  "SKEW",    // 7. LFHV but on 20% domain only.
};

class Update {
  vector<unsigned> arr;
  FILE *in;
  unsigned next_smallest;
  unsigned max_value;
  int U, W;
  std::function<void(unsigned)> update_max_cb;

  mt19937 gen;
  uniform_int_distribution<> dis;

 public:

  Update(char *fn, int u, int w, std::function<void(unsigned)> on_update_max):
      max_value(0), U(u), W(w), update_max_cb(on_update_max), gen(81188) {

    next_smallest = 1080000000;
    in = fopen(fn, "rb");
    if (!in) { fprintf(stderr,"Error opening file %s\n", fn); exit(1); }

    if (U == 6) {
      load(100000);
    } else {
      load();
      if (U == 3) prepare_queue();
      // else if (U == 7) N = 100000;
      // else if (U == 8) N = 10;
    }

    if (U == 5) {
      prepare_deletion(0);
      // MAXQ = min((long long) MAXQ);
    }
  }

  unsigned* get_arr() { return &arr[0]; }
  int get_n() { return min(1000000000, int((W == 0) ? arr.size() : (arr.size() / 2))); }
  int size() { return arr.size(); }
  void clear() { arr.clear(); }
  void prepare_deletion(int N) { shuffle(arr.begin(), arr.begin() + N, gen); }
  unsigned get_next_smallest() { return next_smallest; }

  bool load(int amt = 1000000000) {
    if (!in) return false; // Already loaded.
    unsigned tmp[1024] = { 0 };
    while (!feof(in)) {
      int need = min(1024, amt - size());
      if (need <= 0) break;
      int N = fread(tmp, sizeof(unsigned), need, in);
      for (int i = 0; i < N; i++) arr.push_back(tmp[i]);
      max_value = max(max_value, *std::max_element(tmp, tmp + N));
    }
    if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); exit(1); }
    if (feof(in)) { fclose(in); in = NULL; }
    update_max_cb(max_value);
    return true;
  }

  void prepare_queue() {
    int N = get_n();
    sort(arr.begin(), arr.end());
    for (int i = 0; i < N; i++) {
      if (arr[i] > next_smallest)
        fprintf(stderr, "%d %d\n", i, arr[i]);
      assert(arr[i] <= next_smallest);
      arr[i] += next_smallest;
    }
  }

  void update_queue(unsigned &to_del, unsigned &to_add) {
    static int qidx = -1;
    if (qidx < 0) qidx = get_n() - 1;
  // if ((i + 1) % 1000000 == 0)
  //   fprintf(stderr, "arr[%d] = %d, next = %d\n", qidx,arr[qidx],next_smallest);
    to_del = arr[qidx];
    to_add = next_smallest;
    arr[qidx--] = next_smallest--;
    assert(next_smallest >= 0);
  }

  int update_delete() {
    static int lastN = get_n();
    assert(lastN > 0);
    return arr[--lastN];
  }

  void update(unsigned &to_del, unsigned &to_add) {
    int N = get_n();
    int k = dis(gen) % N;
    to_del = arr[k];
    int l = N + dis(gen) % N;
    assert(l < size());
    swap(arr[k], arr[l]);
    to_add = arr[k];
  }

  void update_skew(unsigned &to_del, unsigned &to_add, int start, int end) {
    int N = end - start;
    int k = start + dis(gen) % N;
    to_del = arr[k];
    int l = N + dis(gen) % N;
    assert(l < size());
    swap(arr[k], arr[l]);
    to_add = arr[k];
  }

  string name() {
    return update_workload[U];
  }

  int code() { return U; }

  void execute(long long i,
      std::function<void(unsigned)> insert,
      std::function<void(unsigned)> erase,
      std::function<void(double)> load_time_cb,
      std::function<void(double)> update_time_cb) {

    switch (U) {
      // NOUP.
      case 0: break;

      // LFHV.
      case 1: if (i % 1000 == 0)
        update_time_cb(time_it([&] {
          unsigned a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        }));
        break;

      // HFHV.
      case 2: if (i % 10 == 0)
        update_time_cb(time_it([&] {
          unsigned a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        }));
        break;

      // QUEUE.
      case 3: if (i % 10 == 0)
        update_time_cb(time_it([&] {
          unsigned a, b;
          REP(j, 10) {
            update_queue(a, b);
            erase(a);
            insert(b);
          }
        }));
        break;

      // TRASH.
      case 4: if (i == 10000)
        update_time_cb(time_it([&] {
          unsigned *arr = get_arr();
          int N = get_n();
          REP(j, 1000000) {
            insert(arr[N + j]);
            // fprintf(stderr, "%d \n", arr[N + j]);
          }
        }));
        break;

      // DELETE.
      case 5: if (i % 1000 == 0)
        update_time_cb(time_it([&] {
          REP(j, 1000) {
            erase(update_delete());
          }
        }));
        break;

      // APPEND SKY SERVER.
      case 6:
        if (arr.size()) {
          clear();
          load_time_cb(time_it([&] { load(100000); }));
          if (arr.size()) {
            update_time_cb(time_it([&] {
              unsigned *arr = get_arr();
              REP(j, size()) insert(arr[j]);
            }));
          }
        }
        break;

      // SKEW.
      case 7:
        update_time_cb(time_it([&] {
          unsigned a, b;
          if (i < 2000) {
            REP(j, 1000) {
              update(a, b);
              erase(a);
              insert(b);
            }
          }
        }));
        break;

      default: break;
    }
  }
};


#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>
#include <random>

#ifndef SKEW_X
  #define SKEW_X 0
  #define SKEW_Y 0
#endif

using namespace std;

const char* workload_name[18] = {
  "SkyServer",     // 0
  "Random",        // 1
  "SeqOver",       // 2
  "SeqInv",        // 3
  "SeqRand",       // 4
  "SeqNoOver",     // 5
  "SeqAlt",        // 6
  "ConsRandom",    // 7
  "ZoomIn",        // 8
  "ZoomOut",       // 9
  "SeqZoomIn",     // 10
  "SeqZoomOut",    // 11
  "Skew",          // 12
  "ZoomOutAlt",    // 13
  "SkewZoomOutAlt",// 14
  "Periodic",      // 15
  "Mixed",         // 16
  "Domain"         // 17
};

class Workload {
  unsigned N;    // The number of elements in arr.
  int W;    // The selected workload to be generated.
  unsigned S;    // The selectivity (unused for some workloads).
  int I;    // The I'th query (internal use only).
  unsigned a, b; // The last query range [a,b].
  unsigned seq_jump;
  long long MAXQ;
  double s;
  vector<unsigned> skyq;
  int skyqi;

  mt19937 gen;
  uniform_int_distribution<> dis;

  int nextInt(int modulo) {
    return dis(gen) % modulo;
  }

  void init_skyq() {
    FILE *in = fopen("data/skyserver.queries", "r");
    if (!in) {
      fprintf(stderr, "Fail loading file data/skyserver.queries\n");
    } else {
      double x, y;
      while (fscanf(in, "%lf %lf", &x, &y) != EOF) {
        unsigned z = unsigned(y * 1000000 * 10);
        skyq.push_back(z);
      }
      fclose(in);
    }
  }

  // Query workload based on the predefined queries from file.
  bool skyserver_w() {
    a = skyq[skyqi++];
    b = a + S;
    if (skyqi >= (int) skyq.size()) {
      skyqi -= skyq.size();
      return false;
    }
    return true;
  }

  // a and b is selected uniformly at random and a < b
  bool random_w() {
    if (S == 0) { // random selectivity
      do {
        a = nextInt(N);
        b = nextInt(N);
        if (a > b) swap(a, b);
      } while (a == b);
    } else {
      a = nextInt(N - S + 1);
      b = max(a + 1, a + S);
    }
    return true;
  }
  
  // a will be incremented by 10 every subsequent query.
  // The range may overlap with the next queries.
  bool seq_over_w() {
    a = 10 + I * seq_jump;
    if (a + 5 > N) return false;
    if (S == 0) {
      b = a + nextInt(N - a) + 1;
    } else {
      b = a + S;
    }
    return b <= N;
  }
  
  // The opposite direction from the seq_over_w.
  bool seq_inv_w() {
    if (!seq_over_w()) return false;
    a = N - a;
    b = N - b;
    swap(a, b);
    return true;
  }
  
  // Half the time is seq_over half the time is random_w.
  bool seq_rand_w(){
    return (I & 1) ? seq_over_w() : random_w();
  }
  
  // Sequential with no overlap with the subsequent query ranges.
  bool seq_no_over_w() {
    static unsigned prevB;
    if (!I) prevB = 0;
    a = prevB + 10;
    if (a + 5 > N) return false;
    if (S == 0) {
      b = a + nextInt(N - a) + 1;
    } else {
      b = a + S;
    }
    prevB = b;
    return b <= N;
  }
  
  // Sequential alternate at the beginning and end.
  bool seq_alt_w() {
    return (I & 1) ? seq_over_w() : seq_inv_w();
  }
  
  // Pick 1000 integers and produce range queries with endpoints using the 1000 picked integers.
  bool cons_rand_w() {
    static unsigned R[1000];
    if (!I) {
      for (int i = 0; i < 1000; i++) {
        R[i] = nextInt(N);
      }
    }
    do {
      a = R[nextInt(1000)];
      b = R[nextInt(1000)];
    } while (a == b);
    if (a > b) swap(a,b);
    return true;
  }
  
  // Start at the [middle - 100500, middle + 100500), 
  // then zoom in by making the query range smaller by 100 on each query.
  bool zoom_in_w() {
    static unsigned L; if (!I) L = N / 3;
    static unsigned R; if (!I) R = 2 * N / 3;
    if (L >= R || R > N) return false;
    a = L; L += 100;  // make the range smaller
    b = R; R -= 100;
    return true;
  }
  
  // Start at the [middle - 500, middle + 500),
  // then zoom out by making the query range larger by 100 each query.
  bool zoom_out_w() {
    static unsigned L; if (!I) L = N / 2 - 500;
    static unsigned R; if (!I) R = N / 2 + 500;
    if (L < 1 || R > N) return false;
    a = L; L -= 100;  // make the range bigger
    b = R; R += 100;
    return true;
  }
  
  // After zooming in on one region, move to next unexplored region to the right.
  bool seq_zoom_in() {
    static unsigned L; if (!I) L = 1;
    static unsigned G = 100000;
    static unsigned R; if (!I) R = G;
    if (L >= R) L += G, R = L + G;
    if (R > N) return false;
    a = L; L += 100;
    b = R; R -= 100;
    return true;
  }
  
  // After zooming out on one ragion, move to the next unexplored region on the right.
  bool seq_zoom_out() {
    static unsigned G = 100000;
    static unsigned L; if (!I) L = G / 2 + 1000;
    static unsigned R; if (!I) R = L + 10;
    if (R > L + G) {
      L = R + G / 2 + 1000;
      R = L + 10;
    }
    if (R > N) return false;
    a = L; L -= 100;
    b = R; R += 100;
    return true;
  }
  
  // X percent of the queries falls within Y percent of the value range and
  // (100-X) percent of the queries falls within the rest of the value range.
  bool skew_w(int X, int Y) {
    int nY = ((long long) N) * Y / 100;
    if (nextInt(100 * 1000000) < X * 1000000) {
      do {
        a = nextInt(nY);
        b = nextInt(nY);
      } while (a == b);
      // fprintf(stderr, "nY1 = %10d, %10d %10d\n", nY, a, b);
    } else {
      do {
        a = nY + nextInt(N - nY);
        b = nY + nextInt(N - nY);
      } while (a == b);
      // fprintf(stderr, "nY2 = %10d, %10d %10d, N = %d\n", nY, a, b, N);
    }
    if (a > b) swap(a, b);
    return true;
  }

  // Start at the [middle - 500, middle + 500),
  // then zoom out by making the query range larger by 100 each query.
  bool zoom_out_alt_w() {
    static unsigned L; if (!I) L = N / 2 - 500;
    static unsigned R; if (!I) R = N / 2 + 500;
    if (L < 1 || R > N) return false;
    if (I & 1) {
      a = L; 
      L -= 100;  // make the range bigger
      b = a + 10;
    } else {
      b = R; 
      R += 100;
      a = b - 10;
    }
    return true;
  }
  
  // Start at the [middle - 500, middle + 500),
  // then zoom out by making the query range larger by 100 each query.
  bool skew_zoom_out_alt_w() {
    static unsigned L; if (!I) L = N - 355000;
    static unsigned R; if (!I) R = N - 350000;
    if (L < 1 || R > N) return false;
    if (I & 1) {
      b = R; 
      R += 20;
      a = b - 10;
    } else {
      a = L; 
      L -= 20;  // make the range bigger
      b = a + 10;
    }
    return true;
  }
  
  bool periodic_w() {
    static long long jump = 1000001;
    a = (I * jump) % N;
    b = a + 10;
    return true;
  }

  bool mixed_w() {
    static unsigned work = 0;
    static unsigned base = 0;
    if (I % 1000 == 0) {
      work = nextInt(15) + 1;
      base = nextInt(20);
    }
    unsigned tW = W; W = work;
    unsigned tI = I; I %= 1000;
    unsigned tN = N; N /= 20;
    unsigned ta, tb;
    bool ok = query(ta, tb);
    W = tW;
    I = tI;
    if (!ok){
      N = tN;
      work = nextInt(15) + 1;
      return mixed_w();
    }
    a = ta + base * N;
    b = tb + base * N;
    N = tN;
    return true;
  }

public : 

  Workload(long long nQ, double sel, int w): N(0), W(w), S(0), I(0), a(0), b(0), seq_jump(20), MAXQ(nQ), s(sel), gen(140384) {
    if (W < 0 || W >= 18) {
      fprintf(stderr,"Workload number %d is not found!\n", W);
      exit(1);
    }
    if (W == 0) init_skyq();

    if (W == 2) {
      int j = N / 100 / MAXQ;
      fprintf(stderr, "j = %d\n", j);
      set_seq_jump(max(1, j));
    } else if (W == 17) {
      set_max(10);
    }
  }

  long long maxq() { return MAXQ; }

  void set_max(unsigned mx) {
    N = mx;
    S = mx * s;
  }

  void inc_max() { N++; }

  void set_seq_jump(unsigned seq_jump_) {
    seq_jump = seq_jump_;
  }

  bool query(unsigned &na, unsigned &nb) {
    assert(N > 0);
    switch (W) {
      case 0 : if (!skyserver_w()) return false; break;
      case 1 : if (!random_w()) return false; break;
      case 2 : if (!seq_over_w()) return false; break;
      case 3 : if (!seq_inv_w()) return false; break;
      case 4 : if (!seq_rand_w()) return false; break;
      case 5 : if (!seq_no_over_w()) return false; break;
      case 6 : if (!seq_alt_w()) return false; break;
      case 7 : if (!cons_rand_w()) return false; break;
      case 8 : if (!zoom_in_w()) return false; break;
      case 9 : if (!zoom_out_w()) return false; break;
      case 10 : if (!seq_zoom_in()) return false; break;
      case 11 : if (!seq_zoom_out()) return false; break;
      case 12 : if (!skew_w(SKEW_X, SKEW_Y)) return false; break;
      case 13 : if (!zoom_out_alt_w()) return false; break;
      case 14 : if (!skew_zoom_out_alt_w()) return false; break;
      case 15 : if (!periodic_w()) return false; break;
      case 16 : if (!mixed_w()) return false; break;
      case 17 : if (!random_w()) return false; break;
      default : fprintf(stderr, "Invalid W = %d\n", W); abort();
    }
    na = a; nb = b; I++;
    return true;
  }

  string name() { return workload_name[W]; }

  double selectivity() { return s; }
};


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

    // fprintf(stderr, "\033[1;31mFAILED\033[0m checksum %llu != %llu\n", m[k], v);
