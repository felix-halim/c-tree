#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <random>

using namespace std;

const char* update_workload[10] = {
  "NOUP",   // 0. Read only queries.
  "LFHV",   // 1. Update 1000 tuples every 1000 queries.
  "HFLV",   // 2. Update 10 tuples every 10 queries.
  "QUEUE",  // 3. Remove the largest value, insert new smaller value than any existing value.
  "TRASH",  // 4. Insert values in the middle of the domain.
  "DELETE", // 5. Delete 1000 tuples every 1000 queries.
  "APPENDSKY", // 6. Insert 100K tuples every query. 
  "APPEND", // 7. Insert 100K tuples every query. 
  "APPEND", // 8. Insert 10 tuples every 10 query. 
  "SKEW",   // 9. LFHV but on 20% domain only.
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

    if (U != 6) {
      load();
      if (U == 3) prepare_queue();
      // else if (U == 7) N = 100000;
      // else if (U == 8) N = 10;
    } else {
      load(100000);
    }

    if (U == 5) {
      prepare_deletion(0);
      // MAXQ = min((long long) MAXQ);
    }
  }

  void prepare_deletion(int N) { shuffle(arr.begin(), arr.begin() + N, gen); }
  unsigned* get_arr() { return &arr[0]; }
  int get_n() { return min(1000000000, (W == 0) ? arr.size() : (arr.size() / 2)); }
  int size() { return arr.size(); }
  void clear() { arr.clear(); }

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
      // if (arr[i] <= next_smallest) {
        arr[i] += next_smallest;
      // }
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

  string name() {
    return update_workload[U];
  }

  int code() { return U; }

  double execute(long long i, std::function<void(unsigned)> insert, std::function<void(unsigned)> erase) {
    switch (U) {
      // NOUP.
      case 0: return 0;

      // LFHV.
      case 1: return (i % 1000) ? 0 :
        time_it([&] {
          unsigned a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        });

      // HFHV.
      case 2: return (i % 10) ? 0 :
        time_it([&] {
          unsigned a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        });

      // QUEUE.
      case 3: return (i % 10) ? 0 :
        time_it([&] {
          unsigned a, b;
          REP(j, 10) {
            update_queue(a, b);
            erase(a);
            insert(b);
          }
        });

      // TRASH.
      case 4: return (i != 10000) ? 0 :
        time_it([&] {
          unsigned *arr = get_arr();
          int N = get_n();
          REP(j, 1000000) {
            insert(arr[N + j]);
            // fprintf(stderr, "%d \n", arr[N + j]);
          }
        });

      // DELETE.
      case 5: return (i % 1000) ? 0 :
        time_it([&] {
          REP(j, 1000) {
            erase(update_delete());
          }
        });

      // APPEND SKY SERVER.
      case 6: return
        time_it([&] {
          // if (MAXQ != -1) {
            clear();
            bool loaded = false;
            // load_time += time_it([&] {
              loaded = load(100000);
              // query_w.set_max(max_element());
              // N += size();
            // });
            if (loaded) {
              unsigned *arr = get_arr();
              REP(j, size()) insert(arr[j]);
            } else {
              // MAXQ = -1;
            }
          // }
        });

      // APPEND.
      case 7: return
        time_it([&] {
          // if (MAXQ != -1) {
            int N = get_n();
            int add = size() / 2 - N;
            if (add > 0) {
              unsigned *arr = get_arr();
              REP(j, min(add, 100000)) insert(arr[N++]);
            } else {
              // MAXQ = -1;
            }
          // }
        });

      // APPEND.
      case 8: return (i % 10) ? 0 :
        time_it([&] {
          int N = get_n();
          // if (MAXQ != -1) {
            int add = size() / 2 - N;
            if (add > 0) {
              unsigned *arr = get_arr();
              REP(j, min(add, 10)) insert(arr[N++]);
            } else {
              // MAXQ = -1;
            }
          // }
        });

      // SKEW.
      case 9: return (i % 1000) ? 0 :
        time_it([&] {
          unsigned a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        });

      default: return 0;
    }
  }
};
