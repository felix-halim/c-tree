#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <random>

using namespace std;

class Update {
  vector<int> arr;
  FILE *in;
  int next_smallest;
  int max_value;
  int U, W, N;
  std::function<void(int)> update_max_cb;

  mt19937 gen;
  uniform_int_distribution<> dis;

  string update_workload[10] {
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

 public:

  Update(char *fn, int u, int w, std::function<void(int)> on_update_max):
      max_value(0), U(u), W(w), update_max_cb(on_update_max), gen(81188) {

    next_smallest = 1080000000;
    in = fopen(fn, "rb");
    if (!in) { fprintf(stderr,"Error opening file %s\n", fn); exit(1); }

    if (U != 6) {
      load();
      N = size();
      if (U == 3) prepare_queue();
      else if (U == 7) N = 100000;
      else if (U == 8) N = 10;
    } else {
      load(100000);
      N = size();
    }

    if (U == 5) {
      prepare_deletion(N);
      // MAXQ = min((long long) MAXQ);
    }
  }

  void prepare_deletion(int N) { shuffle(arr.begin(), arr.begin() + N, gen); }
  int max_element() { return max_value; }
  int* get_arr() { return &arr[0]; }
  int get_n() {
    fprintf(stderr, "W = %d\n", W);
   return (W == 0) ? N : (N / 2); }
  int size() { return arr.size(); }
  void clear() { arr.clear(); }

  int get_next_smallest() { return next_smallest; }

  bool load(int amt = 1000000000) {
    if (!in) return false; // Already loaded.
    int tmp[1024] = { 0 };
    while (!feof(in)) {
      int need = min(1024, amt - size());
      if (need <= 0) break;
      int N = fread(tmp, sizeof(int), need, in);
      for (int i = 0; i < N; i++) arr.push_back(tmp[i]);
      max_value = max(max_value, *std::max_element(tmp, tmp + N));
    }
    if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); exit(1); }
    if (feof(in)) { fclose(in); in = NULL; }
    update_max_cb(max_value);
    return true;
  }

  void prepare_queue() {
    int N = size();
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

  void update_queue(int &to_del, int &to_add) {
    static int qidx = -1;
    if (qidx < 0) qidx = N - 1;
  // if ((i + 1) % 1000000 == 0)
  //   fprintf(stderr, "arr[%d] = %d, next = %d\n", qidx,arr[qidx],next_smallest);
    to_del = arr[qidx];
    to_add = next_smallest;
    arr[qidx--] = next_smallest--;
    assert(next_smallest >= 0);
  }

  int update_delete() {
    static int lastN = N;
    assert(lastN > 0);
    return arr[--lastN];
  }

  void update(int &to_del, int &to_add) {
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

  double execute(long long i, std::function<void(int)> insert, std::function<void(int)> erase) {
    switch (U) {
      // NOUP.
      case 0: return 0;

      // LFHV.
      case 1: return (i % 1000) ? 0 :
        time_it([&] {
          int a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        });

      // HFHV.
      case 2: return (i % 10) ? 0 :
        time_it([&] {
          int a, b;
          REP(j, 1000) {
            update(a, b);
            erase(a);
            insert(b);
          }
        });

      // QUEUE.
      case 3: return (i % 10) ? 0 :
        time_it([&] {
          int a, b;
          REP(j, 10) {
            update_queue(a, b);
            erase(a);
            insert(b);
          }
        });

      // TRASH.
      case 4: return (i != 10000) ? 0 :
        time_it([&] {
          int *arr = get_arr();
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
              N += size();
            // });
            if (loaded) {
              int *arr = get_arr();
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
            int add = size() / 2 - N;
            if (add > 0) {
              int *arr = get_arr();
              REP(j, min(add, 100000)) insert(arr[N++]);
            } else {
              // MAXQ = -1;
            }
          // }
        });

      // APPEND.
      case 8: return (i % 10) ? 0 :
        time_it([&] {
          // if (MAXQ != -1) {
            int add = size() / 2 - N;
            if (add > 0) {
              int *arr = get_arr();
              REP(j, min(add, 10)) insert(arr[N++]);
            } else {
              // MAXQ = -1;
            }
          // }
        });

      // SKEW.
      case 9: return (i % 1000) ? 0 :
        time_it([&] {
          int a, b;
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
