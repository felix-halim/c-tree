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
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:

  Update(char *fn): max_value(0), gen(81188) {
    next_smallest = 1080000000;
    in = fopen(fn, "rb");
    if (!in) { fprintf(stderr,"Error opening file %s\n", fn); exit(1); }
  }

  void prepare_deletion(int N) { shuffle(arr.begin(), arr.begin() + N, gen); }
  int max_element() { return max_value; }
  int* get_arr() { return &arr[0]; }
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
    return true;
  }

  void prepare_queue(int N) {
    assert(size() >= N);
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

  void update_queue(int N, int &to_del, int &to_add) {
    static int qidx = -1;
    if (qidx < 0) qidx = N - 1;
  // if ((i + 1) % 1000000 == 0)
  //   fprintf(stderr, "arr[%d] = %d, next = %d\n", qidx,arr[qidx],next_smallest);
    to_del = arr[qidx];
    to_add = next_smallest;
    arr[qidx--] = next_smallest--;
    assert(next_smallest >= 0);
  }

  int update_delete(int N) {
    static int lastN = N;
    assert(lastN > 0);
    return arr[--lastN];
  }

  void update(int N, int &to_del, int &to_add) {
    int k = dis(gen) % N;
    to_del = arr[k];
    int l = N + dis(gen) % N;
    assert(l < size());
    swap(arr[k], arr[l]);
    to_add = arr[k];
  }
};

static const char *update_workload[] = {
  "NOUP",   // 0. Read only queries.
  "LFHV",   // 1. Update 1000 tuples every 1000 queries.
  "HFLV",   // 2. Update 10 tuples every 10 queries.
  "QUEUE",  // 3. Remove the largest value, insert new smaller value than any existing value.
  "TRASH",  // 4. Insert values in the middle of the domain.
  "DELETE", // 5. Delete 1000 tuples every 1000 queries.
  "APPENDSKY", // 6. Insert 100K tuples every query. 
  "APPEND", // 7. Insert 100K tuples every query. 
  "APPEND", // 8. Insert 10 tuples every 10 query. 
};
