/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <set>
#include <vector>

#include "parallel_sort.h"
#include "random.h"
#include "ska_sort.h"
#include "time_it.h"
#include "vergesort.h"

using namespace std;

static long long *allocate_elements(int size);
static void deallocate(long long *ptr);

Random rng(140384);

#ifdef DBG
#define MAXN (1 << 20)
#define BSIZE (1 << 10)
#define assert_dbg(x) assert(x)
#else
#define MAXN (1 << 26)
#define BSIZE (1 << 12)
#define assert_dbg(x) ((void)0)
#endif

static int nswap = 0;

class Bucket {
 public:
  long long *arr;
  int n;
  int sorted;

  Bucket(long long *a, int sz) : arr(a), n(sz), sorted(0) {
    assert(n <= BSIZE);
    if (!n) return;

    int flipped = 0;
    for (int i = 1; i < n && !flipped; i++) {
      flipped = arr[i - 1] > arr[i];
    }
    sorted = !flipped;
  }

  long long get_random_pivot() const { return arr[rng.nextInt(n)]; }

  long long *data() { return arr; }

  int size() const { return n; }

  bool empty() const { return n == 0; }

  int slack() const { return BSIZE - n; }

  long long first() const { return arr[0]; }

  long long last() const { return arr[n - 1]; }

  void moveToFromIdx(Bucket &to, int fromIdx) {
    assert(arr != to.arr);
    assert(n > fromIdx);
    assert(to.n + n - fromIdx <= BSIZE);

    memcpy(to.arr + to.n, arr + fromIdx, (n - fromIdx) * sizeof(long long));
    to.n += n - fromIdx;
    n = fromIdx;
  }

  void copyTo(long long *target) { memcpy(target, arr, sizeof(long long) * n); }

  int sample_nhi(long long p) const {
    int cnt = 0;
    for (int i = 0; i < n && i < 10; i++) {
      cnt += (arr[rng.nextInt(n)] >= p);
    }
    return cnt;
  }

  void mark_hi(long long P, int *hi, int &nhi) const {
    for (int i = 0; i < n; i++) {
      hi[nhi] = i;
      nhi += (arr[i] >= P);
    }
  }

  void mark_lo(long long P, int *lo, int &nlo) const {
    for (int i = 0; i < n; i++) {
      lo[nlo] = i;
      nlo += (arr[i] < P);
    }
  }

  // Returns a new bucket containing elements >= p from this bucket preserving
  // its order. Elements >= p is stably removed from this bucket.
  Bucket split(long long p) {
    assert(n > 0 && n <= BSIZE);
    Bucket b(allocate_elements(BSIZE), 0);
    long long *a = arr;
    int i = n;
    n = 0;
    for (; i--; a++) {
      int is_less = *a < p;
      arr[n] = *a;
      n += is_less;
      b.arr[b.n] = *a;
      b.n += !is_less;
    }
    // assert(!b.empty());
    return b;
  }

  void sort(long long *tmp) { vergesort::vergesort(arr, arr + n); }

  bool is_less_than(long long p) {
    for (int i = 0; i < n; i++) {
      if (arr[i] >= p) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) {
    for (int i = 0; i < n; i++) {
      if (arr[i] < p) {
        return false;
      }
    }
    return true;
  }
};

#define block_size 64

static void swap_offsets(long long *first, long long *last, unsigned char *hi,
                         unsigned char *lo, int num) {
  if (num > 0) {
    long long *l = first + hi[0];
    long long *r = last - lo[0];
    long long tmp(*l);
    *l = *r;
    for (int i = 1; i < num; ++i) {
      l = first + hi[i];
      *r = (*l);
      r = last - lo[i];
      *l = (*r);
    }
    *r = (tmp);
  }
}

static pair<int, int> partition_branchless(long long *L, int nL, long long *R,
                                           int nR, long long pivot) {
  long long *first = L;
  long long *last = R + nR;
  auto comp = std::less<long long>();

  // TODO: Unroll.
  while (first < L + nL && comp(*first, pivot)) first++;
  while (R < last && !comp(*(last - 1), pivot)) last--;

  // The branchless partitioning is derived from "BlockQuicksort: How Branch
  // Mispredictions don’t affect Quicksort" by Stefan Edelkamp and Armin
  // Weiss.
  unsigned char hi[block_size];
  unsigned char lo[block_size];
  int nhi = 0, nlo = 0, ihi = 0, ilo = 0;

  while (first + block_size <= L + nL && last - block_size >= R) {
    if (nhi == 0) {
      ihi = 0;
      long long *it = first;
      for (unsigned char i = 0; i < block_size;) {
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
      }
    }
    if (nlo == 0) {
      ilo = 0;
      long long *it = last;
      for (unsigned char i = 0; i < block_size;) {
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
      }
    }

    // Swap elements and update block sizes and first/last boundaries.
    int num = std::min(nhi, nlo);
    swap_offsets(first, last, hi + ihi, lo + ilo, num);
    nhi -= num;
    nlo -= num;
    ihi += num;
    ilo += num;
    if (nhi == 0) first += block_size;
    if (nlo == 0) last -= block_size;
  }

  int l_size = min(block_size, int(L + nL - first));
  if (first < L + nL && !nhi) {
    ihi = 0;
    long long *it = first;
    for (unsigned char i = 0; i < l_size;) {
      hi[nhi] = i++;
      nhi += !comp(*it++, pivot);
    }
  }

  int r_size = min(block_size, int(last - R));
  if (last > R && !nlo) {
    ilo = 0;
    long long *it = last;
    for (unsigned char i = 0; i < r_size;) {
      lo[nlo] = ++i;
      nlo += comp(*--it, pivot);
    }
  }

  int num = std::min(nhi, nlo);
  swap_offsets(first, last, hi + ihi, lo + ilo, num);
  nhi -= num;
  nlo -= num;
  ihi += num;
  ilo += num;
  if (nhi == 0) first += l_size;
  if (nlo == 0) last -= r_size;
  assert(first <= L + nL);
  assert(last >= R);

  return make_pair(first - L, R + nR - last);
}

class Chain {
 public:
  vector<Bucket> buckets;

  void append(Bucket b) {
    if (b.empty()) {
      deallocate(b.arr);
      return;
    }

    if (!buckets.empty() && buckets.back().slack()) {
      Bucket &tail = buckets.back();
      if (b.slack() < tail.slack()) {
        swap(b, tail);
      }
      if (tail.slack()) {
        b.moveToFromIdx(tail, b.size() - std::min(b.size(), tail.slack()));
      }
      if (b.empty()) {
        deallocate(b.arr);
        return;
      }
    }

    buckets.push_back(b);
  }

  bool is_less_than(long long p) {
    for (Bucket b : buckets) {
      if (!b.is_less_than(p)) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) {
    for (Bucket b : buckets) {
      if (!b.is_at_least(p)) {
        return false;
      }
    }
    return true;
  }

  bool distinct() {
    set<long long> s;
    for (Bucket b : buckets) {
      auto it = s.insert(b.arr[0]);
      if (!it.second) return false;
    }
    return true;
  }

  bool quick_split(Bucket &b, long long p, Chain &left_chain,
                   Chain &right_chain) {
    if (b.sorted) {
      if (b.first() >= p) {
        right_chain.append(b);
        return true;
      }
      if (b.last() < p) {
        left_chain.append(b);
        return true;
      }
    }
    return false;
  }

  pair<Chain, Chain> split_chain() {
    assert(!buckets.empty());
    long long p = buckets[rng.nextInt(int(buckets.size()))].get_random_pivot();

    Chain left_chain, right_chain;

    // for (Bucket b : buckets) {
    //   if (left_chain.buckets.size() > right_chain.buckets.size()) {
    //     right_chain.append(b);
    //   } else {
    //     left_chain.append(b);
    //   }
    // }
    // return make_pair(left_chain, right_chain);

    int i = 0, j = int(buckets.size()) - 1;
    Bucket *L = nullptr;
    Bucket *R = nullptr;
    int ihi = 0, ilo = 0;
    while (i <= j) {
      // fprintf(stderr, "%d %d\n", i, j);
      if (L == nullptr) {
        L = &buckets[i++];
        ihi = 0;
        if (quick_split(*L, p, left_chain, right_chain)) {
          L = nullptr;
        }
      } else if (R == nullptr) {
        R = &buckets[j--];
        ilo = 0;
        if (quick_split(*R, p, left_chain, right_chain)) {
          R = nullptr;
        }
      } else {
        L->sorted = R->sorted = false;
        auto m = partition_branchless(L->arr + ihi, L->n - ihi, R->arr,
                                      R->n - ilo, p);
        // fprintf(stderr, "m %d %d\n", m.first, m.second);
        ihi += m.first;
        assert(ihi <= L->n);
        ilo += m.second;
        assert(ilo <= R->n);
        if (ihi == L->n) {
          assert_dbg(left_chain.is_less_than(p));
          left_chain.append(*L);
          assert_dbg(left_chain.is_less_than(p));
          L = nullptr;
        }
        if (ilo == R->n) {
          assert_dbg(right_chain.is_at_least(p));
          right_chain.append(*R);
          assert_dbg(right_chain.is_at_least(p));
          R = nullptr;
        }
      }
    }
    if (L) {
      L->sorted = false;
      right_chain.append(L->split(p));
      assert_dbg(right_chain.is_at_least(p));
      left_chain.append(*L);
      assert_dbg(left_chain.is_less_than(p));
    }
    if (R) {
      R->sorted = false;
      right_chain.append(R->split(p));
      assert_dbg(right_chain.is_at_least(p));
      left_chain.append(*R);
      assert_dbg(left_chain.is_less_than(p));
    }
    reverse(right_chain.buckets.begin(), right_chain.buckets.end());
    return make_pair(left_chain, right_chain);
  }
};

void ctree_sort(long long arr[], long long tmp[], int N) {
  Chain c;
  double it = time_it([&]() {
    for (int i = 0; i < N; i += BSIZE) {
      int sz = min(BSIZE, N - i);
      long long *a = allocate_elements(BSIZE);
      memcpy(a, arr + i, sz * sizeof(long long));
      c.append(Bucket(a, sz));
    }
  });
  nswap = 0;
  size_t csize = c.buckets.size();

  assert_dbg(c.distinct());

  long long *bucket_tmp = allocate_elements(BSIZE);
  vector<Chain> stk;
  stk.push_back(c);
  while (stk.size()) {
    c = stk.back();
    stk.pop_back();

    // fprintf(stderr, "Stack = %lu\n", stk.size());
    if (c.buckets.size() == 1) {
      Bucket b = c.buckets[0];
      b.sort(bucket_tmp);
      b.copyTo(arr + N - b.size());

      // fprintf(stderr, "block %d\n", b.n);
      // for (int i = b.n - 1; i >= 0; i--) {
      //   int j = N - b.n + i;
      //   if (arr[j] != j) {
      //     fprintf(stderr, "arr[%d] = %lld, %lld\n", j, arr[j], b.arr[i]);
      //     assert(arr[j] == j);
      //   }
      // }

      N -= b.size();
      continue;
    }

    auto p = c.split_chain();
    if (!p.first.buckets.empty()) {
      assert_dbg(p.first.distinct());
      assert(!p.first.buckets[0].empty());
      stk.push_back(p.first);
    }
    if (!p.second.buckets.empty()) {
      assert_dbg(p.second.distinct());
      assert(!p.second.buckets[0].empty());
      stk.push_back(p.second);
    }
  }
  assert(N == 0);

  fprintf(stderr, "BSIZE = %d, chains = %lu, in %.3lf, nswap = %9d; ", BSIZE,
          csize, it, nswap);
}

long long arr[MAXN], tmp[MAXN * 10];
vector<long long *> free_ptrs;
int n_allocated;

static long long *allocate_elements(int size) {
  if (free_ptrs.empty()) {
    long long *t = tmp + n_allocated;
    n_allocated += size;
    assert(n_allocated < MAXN * 10);
    return t;
  }
  long long *ret = free_ptrs.back();
  free_ptrs.pop_back();
  return ret;
}

static void deallocate(long long *ptr) { free_ptrs.push_back(ptr); }

void test_sort(const char *name,
               std::function<void(long long[], int)> sort_algo) {
  fprintf(stderr, "%s (N = %d):\n", name, MAXN);
  Random random(140384);
  for (int percent_sorted = 100000, m = 1;;) {
    percent_sorted = max(0, percent_sorted);
    int misplaced = 0;
    for (int i = 0; i < MAXN; i++) {
      arr[i] = i;
      if (random.nextInt(100000) < 100000 - percent_sorted) {
        swap(arr[i], arr[random.nextInt(i + 1)]);
        misplaced++;
      }
    }

    n_allocated = 0;
    free_ptrs.clear();
    double sort_time = time_it([&]() { sort_algo(arr, MAXN); });

    long long chk = 0;
    for (int i = 0; i < MAXN; i++) {
      chk = chk * 13 + arr[i];
    }
    fprintf(stderr, "%7.3lf%% sorted: %10.6lf s, nrswap = %8d, chk = %20lld\n",
            percent_sorted / 1000.0, sort_time, misplaced, chk);

    if (percent_sorted <= 0) break;
    percent_sorted -= m;
    m *= 2;
  }
  fprintf(stderr, "\n");
}

int main() {
  for (int i = 0; i < MAXN; i++) {
    tmp[i] = rng.nextInt();
  }
  test_sort("std::sort", [](long long arr[], int n) { sort(arr, arr + n); });
  test_sort("ctree_sort",
            [](long long arr[], int n) { ctree_sort(arr, tmp, n); });
  test_sort("std::parallel_sort",
            [](long long arr[], int n) { parallel_sort(arr, n); });
  test_sort("vergesort",
            [](long long arr[], int n) { vergesort::vergesort(arr, arr + n); });
  test_sort("ska_sort", [](long long arr[], int n) { ska_sort(arr, arr + n); });
}
