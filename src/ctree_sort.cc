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

long long *allocate_elements(int size);

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

static void ctree_sort2(long long *arr, long long *tmp, int n) {
  if (n < 30) {
    sort(arr, arr + n);
    return;
  }

  long long P = arr[rng.nextInt(n)];
  int nlo = 0, nhi = n - 1;
  for (int i = 0; i < n; i++) {
    tmp[nlo] = tmp[nhi] = arr[i];
    int j = arr[i] < P;
    nhi -= 1 - j;
    nlo += j;
  }

  ctree_sort2(tmp, arr, nlo);
  ctree_sort2(tmp + nlo, arr + nlo, (n - nlo));

  memcpy(arr, tmp, sizeof(long long) * nlo);
  memcpy(arr + nlo, tmp + nlo, sizeof(long long) * (n - nlo));
}

static int nswap = 0;

class Bucket {
 public:
  long long *arr;
  int n;
  int n_flipped;

  Bucket(int sz) : arr(allocate_elements(BSIZE)), n(sz) { assert(n <= BSIZE); }

  Bucket(long long *a, int sz)
      : arr(allocate_elements(BSIZE)), n(sz), n_flipped(-1) {
    assert(n <= BSIZE);
    memcpy(arr, a, n * sizeof(long long));
  }

  bool is_asc() const { return n_flipped == 0; }

  bool is_desc() const { return n == 0 || n_flipped == n - 1; }

  bool is_sorted() const { return !is_asc() && !is_desc(); }

  long long get_random_pivot() const { return arr[rng.nextInt(n)]; }

  long long *data() { return arr; }

  int size() const { return n; }

  bool empty() const { return n == 0; }

  int slack() const { return BSIZE - n; }

  long long first() const { return arr[0]; }

  long long last() const { return arr[n - 1]; }

  int count_n_flipped() {
    n_flipped = 0;
    for (int i = 1; i < n; i++) {
      n_flipped += arr[i - 1] > arr[i];
    }
    return n_flipped;
  }

  void moveToFromIdx(Bucket &to, int fromIdx) {
    assert(arr != to.arr);
    assert_dbg(check_n_flipped());
    assert_dbg(to.check_n_flipped());
    // Ensure there's something to move.
    assert(n > fromIdx);
    // Ensure the receiver has enough space.
    assert(to.n + n - fromIdx <= BSIZE);

    if (to.is_asc()) {
      // TODO: maintain sortedness if possible.
      // fprintf(stderr, "here\n");
    }

    memcpy(to.arr + to.n, arr + fromIdx, (n - fromIdx) * sizeof(long long));
    to.n += n - fromIdx;

    if (is_asc() && to.empty()) {
      to.n_flipped = 0;
    } else {
      // TOOD: may be expensive.
      to.count_n_flipped();
    }

    n = fromIdx;
    if (n_flipped > 0) {
      // TODO: may be expensive.
      count_n_flipped();
    }

    assert_dbg(check_n_flipped());
    assert_dbg(to.check_n_flipped());
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
    // fprintf(stderr, "s");
    assert(n > 0 && n <= BSIZE);
    if (is_desc()) {
      reverse(arr, arr + n);
      n_flipped = 0;
    }
    if (is_asc()) {
      assert_dbg(check_n_flipped());
      int old_n = n;
      n = int(lower_bound(arr, arr + n, p) - arr);
      return Bucket(arr + n, old_n - n);
    }
    Bucket b(0);
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
    count_n_flipped();
    assert(!b.empty());
    b.count_n_flipped();
    assert_dbg(check_n_flipped());
    assert_dbg(b.check_n_flipped());
    return b;
  }

  void sort(long long *tmp) {
    if (is_desc()) {
      reverse(arr, arr + n);
    } else if (!is_asc()) {
      if (n_flipped * 4 > n && n_flipped * 4 < n * 3) {
        ctree_sort2(arr, tmp, n);
      } else {
        vergesort::vergesort(arr, arr + n);
        // std::sort(arr, arr + n);
      }
    }
    // vergesort::vergesort(arr, arr + n);

    n_flipped = 0;
    assert_dbg(count_n_flipped() == 0);
  }

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

  bool check_n_flipped() {
    if (is_asc()) {
      return count_n_flipped() == 0;
    }
    if (is_desc()) {
      return n == 0 || count_n_flipped() == n - 1;
    }
    return true;
  }
};

class Chain {
 public:
  vector<Bucket> buckets;

  void append(Bucket b) {
    if (b.empty()) return;

    if (!buckets.empty() && buckets.back().slack()) {
      Bucket &tail = buckets.back();
      if (b.slack() < tail.slack()) {
        swap(b, tail);
      }
      if (tail.slack()) {
        b.moveToFromIdx(tail, b.size() - std::min(b.size(), tail.slack()));
      }
      if (b.empty()) return;
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

  bool check_n_flipped() {
    for (Bucket b : buckets) {
      if (!b.check_n_flipped()) {
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

  pair<Chain, Chain> split_chain() {
    assert(!buckets.empty());
    long long p = buckets[rng.nextInt(int(buckets.size()))].arr[0];

    Chain left_chain;
    Chain right_chain;

    // for (Bucket b : buckets) {
    //   if (left_chain.buckets.size() > right_chain.buckets.size()) {
    //     right_chain.append(b);
    //   } else {
    //     left_chain.append(b);
    //   }
    // }
    // return make_pair(left_chain, right_chain);

    vector<pair<int, Bucket>> pending;
    for (Bucket &b : buckets) {
      assert_dbg(b.check_n_flipped());

      if (b.is_desc()) {
        reverse(b.arr, b.arr + b.n);
        b.n_flipped = 0;
      }

      if (b.is_asc()) {
        if (b.first() >= p) {
          right_chain.append(b);
          continue;
        }
        if (b.last() < p) {
          left_chain.append(b);
          continue;
        }
      }

      pending.push_back(make_pair(b.sample_nhi(p), b));
    }

    assert_dbg(left_chain.is_less_than(p));
    assert_dbg(left_chain.check_n_flipped());

    assert_dbg(right_chain.is_at_least(p));
    assert_dbg(right_chain.check_n_flipped());

    sort(pending.begin(), pending.end(),
         [](const auto &a, const auto &b) { return a.first < b.first; });

    // if (int(pending.size()) > 0) {
    //   fprintf(stderr, "npending = %lu, nswap = %d\n", pending.size(), nswap);
    //   for (auto p : pending) {
    //     fprintf(stderr, "%d ", p.first);
    //   }
    //   fprintf(stderr, "\n");
    // }

    int lo[BSIZE], hi[BSIZE];
    int nhi = 0, nlo = 0;
    int i = 0, j = int(pending.size());
    Bucket *L, *R = nullptr;

    while (true) {
      if (nhi == 0) {
        if (i < j) {
          L = &pending[i++].second;
          L->mark_hi(p, hi, nhi);

          if (nhi == 0) {
            left_chain.append(*L);
            continue;
          }

          if (nhi == L->size()) {
            right_chain.append(*L);
            nhi = 0;
            continue;
          }

          // TODO: separate out minor outliers, keep sortedness.
        } else {
          break;
        }
      }

      if (nlo == 0) {
        if (i < j) {
          R = &pending[--j].second;
          R->mark_lo(p, lo, nlo);

          if (nlo == 0) {
            right_chain.append(*R);
            continue;
          }

          if (nlo == R->size()) {
            left_chain.append(*R);
            nlo = 0;
            continue;
          }

        } else {
          break;
        }
      }

      int m = std::min(nhi, nlo);
      // fprintf(stderr, "m = %d\n", m);
      assert(m > 0);

      nswap += m;

      long long *Lp = L->arr;
      long long *Rp = R->arr;
      int *hip = hi + nhi - 1;
      int *lop = lo + nlo - 1;
      L->n_flipped = -1;
      R->n_flipped = -1;

      for (int k = m; k--;) std::swap(Lp[*(hip--)], Rp[*(lop--)]);

      nhi -= m;
      nlo -= m;

      if (nhi == 0) {
        assert_dbg(left_chain.is_less_than(p));
        left_chain.append(*L);
        assert_dbg(left_chain.is_less_than(p));
      }

      if (nlo == 0) {
        assert_dbg(right_chain.is_at_least(p));
        right_chain.append(*R);
        assert_dbg(right_chain.is_at_least(p));
      }
    }

    assert_dbg(left_chain.is_less_than(p));
    assert_dbg(left_chain.check_n_flipped());

    assert_dbg(right_chain.is_at_least(p));
    assert_dbg(right_chain.check_n_flipped());

    if (nhi) {
      assert_dbg(L->check_n_flipped());
      int n = L->size();
      Bucket b = L->split(p);
      assert(b.n + L->size() == n);
      assert(b.n == nhi);
      assert_dbg(L->check_n_flipped());
      assert_dbg(b.check_n_flipped());
      if (!L->empty()) {
        left_chain.append(*L);
      }
      if (!b.empty()) {
        right_chain.append(b);
      }
    }

    if (nlo) {
      assert_dbg(R->check_n_flipped());
      int n = R->size();
      Bucket b = R->split(p);
      assert(b.n + R->size() == n);
      assert(R->size() == nlo);
      assert_dbg(R->check_n_flipped());
      assert_dbg(b.check_n_flipped());
      assert(!R->empty());
      assert(!b.empty());
      if (!R->empty()) {
        left_chain.append(*R);
      }
      if (!b.empty()) {
        right_chain.append(b);
      }
    }

    assert_dbg(left_chain.is_less_than(p));
    assert_dbg(left_chain.check_n_flipped());

    assert_dbg(right_chain.is_at_least(p));
    assert_dbg(right_chain.check_n_flipped());

    // fprintf(stderr, "split %lu %lu\n", left_chain.buckets.size(),
    //         right_chain.buckets.size());

    return make_pair(left_chain, right_chain);
  }
};

void ctree_sort(long long arr[], long long tmp[], int N) {
  Chain c;
  double it = time_it([&]() {
    for (int i = 0; i < N; i += BSIZE) {
      c.append(Bucket(arr + i, min(BSIZE, N - i)));
      c.buckets.back().count_n_flipped();
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

long long arr[MAXN], tmp[MAXN * 20];
int n_allocated;

long long *allocate_elements(int size) {
  long long *t = tmp + n_allocated;
  n_allocated += size;
  assert(n_allocated < MAXN * 20);
  return t;
}

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
