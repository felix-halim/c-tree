#include <algorithm>
#include <vector>

#include "../random.h"
#include "vergesort.h"

#ifdef DBG
#define assert_dbg(x) assert(x)
#else
#define assert_dbg(x) ((void)0)
#endif

static int nswap;
static Random rng(140384);

static const int tmp_max = 1 << 27;
static long long *tmp;
static std::vector<long long *> free_ptrs;
static int tmp_size;

static long long *allocate_elements(int size) {
  if (free_ptrs.empty()) {
    long long *t = tmp + tmp_size;
    tmp_size += size;
    assert(tmp_size < tmp_max);
    return t;
  }
  long long *ret = free_ptrs.back();
  free_ptrs.pop_back();
  return ret;
}

static void deallocate_elements(long long *ptr) { free_ptrs.push_back(ptr); }

static const int block_size = 64;

static void swap_offsets(long long *first, long long *last, unsigned char *hi,
                         unsigned char *lo, int num) {
  if (num > 0) {
    long long *l = first + hi[0];
    long long *r = last - lo[0];
    long long t(*l);
    *l = *r;
    for (int i = 1; i < num; ++i) {
      l = first + hi[i];
      *r = (*l);
      r = last - lo[i];
      *l = (*r);
    }
    *r = (t);
    nswap += num;
  }
}

static std::pair<int, int> partition_branchless(long long *L, int nL,
                                                long long *R, int nR,
                                                long long pivot) {
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

  int l_size = std::min(block_size, int(L + nL - first));
  if (first < L + nL && !nhi) {
    ihi = 0;
    long long *it = first;
    for (unsigned char i = 0; i < l_size;) {
      hi[nhi] = i++;
      nhi += !comp(*it++, pivot);
    }
  }

  int r_size = std::min(block_size, int(last - R));
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

  return std::make_pair(first - L, R + nR - last);
}

template <int BCAP>
class Bucket {
  long long *arr;
  int n;
  int sorted;  // 0: no, 1: yes, 2: unknown

 public:
  Bucket(long long *a, int sz) : arr(a), n(sz), sorted(2) { assert(n <= BCAP); }

  int capacity() const { return BCAP; }
  int size() const { return n; }
  bool empty() const { return n == 0; }
  int slack() const { return BCAP - n; }
  long long first() const { return arr[0]; }
  long long last() const { return arr[n - 1]; }
  bool at_least_half_empty() const { return n * 2 <= BCAP; }
  long long get_random_pivot() const { return arr[rng.nextInt(n)]; }

  void deallocate() const { deallocate_elements(arr); }

  bool is_sorted() {
    if (sorted == 2) {
      int flipped = 0;
      for (int i = 1; i < n && !flipped; i++) {
        flipped = arr[i - 1] > arr[i];
      }
      sorted = !flipped;
    }
    return sorted;
  }

  void sort() {
    if (sorted == 1) return;
    // std::sort(arr, arr + n);
    vergesort::vergesort(arr, arr + n);
    sorted = 1;
  }

  void copyTo(long long *target) const {
    memcpy(target, arr, sizeof(long long) * n);
  }

  void partition(int &ihi, Bucket *R, int &ilo, long long p) {
    int prev_nswap = nswap;
    auto m = partition_branchless(arr + ihi, n - ihi, R->arr, R->n - ilo, p);
    ihi += m.first;
    assert(ihi <= n);
    ilo += m.second;
    assert(ilo <= R->n);
    if (prev_nswap != nswap) {
      sorted = R->sorted = 2;
    }
  }

  // Returns a new bucket containing elements >= p from this bucket preserving
  // its order. Elements >= p is stably removed from this bucket.
  Bucket split(long long p) {
    assert(n > 0 && n <= BCAP);
    Bucket b(allocate_elements(BCAP), 0);
    long long *a = arr;
    long long *x = arr;
    long long *y = b.arr;
    for (int i = n; i--;) {
      int is_less = *a < p;
      *x = *y = *a++;
      x += is_less;
      y += !is_less;
    }
    n = int(x - arr);
    b.n = int(y - b.arr);
    return b;
  }

  bool is_less_than(long long p) const {
    for (int i = 0; i < n; i++) {
      if (arr[i] >= p) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) const {
    for (int i = 0; i < n; i++) {
      if (arr[i] < p) {
        return false;
      }
    }
    return true;
  }
};

template <class T>
class Chain {
  std::vector<T> buckets;

 public:
  bool empty() const { return buckets.empty(); }

  int size() const { return int(buckets.size()); }

  long long get_random_pivot() const {
    return buckets[rng.nextInt(size())].get_random_pivot();
  }

  void trim_last() {
    if (!buckets.empty() && !buckets.back().size()) {
      buckets.back().deallocate();
      buckets.pop_back();
    }
  }

  T &ref(int i) { return buckets[i]; }

  T &get_bucket_to_append(std::function<T()> new_bucket) {
    if (buckets.empty() || buckets.back().at_least_half_empty()) {
      buckets.push_back(new_bucket());
    }
    return buckets.back();
  }

  void append(T b) {
    if (b.empty()) {
      b.deallocate();
      return;
    }
    trim_last();
    buckets.push_back(b);
  }

  bool is_less_than(long long p) const {
    for (const auto &b : buckets) {
      if (!b.is_less_than(p)) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) const {
    for (const auto &b : buckets) {
      if (!b.is_at_least(p)) {
        return false;
      }
    }
    return true;
  }

  bool distinct() const {
    std::vector<long long> s;
    for (const auto &b : buckets) {
      s.push_back(b.arr[0]);
    }
    std::sort(s.begin(), s.end());
    for (int i = 1; i < int(s.size()); i++) {
      if (s[i - 1] == s[i]) return false;
    }
    return true;
  }

  bool quick_split(T &b, long long p, Chain<T> &left_chain,
                   Chain<T> &right_chain) const {
    if (b.is_sorted()) {
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

  std::pair<Chain<T>, Chain<T>> split_chain(std::function<T()> new_bucket) {
    assert(!empty());
    long long p = get_random_pivot();

    Chain<T> left_chain, right_chain;
    for (int nth = 0, k = 0; nth < size(); nth++) {
      T &L = left_chain.get_bucket_to_append(new_bucket);
      T &R = right_chain.get_bucket_to_append(new_bucket);
      T &b = ref(nth);

      long long *x = L.arr + L.n;
      long long *y = R.arr + R.n;
      long long *z = b.arr + k;
      int n =
          std::min(b.n - k, std::min(L.capacity() - L.n, R.capacity() - R.n));
      k += n;
      for (int j = n; j--; z++) {
        int is_less = *z < p;
        *x = *z;
        x += is_less;
        *y = *z;
        y += !is_less;
      }
      L.n = int(x - L.arr);
      R.n = int(y - R.arr);

      if (k < b.n) {
        nth--;
      } else {
        b.deallocate();
        k = 0;
      }
    }

    left_chain.trim_last();
    right_chain.trim_last();

    // fprintf(stderr, "%lu -> %lu %lu\n", buckets.size(),
    // left_chain.buckets.size(), right_chain.buckets.size());
    return std::make_pair(left_chain, right_chain);
  }

  std::pair<Chain<T>, Chain<T>> split_chain_noop(std::function<T()> ignored) {
    Chain<T> left_chain, right_chain;
    for (int i = 0; i < size(); i++) {
      T &b = ref(i);
      if (left_chain.size() > right_chain.size()) {
        right_chain.append(b);
      } else {
        left_chain.append(b);
      }
    }
    return std::make_pair(left_chain, right_chain);
  }

  std::pair<Chain<T>, Chain<T>> split_chain_nb(std::function<T()> new_bucket) {
    long long p = get_random_pivot();
    Chain<T> left_chain, right_chain;
    int i = 0, j = size() - 1;
    T *L = nullptr;
    T *R = nullptr;
    int ihi = 0, ilo = 0;
    while (i <= j) {
      // fprintf(stderr, "%d %d\n", i, j);
      if (L == nullptr) {
        L = &ref(i++);
        ihi = 0;
        if (quick_split(*L, p, left_chain, right_chain)) {
          L = nullptr;
        }
      } else if (R == nullptr) {
        R = &ref(j--);
        ilo = 0;
        if (quick_split(*R, p, left_chain, right_chain)) {
          R = nullptr;
        }
      } else {
        L->partition(ihi, R, ilo, p);
        if (ihi == L->size()) {
          assert_dbg(left_chain.is_less_than(p));
          left_chain.append(*L);
          assert_dbg(left_chain.is_less_than(p));
          L = nullptr;
        }
        if (ilo == R->size()) {
          assert_dbg(right_chain.is_at_least(p));
          right_chain.append(*R);
          assert_dbg(right_chain.is_at_least(p));
          R = nullptr;
        }
      }
    }
    if (L) {
      right_chain.append(L->split(p));
      assert_dbg(right_chain.is_at_least(p));
      left_chain.append(*L);
      assert_dbg(left_chain.is_less_than(p));
    }
    if (R) {
      right_chain.append(R->split(p));
      assert_dbg(right_chain.is_at_least(p));
      left_chain.append(*R);
      assert_dbg(left_chain.is_less_than(p));
    }
    std::reverse(right_chain.buckets.begin(), right_chain.buckets.end());
    return std::make_pair(left_chain, right_chain);
  }
};

template <int BCAP>
static void ctreesort(long long arr[], int N) {
  Chain<Bucket<BCAP>> c;

  double it = time_it([&]() {
    if (!tmp) {
      tmp = new long long[tmp_max];
    }
    tmp_size = 0;
    free_ptrs.clear();

    for (int i = 0; i < N; i += BCAP) {
      int sz = std::min(BCAP, N - i);
      long long *a = allocate_elements(BCAP);
      memcpy(a, arr + i, sz * sizeof(long long));
      c.append(Bucket<BCAP>(a, sz));
    }
  });

  nswap = 0;
  size_t csize = c.size();

  assert_dbg(c.distinct());

  std::vector<Chain<Bucket<BCAP>>> stk;
  stk.push_back(c);
  while (stk.size()) {
    c = stk.back();
    stk.pop_back();

    // fprintf(stderr, "Stack = %lu\n", stk.size());
    if (c.size() <= 3) {
      int R = N;
      for (int i = c.size() - 1; i >= 0; i--) {
        auto &b = c.ref(i);
        N -= b.size();
        b.copyTo(arr + N);
        b.deallocate();
      }
      vergesort::vergesort(arr + N, arr + R);
      continue;
    }

    auto p = c.split_chain_nb(
        []() { return Bucket<BCAP>(allocate_elements(BCAP), 0); });
    if (p.first.size()) {
      assert_dbg(p.first.distinct());
      assert(!p.first.ref(0).empty());
      stk.push_back(p.first);
    }
    if (p.second.size()) {
      assert_dbg(p.second.distinct());
      assert(!p.second.ref(0).empty());
      stk.push_back(p.second);
    }
  }
  assert(N == 0);

  fprintf(stderr, "BCAP = %d, chains = %lu, in %.3lf, nswap = %9d; ", BCAP,
          csize, it, nswap);
}
