/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <thread>
#include <vector>

#include "random.h"
#include "ska_sort.h"
#include "time_it.h"
#include "vergesort.h"

using namespace std;

long long *allocate_elements(int size);

Random rng(140384);

#define BSIZE (1 << 12)

class Bucket {
 public:
  long long *arr;
  int n;
  Bucket *next;
  int n_flipped;
  bool original;

  Bucket(long long *a, int sz, int nf)
      : arr(a), n(sz), next(nullptr), n_flipped(nf), original(sz > 0) {
    for (int i = 1; i < sz; i++) {
      n_flipped += arr[i - 1] > arr[i];
    }
  }

  bool is_asc() const { return n_flipped == 0; }

  bool is_desc() const { return n < 2 || n_flipped == n - 1; }

  long long get_random_pivot() const { return arr[rng.nextInt(n)]; }

  long long *data() { return arr; }

  int size() const { return n; }

  int slack() const { return BSIZE - n; }

  void moveToFromIdx(Bucket *to, int fromIdx) {
    // Ensure there's something to move.
    assert(n > fromIdx);
    // Ensure the receiver has enough space.
    assert(to->n + n - fromIdx <= BSIZE);
    memmove(to->arr + to->n, arr + fromIdx, (n - fromIdx) * sizeof(long long));
    to->n += n - fromIdx;
    to->n_flipped = is_asc() ? 0 : -1;
    n = fromIdx;
  }

  void safe_edit() {
    if (original) {
      long long *src = arr;
      arr = allocate_elements(BSIZE);
      memcpy(arr, src, n * sizeof(long long));
      original = false;
    }
  }

  int partition(long long v) {
    assert(n > 0 && n <= BSIZE);
    safe_edit();
    n_flipped = -1;
    return int(
        std::partition(arr, arr + n, [&](long long x) { return (x < v); }) -
        arr);
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
};

static void fusion(long long *Lp, long long *Rp, int *hi, int *lo, int &nhi,
                   int &nlo) {
  int m = std::min(nhi, nlo);
  assert(m > 0);
  int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
  nhi -= m;
  nlo -= m;
  while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
}

class Chain {
 public:
  Bucket *head, *tail;
  int length;

  Chain() : head(nullptr), tail(nullptr), length(0) {}

  void append(Bucket *b) {
    assert(!b->next);

    if (!b->size()) {
      delete b;
      return;
    }

    if (!head) {        // the root chain is empty.
      head = tail = b;  // b is the head of the chain.
      length = 1;
      return;
    }

    assert(!tail->next);
    assert(head == tail || !head->slack());
    if (tail->slack()) {
      b->moveToFromIdx(tail, b->size() - std::min(b->size(), tail->slack()));
    }

    if (!b->size()) {
      delete b;
      return;
    }

    tail->next = b;
    tail = b;
    length++;
  }

  long long get_random_pivot() {
    long long p = head->arr[0];
    int cnt = 1;
    for (Bucket *b = head->next; b; b = b->next) {
      if (rng.nextInt(++cnt) == 0) {
        p = b->arr[0];
      }
    }
    return p;
  }

  Bucket *quick_split(Bucket *b, long long p, Chain &left_chain,
                      Chain &right_chain) {
    assert(b->n);

    if (b->is_desc()) {
      b->safe_edit();
      reverse(b->arr, b->arr + b->n);
      b->n_flipped = 0;
    }

    if (b->is_asc()) {
      if (b->arr[0] >= p) {
        right_chain.append(b);
        return nullptr;
      }
      if (b->arr[b->n - 1] < p) {
        left_chain.append(b);
        return nullptr;
      }
      return b;
    }
    return b;
  }

  // TODO: split chain should use vector and selective split up to X buffers
  pair<Chain, Chain> split_chain() {
    long long p = get_random_pivot();

    Chain left_chain;
    Chain right_chain;
    Bucket *Lb = nullptr, *Rb = nullptr;
    int hi[BSIZE], lo[BSIZE];
    int nhi = 0, nlo = 0;

    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        fusion(Lb->data(), Rb->data(), hi, lo, nhi, nlo);
        Lb->safe_edit();
        Rb->safe_edit();
        Lb->n_flipped = Rb->n_flipped = -1;
        if (!nhi) {
          left_chain.append(Lb);
          Lb = nullptr;
        }
        if (!nlo) {
          right_chain.append(Rb);
          Rb = nullptr;
        }
      } else if (!Lb) {
        if (!head) break;
        Lb = head;
        head = head->next;
        Lb->next = nullptr;
        Lb = quick_split(Lb, p, left_chain, right_chain);
      } else if (!nhi) {
        assert(Lb);
        Lb->mark_hi(p, hi, nhi);
        if (!nhi) {
          left_chain.append(Lb);
          Lb = nullptr;
        }
      } else if (!Rb) {
        if (!head) break;
        Rb = head;
        head = head->next;
        Rb->next = nullptr;
        Rb = quick_split(Rb, p, left_chain, right_chain);
      } else if (!nlo) {
        assert(Rb);
        Rb->mark_lo(p, lo, nlo);
        if (!nlo) {
          right_chain.append(Rb);
          Rb = nullptr;
        }
      } else {
        assert(0);
      }
    }

    if (Rb) {
      assert(!Lb);
      Lb = Rb;
    }
    if (Lb) {
      if (Lb->size()) {
        int i = Lb->partition(p);
        if (i == 0) {
          right_chain.append(Lb);
        } else if (i == Lb->size()) {
          left_chain.append(Lb);
        } else {
          Rb = new Bucket(allocate_elements(BSIZE), 0, -1);
          Lb->moveToFromIdx(Rb, i);
          assert(Rb->n);
          assert(Lb->n);
          left_chain.append(Lb);
          right_chain.append(Rb);
        }
      } else {
        delete Lb;
      }
    }

    // assert(left_chain);
    assert(right_chain.length);
    return make_pair(left_chain, right_chain);
  }
};

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

void ctree_sort(long long arr[], long long tmp[], int N) {
  long long *bucket_tmp = allocate_elements(BSIZE);
  Chain c;
  double it = time_it([&]() {
    for (int i = 0; i < N; i += BSIZE) {
      c.append(new Bucket(arr + i, min(BSIZE, N - i), 0));
    }
  });
  fprintf(stderr, "BSIZE = %d, chains = %d, in %.3lf; ", BSIZE, c.length, it);

  vector<Chain> stk;
  stk.push_back(c);
  while (stk.size()) {
    c = stk.back();
    stk.pop_back();

    // fprintf(stderr, "stdk sz = %u\n", stk.size());
    if (!c.head->next) {
      if (!c.head->is_asc()) {
        // sort(c.head->arr, c.head->arr + c.head->n);
        ctree_sort2(c.head->arr, bucket_tmp, c.head->n);
      }
      if (arr + N - c.head->n != c.head->arr) {
        memcpy(arr + N - c.head->n, c.head->arr, sizeof(long long) * c.head->n);
      }
      // fprintf(stderr, "block %d\n", c.head->n);
      // for (int i = c.head->n - 1; i >= 0; i--) {
      //   int j = N - c.head->n + i;
      //   if (arr[j] != j) {
      //     fprintf(stderr, "arr[%d] = %lld, %lld\n", j, arr[j],
      //     c.head->arr[i]);
      //     assert(arr[j] == j);
      //   }
      // }
      N -= c.head->n;
      delete c.head;
      continue;
    }

    auto p = c.split_chain();
    if (p.first.length) {
      assert(p.first.head->n);
      stk.push_back(p.first);
    }
    if (p.second.length) {
      assert(p.second.head->n);
      stk.push_back(p.second);
    }
  }
  assert(N == 0);
}

void merge(long long a[], int na, long long b[], int nb, long long c[]) {
  int i = 0, j = 0, k = 0;
  while (i < na && j < nb) {
    if (a[i] < b[j]) {
      c[k++] = a[i++];
    } else {
      c[k++] = b[j++];
    }
  }
  while (i < na) {
    c[k++] = a[i++];
  }
  while (j < nb) {
    c[k++] = b[j++];
  }
  assert(k == na + nb);
}

void psort(long long arr[], long long tmp[], int N) {
  int q = N / 4;
  thread t0([&]() { sort(arr, arr + q); });
  thread t1([&]() { sort(arr + q, arr + q * 2); });
  thread t2([&]() { sort(arr + q * 2, arr + q * 3); });
  thread t3([&]() { sort(arr + q * 3, arr + N); });
  t0.join();
  t1.join();
  t2.join();
  t3.join();

  thread t4([&]() { merge(arr, q, arr + q, q, tmp); });
  thread t5(
      [&]() { merge(arr + 2 * q, q, arr + 3 * q, N - 3 * q, tmp + 2 * q); });
  t4.join();
  t5.join();

  merge(tmp, 2 * q, tmp + 2 * q, N - 2 * q, arr);
}

#define MAXN (1 << 25)

long long arr[MAXN], tmp[MAXN * 4];
int n_allocated;

long long *allocate_elements(int size) {
  long long *t = tmp + n_allocated;
  n_allocated += size;
  assert(n_allocated < MAXN * 4);
  return t;
}

void test_sort(const char *name,
               std::function<void(long long[], int)> sort_algo) {
  fprintf(stderr, "%s (N = %d):\n", name, MAXN);
  Random rng(140384);
  for (int percent_sorted = 100; percent_sorted >= 50; percent_sorted -= 10) {
    for (int i = 0; i < MAXN; i++) {
      arr[i] = i;
      if (rng.nextInt(100) < 100 - percent_sorted) {
        swap(arr[i], arr[rng.nextInt(i + 1)]);
      }
    }

    n_allocated = 0;
    double sort_time = time_it([&]() { sort_algo(arr, MAXN); });

    long long chk = 0;
    for (int i = 0; i < MAXN; i++) {
      chk = chk * 13 + arr[i];
    }
    fprintf(stderr, "%3d%% sorted: %10.6lf s, chk = %20lld\n", percent_sorted,
            sort_time, chk);
  }
  fprintf(stderr, "\n");
}

int main() {
  test_sort("ctree_sort",
            [](long long arr[], int n) { ctree_sort(arr, tmp, n); });
  test_sort("std::sort", [](long long arr[], int n) { sort(arr, arr + n); });
  // test_sort("std::psort", [](long long arr[], int n) { psort(arr, tmp, n);
  // });
  test_sort("vergesort",
            [](long long arr[], int n) { vergesort::vergesort(arr, arr + n);
            });
  // test_sort("ska_sort", [](long long arr[], int n) { ska_sort(arr, arr + n);
  // });
}
