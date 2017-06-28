/*
To run, go to the root folder and execute:

make ctree_sort
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <vector>

#include "random.h"
#include "ska_sort.h"
#include "time_it.h"

using namespace std;

Random rng(140384);

#define BSIZE 4096

class Bucket {
 public:
  long long arr[BSIZE];
  int n;
  Bucket *next;

  Bucket() : n(0), next(nullptr) {}

  long long get_random_pivot() { return arr[rng.nextInt(n)]; }

  long long *data() { return arr; }

  int size() { return n; }

  int slack() const { return BSIZE - n; }

  void moveToFromIdx(Bucket *to, int fromIdx) {
    // Ensure there's something to move.
    assert(n > fromIdx);
    // Ensure the receiver has enough space.
    assert(to->n + n - fromIdx <= BSIZE);
    memmove(to->arr + to->n, arr + fromIdx, (n - fromIdx) * sizeof(long long));
    to->n += n - fromIdx;
    n = fromIdx;
  }

  int partition(long long const &v) {
    assert(n > 0 && n <= BSIZE);
    return int(
        std::partition(arr, arr + n, [&](long long x) { return (x < v); }) -
        arr);
  }

  void mark_hi(long long const &P, int *hi, int &nhi) {
    for (int i = 0; i < n; i++) {
      hi[nhi] = i;
      nhi += (arr[i] >= P);
    }
  }

  void mark_lo(long long const &P, int *lo, int &nlo) {
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

  pair<Chain, Chain> split_chain() {
    long long p = head->get_random_pivot();

    Chain left_chain;
    Chain right_chain;
    Bucket *Lb = nullptr, *Rb = nullptr;
    int hi[BSIZE], lo[BSIZE];
    int nhi = 0, nlo = 0;

    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        fusion(Lb->data(), Rb->data(), hi, lo, nhi, nlo);
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
          Rb = new Bucket();
          Lb->moveToFromIdx(Rb, i);
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

void ctree_sort(long long arr[], int N) {
  long long tmp[BSIZE];
  Chain c;
  for (int i = 0; i < N;) {
    Bucket *b = new Bucket();
    b->n = min(BSIZE, N - i);
    memcpy(b->arr, arr + i, sizeof(long long) * b->n);
    i += b->n;
    c.append(b);
  }
  fprintf(stderr, "chains = %d\n", c.length);

  vector<Chain> stk;
  stk.push_back(c);
  while (stk.size()) {
    c = stk.back();
    stk.pop_back();
    if (!c.head->next) {
      ctree_sort2(c.head->arr, tmp, c.head->n);
      // sort(c.head->arr, c.head->arr + c.head->n);
      memcpy(arr + N - c.head->n, c.head->arr, sizeof(long long) * c.head->n);
      N -= c.head->n;
      delete c.head;
      continue;
    }

    auto p = c.split_chain();
    if (p.first.length) {
      stk.push_back(p.first);
    }
    if (p.second.length) {
      stk.push_back(p.second);
    }
  }
  assert(N == 0);
}

#define MAXN 100000000

long long arr[MAXN];

void test_sort(const char *name,
               std::function<void(long long[], int)> sort_algo) {
  Random rng(140384);
  for (int i = 0; i < MAXN; i++) {
    arr[i] = rng.nextLong();
  }

  double sort_time = time_it([&]() { sort_algo(arr, MAXN); });

  long long chk = 0;
  for (int i = 0; i < MAXN; i++) {
    chk = chk * 13 + arr[i];
  }
  fprintf(stderr, "%s in %.6lf, chk = %lld\n", name, sort_time, chk);
  assert(chk == 894150326305700963LL);
}

int main() {
  test_sort("std::sort", [](long long arr[], int n) { sort(arr, arr + n); });
  test_sort("ctree_sort", [](long long arr[], int n) { ctree_sort(arr, n); });
  test_sort("ska_sort", [](long long arr[], int n) { ska_sort(arr, arr + n); });
}
