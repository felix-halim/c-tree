#include <cassert>
#include <cstdio>

#include <functional>
#include <algorithm>
#include <vector>

#include "random.h"
#include "time_it.h"

using namespace std;

Random rng(140384);

#define BSIZE 4096

class Bucket {
public:
  long long arr[BSIZE];
  int n;
  Bucket *next, *tail;

  Bucket():n(0), next(nullptr), tail(nullptr) {}

  long long* data() { return arr; }
  int size() { return n; }

  long long remove_rng_data() {
    int j = rng.nextInt(this->n);
    long long ret = arr[j];
    arr[j] = arr[--this->n];
    return ret;
  }

  int slack() const { return BSIZE - this->n; }

  void moveToFromIdx(Bucket *to, int fromIdx) {
    assert(this->n > fromIdx);            // make sure there is something to move
    assert(to->n + this->n - fromIdx <= BSIZE);    // make sure the receiver has enough space
    memmove(to->arr + to->n, arr+fromIdx, (this->n - fromIdx) * sizeof(long long));
    to->n += this->n - fromIdx;
    this->n = fromIdx;
  }

  void add_chain(Bucket *b) {
    if (next) {
      assert(tail);
      tail->next = b;
    } else {
      next = b;
    }
    tail = b;
    b->next=(nullptr);
    b->tail=(nullptr);
  }

  int partition(long long const &v) {
    assert(this->n > 0 && this->n <= BSIZE);
    return int(std::partition(arr, arr + this->n, [&](long long x){ return (x < v); }) - arr);
  }
};

long long get_rng_pivot(Bucket *b) {  // pick the pivot near the median
  return b->arr[rng.nextInt(b->n)];
}

void mark_hi(long long* D, int N, long long const &P, int *hi, int &nhi){
  for (int i=0; i < N; i++){
    hi[nhi] = i;
    nhi += (D[i] >= P);
  }
}

void mark_lo(long long* D, int N, long long const &P, int *lo, int &nlo){
  for (int i=0; i < N; i++){
    lo[nlo] = i;
    nlo += (D[i] < P);
  }
}

void fusion(long long *Lp, long long *Rp, int *hi, int *lo, int &nhi, int &nlo){
  int m = std::min(nhi, nlo); assert(m > 0);
  int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
  nhi -= m; nlo -= m;
  while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
}


bool add_to_chain(Bucket *&chain, Bucket *b) {
  if (!b->size()) {
    delete b;
    return true;
  }

  if (!chain) { // the root chain is empty.
    chain = b;  // b is the head of the chain.
    chain->next=(nullptr);
    chain->tail=(nullptr);
  } else {
    assert(!chain->tail || !((Bucket*) chain->tail)->next);
    if (chain->slack()) {
      b->moveToFromIdx(chain, b->size() - std::min(b->size(), chain->slack()));
    } else if (chain->tail) {
      if (((Bucket*) chain->tail)->slack()) {
        b->moveToFromIdx((Bucket*) chain->tail, b->size() - std::min(b->size(), ((Bucket*) chain->tail)->slack()));
      }
    }
    if (b->size()) {
      chain->add_chain(b);
    } else {
      delete b;
      return true;
    }
  }
  return false;
}

pair<Bucket*, Bucket*> split_chain(Bucket* b) {
  long long p = get_rng_pivot(b);

  Bucket *left_chain = nullptr;
  Bucket *right_chain = nullptr;
  Bucket *Lb = nullptr, *Rb = nullptr;
  int hi[BSIZE], lo[BSIZE];
  int nhi = 0, nlo = 0;

  while (true) {
    if (nhi && nlo) {
      assert(Lb && Rb);
      fusion(Lb->data(), Rb->data(), hi, lo, nhi, nlo);
      if (!nhi) { add_to_chain(left_chain, Lb); Lb = nullptr; }
      if (!nlo) { add_to_chain(right_chain, Rb); Rb = nullptr; }
    } else if (!Lb) {
      if (!b) break;
      Lb = b;
      b = (Bucket*) b->next;
    } else if (!nhi) {
      assert(Lb);
      mark_hi(Lb->data(), Lb->size(), p, hi, nhi);
      if (!nhi) { add_to_chain(left_chain, Lb); Lb = nullptr; }
    } else if (!Rb) {
      if (!b) break;
      Rb = b;
      b = (Bucket*) b->next;
    } else if (!nlo) {
      assert(Rb);
      mark_lo(Rb->data(), Rb->size(), p, lo, nlo);
      if (!nlo) { add_to_chain(right_chain, Rb); Rb = nullptr; }
    } else {
      assert(0);
    }
  }

  if (Rb) { assert(!Lb); Lb = Rb; }
  if (Lb) {
    if (Lb->size()) {
      int i = Lb->partition(p);
      if (i == 0) {
        add_to_chain(right_chain, Lb);
      } else if (i == Lb->size()) {
        add_to_chain(left_chain, Lb);
      } else {
        Rb = new Bucket();
        Lb->moveToFromIdx(Rb, i);
        add_to_chain(left_chain, Lb);
        add_to_chain(right_chain, Rb);
      }
    } else {
      delete Lb;
    }
  }

  // assert(left_chain);
  assert(right_chain);
  return make_pair(left_chain, right_chain);
}

void ctree_sort(long long arr[], int N) {
  Bucket *head = nullptr;
  int nbuckets = 0;
  for (int i = 0; i < N; ) {
    int cnt = min(BSIZE, N - i);
    if (!head) {
      head = new Bucket();
      head->tail = head;
    } else {
      head->tail->next = new Bucket();
      head->tail = head->tail->next;
    }
    memcpy(head->tail->arr, arr + i, sizeof(long long) * cnt);
    head->tail->n = cnt;
    i += cnt;
    nbuckets++;
  }
  fprintf(stderr, "nbuckets = %d\n", nbuckets);

  vector<Bucket*> stk;
  stk.push_back(head);
  while (stk.size()) {
    Bucket* chain = stk.back();
    stk.pop_back();
    if (!chain->next) {
      sort(chain->arr, chain->arr + chain->n);
      memcpy(arr + N - chain->n, chain->arr, sizeof(long long) * chain->n);
      N -= chain->n;
      delete chain;
      continue;      
    }

    auto p = split_chain(chain);
    if (p.first) {
      stk.push_back(p.first);
    }
    if (p.second) {
      stk.push_back(p.second);
    }
  }
  assert(N == 0);
}

#define MAXN 100000000

long long arr[MAXN];

void test_sort(const char* name, std::function<void(long long[], int)> sort_algo) {
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
  test_sort("std::sort", [](long long arr[], int n){ sort(arr, arr + n); });
  test_sort("ctree_sort", [](long long arr[], int n){ ctree_sort(arr, n); });
}
