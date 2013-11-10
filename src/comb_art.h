#ifndef _COMB_H_
#define _COMB_H_

#include <cassert>
#include <cstdio>
#include <cstring>

#include <functional>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

#include "art_best.h"

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

#define CRACK_AT (BUCKET_SIZE >> 5)
#define DECRACK_AT (BUCKET_SIZE >> 6)

// Dynamically resize COMB bucket sizes.
// If number of cracks > 32, it splits to two smaller buckets.
class CrackBucket {
  int N;        // Number of data elements in this bucket pointed by D.
  int cap;      // last indexed position
  int *D;      // the data elements
  CrackBucket* next_b;   // buckets can be chained like a linked list of buckets
                        // the value of next is -1 if there is no next chain
                        // otherwise the index of the bucket [0, num_of_buckets)
  CrackBucket* tail_b;   // pointer to the last bucket in the chain.

  int* partition(int *F, int *L, int const &v) {
    while (true) {
      while (true)
        if (F == L) return F;
        else if (*F < v) ++F;
        else break;
      --L;
      while (true)
        if (F == L) return F;
        else if (!bool(*L < v)) --L;
        else break;
      std::iter_swap(F, L);
      ++F;
    }
    assert(false);
  }

  // partitions roughly in the middle satisfying the DECRACK_AT
  int* rough_middle_partition(int *L, int *R, Random &rng, int gap) {
    int ntry = 10;
    for (int *i=L, *j=R, *p; ntry--; ){
      std::iter_swap(j-1, i + rng.nextInt(j-i));
      std::iter_swap(j-1, p = partition(i, j-1, *(j-1)));
      if (p-L <= gap) i = p;
      else if (R-p <= gap) j = p;
      else return p;
    }
    int *M = L+((R-L)>>1);
    std::nth_element(L,M,R);
    return M;
  }

public:
  CrackBucket(int c): cap(c), next_b(nullptr), tail_b(nullptr) {
    D = new int[c];
    this->N = 0;
  }

  CrackBucket(int *arr, int n): cap(n), next_b(nullptr), tail_b(nullptr) {
    D = new int[n];
    memcpy(D, arr, sizeof(int) * n);
    this->N = n;
  }

  ~CrackBucket() {
    delete[] D;
  }

  int size() const { return N; }
  int slack() const { return cap - this->N; }
  void set_next(CrackBucket *b){ next_b = b; }
  void set_tail(CrackBucket *b){ tail_b = b; }
  int remove_first() { assert(N > 0); int ret = D[0]; D[0] = D[--N]; return ret; }
  int capacity() const { return cap; }
  CrackBucket* next() const { return next_b; }
  CrackBucket* tail() const { return tail_b; }
  int& data(int i) { assert(i >= 0 && i < this->N); return D[i]; }
  void set_data(int i, int v) { assert(i >= 0 && i < this->N); D[i] = v; }
  int* data() { return D; }

  CrackBucket* split(Random &rng) {
    int M = rough_middle_partition(D + 1, D + N, rng, N / 3) - D;
    CrackBucket *rb = new CrackBucket(D + M, N - M);
    N = M;
    return rb;
  }

  void add_chain(CrackBucket *next) {
    if (next_b) {
      assert(tail_b);
      tail_b->next_b = next;
    } else {
      next_b = next;
    }
    tail_b = next;
    next->set_next(nullptr);
    next->set_tail(nullptr);
  }

  void insert(int const &v) {
    if (this->N < cap) {
      D[this->N++] = v;
      return;
    }
    assert(!tail_b || !tail_b->next());
    if (!(!next_b || tail_b)) fprintf(stderr, "%p %p\n", next_b, tail_b);
    assert(!next_b || tail_b);
    if (!next_b || tail_b->N == tail_b->capacity()) {
      // fprintf(stderr, "c");
      add_chain(new CrackBucket(cap));
      // fprintf(stderr, "d");
    }
    assert(tail_b && tail_b->N < tail_b->capacity());
    tail_b->D[tail_b->N++] = v;
  }

  int remove_random_data(Random &rng) {
    int j = rng.nextInt(this->N);
    int ret = D[j];
    D[j] = D[--this->N];
    return ret;
  }

  void swap_random_data_with(int &R, Random &rng) {
    assert(this->N > 0);
    swap(R, D[rng.nextInt(this->N)]);
  }

  void append(int value) {
    assert(this->N < cap);
    D[this->N++] = value;
  }

  void each(std::function<void(int)> callback) {
    for (int i = 0; i < N; i++) {
      callback(D[i]);
    }
  }

  int bulk_insert(int const *v, int length) {
    assert(this->N == 0);
    memcpy(D, v, sizeof(int) * length);
    return this->N = length;
  }

  // partition this bucket based on value v, destroying all cracker indices
  int partition(int const &v) {
    assert(this->N > 0 && this->N <= cap);
    return std::partition(D, D + this->N, [&](int x){ return (x < v); }) - D;
  }

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified CrackBucket "to", destroying all cracker indices
  void moveToFromIdx(CrackBucket *to, int fromIdx) {
    assert(this->N > fromIdx);            // make sure there is something to move
    assert(to->N + this->N - fromIdx <= to->cap);    // make sure the receiver has enough space
    memmove(to->D + to->N, D+fromIdx, (this->N - fromIdx) * sizeof(int));
    to->N += this->N - fromIdx;
    this->N = fromIdx;
  }
};

inline bool isPointer(uintptr_t v) {
  return !(v & 1);
}

inline uintptr_t getData(uintptr_t v) {
  return v >> 1;
}

inline uintptr_t getLeafValue(Node* node) {
  // The the value stored in the pseudo-leaf
  uintptr_t v = getData((uintptr_t) node);
  return isPointer(v) ? ((CrackBucket*) v)->data(0) : getData(v);
}

class Comb {
  Random rng;  // The random number generator.
  Node* tree = NULL;

  // add a CrackBucket (bidx) to the root chain 'ridx'
  template <typename B>
  bool add_to_chain(B *&chain, B *b) {
    if (!b->size()) {
      delete b;
      return true;
    }

    if (!chain) { // the root chain is empty.
      chain = b;  // b is the head of the chain.
      chain->set_next(nullptr);
      chain->set_tail(nullptr);
    } else {
      assert(!chain->tail() || !((B*) chain->tail())->next());
      if (chain->slack()) {
        b->moveToFromIdx(chain, b->size() - std::min(b->size(), chain->slack()));
      } else if (chain->tail()) {
        if (((B*) chain->tail())->slack()) {
          b->moveToFromIdx((B*) chain->tail(), b->size() - std::min(b->size(), ((B*) chain->tail())->slack()));
        }
      }
      if (b->size()) {
        // fprintf(stderr, "a");
        chain->add_chain(b);
        // fprintf(stderr, "b");
      } else {
        delete b;
        return true;
      }
    }
    return false;
  }

  template <typename B>
  int get_random_pivot(B *b, Random &rng) {  // pick the pivot near the median
    assert(b->next());
    int R[11]; // Randomly pick 11 elements from b.
    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "%p sz = %d\n", b, b->size());
      if (b->size()) {
        R[i] = b->remove_random_data(rng);
      } else {
        b = b->next();
        assert(b);
        i--;
      }
    }

    // Replace R with the next buckets in the chain using reservoir sampling.
    B *Nb = (B*) b->next();
    for (int i = 1; Nb; i++) {
      if (!Nb->size()) {
        i--;
      } else {
        int j = rng.nextInt(i);
        if (j < 11) Nb->swap_random_data_with(R[j], rng);
      }
      Nb = (B*) Nb->next();
    }

    std::nth_element(R, R + 5, R + 11);
    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      if (i != 5) b->append(R[i]);
    }
    // fprintf(stderr, "pivot = %d\n", pivot);
    return R[5];
  }

  void mark_hi(int* D, int N, int const &P, int *hi, int &nhi){
    for (int i=0; i < N; i++){
      hi[nhi] = i;
      nhi += (D[i] >= P);
    }
  }

  void mark_lo(int* D, int N, int const &P, int *lo, int &nlo){
    for (int i=0; i < N; i++){
      lo[nlo] = i;
      nlo += (D[i] < P);
    }
  }

  void fusion(int *Lp, int *Rp, int *hi, int *lo, int &nhi, int &nlo){
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  // fusion
  template <typename B>
  B* stochastic_split_chain(B *b, Random &rng) {
    // fprintf(stderr, "bsz = %d\n", b->size());
    const int &bp = b->remove_first();
    const int &p = get_random_pivot(b, rng);

    B *left_chain = nullptr;
    B *right_chain = nullptr;
    B *Lb = nullptr, *Rb = nullptr;
    int hi[BUCKET_SIZE], lo[BUCKET_SIZE];
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
        b = (B*) b->next();
      } else if (!nhi) {
        assert(Lb);
        mark_hi(Lb->data(), Lb->size(), p, hi, nhi);
        if (!nhi) { add_to_chain(left_chain, Lb); Lb = nullptr; }
      } else if (!Rb) {
        if (!b) break;
        Rb = b;
        b = (B*) b->next();
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
          Rb = new B(Lb->capacity());
          Lb->moveToFromIdx(Rb, i);
          add_to_chain(left_chain, Lb);
          add_to_chain(right_chain, Rb);
        }
      } else {
        delete Lb;
      }
    }

    assert(left_chain);
    int first = left_chain->data(0);
    left_chain->set_data(0, bp);
    left_chain->insert(first);

    assert(right_chain);
    first = right_chain->data(0);
    right_chain->set_data(0, p);
    right_chain->insert(first);

//    assert(check());
    return right_chain;
  }

public:

  int n_buckets; // Number of leaf buckets.

  void insert_root(CrackBucket *b) {
    uint8_t key[8];
    assert(b->data(0) >= 0);
    loadKey(b->data(0), key);
    n_buckets++;
    // fprintf(stderr, "add root %d\n", b->data(0));
    ::insert(&tree, key, 0, (uintptr_t) b, 8);
  }

  void insert_root_value(uintptr_t value) {
    uint8_t key[8];
    assert(value >= 0);
    assert(value < (1ULL << 40));
    // fprintf(stderr, "insert root value %llu\n", (unsigned long long) value);
    loadKey(value, key);
    ::insert(&tree, key, 0, (value << 1) | 1, 8);
  }

  void load(int const *arr, int n) {
    CrackBucket *root = new CrackBucket(BUCKET_SIZE);
    root->append(0); // Dummy leftmost bucket.
    assert(root);
    assert(!( ((uintptr_t) root) & 3));
    insert_root(root);

    int i = 0;
    while (i + BUCKET_SIZE <= n) {
      CrackBucket *b = new CrackBucket(BUCKET_SIZE);
      b->bulk_insert(arr + i, BUCKET_SIZE);
      i += BUCKET_SIZE;
      if (root) {
        ((CrackBucket*) root)->add_chain(b);
      } else {
        root = b;
      }
    }

    while (i < n) {
      insert(arr[i++]);
    }
  }

  Node* find_bucket(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    Node* ret = lower_bound_prev(tree,key,8,0,8);
    assert(ret);
    assert(isLeaf(ret));
    return ret;
  }

  void insert(int const &value) {
    // fprintf(stderr, "ins %d\n", value);
    uintptr_t v = getData((uintptr_t) find_bucket(value));
    if (isPointer(v)) {
      ((CrackBucket*) v)->insert(value);
    } else {
      insert_root_value(value);
    }
  }

  CrackBucket* transition_to_art(uintptr_t v) {
    assert(isPointer(v));
    // Transition to root array.
    // fprintf(stderr, "T");
    CrackBucket *b = (CrackBucket*) v;

    if (b->size() < 128) {
      bool ok = erase_root(b->data(0));
      assert(ok);
      // fprintf(stderr, "completing %d\n", b->size());
      b->each([&](int x) { insert_root_value(x); });
      // fprintf(stderr, "completing2 %d\n", b->size());
      delete b;
      n_buckets--;
      assert(n_buckets >= 0);
      if (!n_buckets) fprintf(stderr, "YAY!\n");
      return nullptr;
    } else {
      // fprintf(stderr, "splitting\n");
      return b->split(rng);
    }
  }

  void artify(int value) {
    uintptr_t v = getData((uintptr_t) find_bucket(value));
    if (isPointer(v)) {
      CrackBucket *lb = (CrackBucket*) v;
      while (lb->next()) {
        CrackBucket *rb = stochastic_split_chain(lb, rng);
        // fprintf(stderr, "stochastic_split_chain %d %d\n", lb->size(), rb->size());
        insert_root(rb);
        lb = (value < rb->data(0)) ? lb : rb;
      }
      while (true) {
        // fprintf(stderr, "trans %p, %d\n", lb, lb->size());
        CrackBucket *rb = transition_to_art((uintptr_t) lb);
        if (!rb) break;
        // fprintf(stderr, "trans %p, %d %d\n", rb, lb->size(), rb->size());
        insert_root(rb);
        lb = (value < rb->data(0)) ? lb : rb;
      }
    }
  }

  bool erase_root(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    // assert(lookup(&tree,key,8,0,8));
    return ::erase(tree,&tree,key,8,0,8);
  }

  bool erase(int const &value) {
    // fprintf(stderr, "ERASE %d\n", value);
    // art_debug = 1;
    if (n_buckets) artify(value);
    return erase_root(value);
  }

  /* TODO: lazy lower_bound */
  int lower_bound(int const &value) {
    static int nth = 0; nth++;
    // art_debug = 1;
    // if (nth % 1000 == 0) assert(check());
    // fprintf(stderr, "lower_bound %d\n", value);
    if (n_buckets) artify(value);

    uint64_t value64 = value;
    uint8_t key[8];
    loadKey(value64, key);

    Node* leaf = ::lower_bound(tree,key,8,0,8);
    return isLeaf(leaf) ? getLeafValue(leaf) : 0;
  }
};

#endif
