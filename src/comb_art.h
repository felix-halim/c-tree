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

#include "art.h"

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

int n_buckets; // Number of leaf buckets.

// Dynamically resize COMB bucket sizes.
// If number of cracks > 32, it splits to two smaller buckets.
class CrackBucket {
  static const unsigned short MAX_CRACK = 64;

  int N;        // Number of data elements in this bucket pointed by D.
  int I;                // last indexed position
  unsigned char nC;     // the number of cracker indices
  unsigned long long S; // sorted bits
  int C[MAX_CRACK-1];   // the cracker indices
  int V[MAX_CRACK-1];     // the cracker value
  int D[BUCKET_SIZE];      // the data elements
  int n_erase;  // number of erase operations performed to this bucket.
  int n_touch;
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

  void piece_set_sorted(int i, bool sorted) {
    assert(i >= 0 && i < MAX_CRACK);
    if (sorted) {
      S |= 1ULL << i;
    } else {
      S &= ~(1ULL << i);
    }
  }

  bool piece_is_sorted(int i) const {
    assert(i >= 0 && i < MAX_CRACK);
    return S & (1ULL << i);
  }

  void piece_set_unsorted_onwards(int i) {
    assert(i >= 0 && i < MAX_CRACK);
    S &= (1ULL << i) - 1;          // destroy sorted bit std::vector from i onwards
  }

  void insert_bit_at(unsigned long long &S, int at) {
    S = ((S<<1) & ~((((1ULL<<at)-1)<<1)|1)) | (S & ((1ULL<<at)-1));
  }

  void remove_bit_at(unsigned long long &S, int at) {
    S = ((S & ~((((1ULL<<at)-1)<<1)|1)) >> 1) | (S & ((1ULL<<at)-1));
  }

  void add_cracker_index(int at, int M) {
    assert(at >= 0 && at <= nC && nC < MAX_CRACK - 1);
    for (int i = nC-1; i >= at; i--) {
      C[i + 1] = C[i];
      V[i + 1] = V[i];
    }
    C[at] = M;
    V[at] = D[M];
    nC++;
    assert(nC < MAX_CRACK);
    assert(at == 0 || C[at - 1] < C[at]);
    assert(at + 1 == nC || C[at] < C[at + 1]);
    insert_bit_at(S, at);
  }

  void remove_cracker_index(int at) {
    assert(at >= 0 && at < nC);
    for (int i = at + 1; i < nC; i++) {
      C[i - 1] = C[i];
      V[i - 1] = V[i];
    }
    nC--;
    // assert(nC>=0);
    remove_bit_at(S,at);
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

  void flush_pending_inserts() {
    // fprintf(stderr, "I = %d, N = %d\n", I, N);
    assert(I <= this->N);                     // Indexed index should be less than the number of elements
    if (!nC){ I = this->N; S = 0; return; }   // no index yet, all the elements are considered "inserted"
    assert(!next_b && !tail_b);               // Indexes only makes sense when there is no chain

    // IMPROVE: bulk insert? (Currently using Merge Completely)
    int minC = nC;
    for (int j = I; I < this->N; j= ++I) {        // insert all pending elements (from I to N)
      int i = nC - 1;
      int tmp = D[j];            // store the pending tuple
      for (; i>=0 && (tmp < V[i]); i--){  // insert by shuffling through the cracker indices C
        int &L = C[i];          // left boundary of this cracker piece
        D[j] = D[L+1];          // replace the pending with the next to cracker boundary
        D[L+1] = D[L];          // shift the cracker boundary to the right
        j = L++;            // reposition the cracker piece separator
      }
      D[j] = tmp;              // the pending tuple is now merged in
      minC = std::min(minC, i+1);      // keep track the lowest piece that is touched
    }

    piece_set_unsorted_onwards(minC);
  }


  // returns a piece [L,R) that contains v
  // it will reorganize the elements so that the R-L range gets smaller overtime
  int get_piece_by_value(int v, int &L, int &R, Random &rng) {
    flush_pending_inserts();
    int i = 0;
    while (i<nC && (v >= V[i])) i++;      // find the cracker indices that covers v
    L = i==0? 0 : C[i-1];            // the left crack boundary
    R = i==nC? this->N : C[i];            // the right crack boundary
    while (R-L > CRACK_AT){            // narrow down the piece using DDR
      int M = rough_middle_partition(D+L+(i?1:1), D+R, rng, DECRACK_AT) - D;
      add_cracker_index(i, M);
//        fprintf(stderr,"CRACKING %d %d, [%d %d]\n",M,D[M],L,R);
      if ((v < D[M])) R=M; else L=M, i++;  // adjust the cracker index i
    }
    assert(i>=0 && i<=nC);
//      assert(check(D[0],false,D[0],false));
    return i;
  }

public:

  CrackBucket(): next_b(nullptr), tail_b(nullptr) {
    this->N = 0;
    clear_indexes();
    n_erase = n_touch = 0;
  }

  CrackBucket(int *arr, int n): next_b(nullptr), tail_b(nullptr) {
    memcpy(D, arr, sizeof(int) * n);
    this->N = n;
    n_erase = n_touch = 0;
  }

  int n_cracks() { return nC; }
  int size() const { return N; }
  int slack() const { return BUCKET_SIZE - this->N; }
  void set_next(CrackBucket *b){ next_b = b; }
  void set_tail(CrackBucket *b){ tail_b = b; }
  void clear_indexes(){ S = nC = I = 0; }
  int touch() { return ++n_touch; }
  int n_updates() { return n_erase; }
  int remove_first() { assert(N > 0); int ret = D[0]; D[0] = D[--N]; return ret; }
  int capacity() const { return BUCKET_SIZE; }
  CrackBucket* next() const { return next_b; }
  CrackBucket* tail() const { return tail_b; }
  int& data(int i) { assert(i >= 0 && i < this->N); return D[i]; }
  void sort() { std::sort(D, D + N); };
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
    if (this->N < BUCKET_SIZE) {
      D[this->N++] = v;
      return;
    }
    assert(!tail_b || !tail_b->next());
    if (!(!next_b || tail_b)) fprintf(stderr, "%p %p\n", next_b, tail_b);
    assert(!next_b || tail_b);
    if (!next_b || tail_b->N == tail_b->capacity()) {
      // fprintf(stderr, "c");
      add_chain(new CrackBucket());
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
    assert(this->N < BUCKET_SIZE);
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
    clear_indexes();
    assert(this->N > 0 && this->N <= BUCKET_SIZE);
    return std::partition(D, D + this->N, [&](int x){ return (x < v); }) - D;
  }

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified CrackBucket "to", destroying all cracker indices
  void moveToFromIdx(CrackBucket *to, int fromIdx) {
    clear_indexes(); to->clear_indexes();    // destroy both buckets' cracker indices
    assert(this->N > fromIdx);            // make sure there is something to move
    assert(to->N + this->N - fromIdx <= BUCKET_SIZE);    // make sure the receiver has enough space
    memmove(to->D + to->N, D+fromIdx, (this->N - fromIdx) * sizeof(int));
    to->N += this->N - fromIdx;
    this->N = fromIdx;
  }

  int lower_pos(int value, Random &rng) {
    int i, L, R;
    return crack(value, i, L, R, true, rng);
  }

  int crack(int v, int &i, int &L, int &R, bool sort_piece, Random &rng) {
    assert(!next());            // it doesn't make sense crack a chained bucket!
    i = get_piece_by_value(v,L,R,rng);    // find the piece [L,R) containing v
    assert(L>=0 && L<=R && R <= this->N);        // range check
    if (!piece_is_sorted(i)){
      if (sort_piece){            // sort the piece if requested
        std::sort(D+L,D+R);
        piece_set_sorted(i,true);
      } else {
        for (int at=L; at<R; at++)
          // if (eq(D[at],v)) return at;
          if (D[at] == v) return at;
        return R;
      }
    }
    for (int pos = L; pos < R; pos++) if ((D[pos] >= v)) return pos;
    return R;
    // return std::lower_bound(D+L, D+R, v) - D;    // find the element v using binary search
  }

  int erase(int v, Random &rng) {
    n_erase++;
    assert(!next_b);
    int i, L, R, at = crack(v, i, L, R, false, rng);
    // if (at >= R || !eq(D[at],v)) return false;  // the element to be erased is not found!
    if (at >= R) {
      fprintf(stderr, "R = %d, N = %d, v = %d, DR = %d\n", R, N, v, D[R]);
      for (int i = 0; i < N; i++) {
        fprintf(stderr, "D[%d] = %d\n", i, D[i]);
      }
      assert(R == N);
      return at;
    }
    if (D[at] != v) return -1;  // the element to be erased is not found!

    // decrack this cracker piece (it becomes too small) or
    // if the deleted element index is a cracker index
    if (nC && R-L+1 <= DECRACK_AT){
      remove_cracker_index((i>0)?(--i):i);
    } else if (i > 0 && at == L){
      at = L+1;
      for (int j=L+2; j<R; j++)    // find a replacement element for the cracker index
        if ((D[j] < D[at])) at = j;  // that is the smallest element in the piece (L,R)
      std::swap(D[L], D[at]);
      piece_set_sorted(i,false);
      V[i-1] = D[L];
    }

    assert(at<R && (D[at] == v));    // the element v must be found!

    // IMPROVE: use pending delete? antimatter?
    piece_set_unsorted_onwards(i);      // unset the sorted bit i onwards
    assert(i==nC || (D[at] < V[i]));
    for (int j=i; j<nC; j++){        // shuffle out the deleted element
      R = C[j]--;
      D[at] = D[R-1];
      D[R-1] = D[R];
      at = R;
    }
    D[at] = D[--this->N];    // the deleted element has been shuffled out from the bucket
    I--;        // adjust the pending index

    if (at == 0) {
      assert(0);
      R = 0==nC? this->N : C[0];            // the right crack boundary
      nth_element(D, D, D + R);
    }
//      assert(check(D[0],false,D[0],false));
    return at;
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
          Rb = new B();
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


  void insert_root(CrackBucket *b) {
    uint8_t key[8];
    assert(b->data(0) >= 0);
    loadKey(b->data(0), key);
    n_buckets++;
    // fprintf(stderr, "add root %d\n", b->data(0));
    ::insert(tree, &tree, key, 0, (uintptr_t) b, 8);
  }

  void insert_root_value(uintptr_t value) {
    uint8_t key[8];
    assert(value >= 0);
    assert(value < (1ULL << 40));
    // fprintf(stderr, "insert root value %llu\n", (unsigned long long) value);
    loadKey(value, key);
    ::insert(tree, &tree, key, 0, (value << 1) | 1, 8);
  }

  void load(int const *arr, int n) {
    CrackBucket *root = new CrackBucket();
    root->append(0); // Dummy leftmost bucket.
    assert(root);
    assert(!( ((uintptr_t) root) & 3));
    insert_root(root);

    int i = 0;
    while (i + BUCKET_SIZE <= n) {
      CrackBucket *b = new CrackBucket();
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

  Node* find_bucket_upper(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    return ::lower_bound(tree,key,8,0,8);
  }

  void insert(int const &value) {
    // fprintf(stderr, "ins %d\n", value);
    Node *n = find_bucket(value);
    assert(isLeaf(n));
    uintptr_t v = getData((uintptr_t) n);
    if (isPointer(v)) {
      ((CrackBucket*) v)->insert(value);
    } else {
      insert_root_value(value);
    }
  }

  void transition_to_art(CrackBucket *b) {
    // fprintf(stderr, "a");

    // fprintf(stderr, "transdel\n");
    bool ok = erase_root(b->data(0));
    assert(ok);
    // fprintf(stderr, "completing %d\n", b->size());
    // fprintf(stderr, "b");

    b->sort();
    for (int i = 0; i < b->size(); i++) {
      insert_root_value(b->data(i));
    }
    // fprintf(stderr, "c");

    // fprintf(stderr, "completing2 %d\n", b->size());
    delete b;
    n_buckets--;
    assert(n_buckets >= 0);
    if (!n_buckets) fprintf(stderr, "YAY!\n");
    // fprintf(stderr, "d");
  }

  bool erase_root(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    // assert(lookup(&tree,key,8,0,8));
    // fprintf(stderr, "erase root %llu\n", value64);
    ::erase(tree,&tree,key,8,0,8);
    // fprintf(stderr, "erase root done %llu\n", value64);
    return true;
  }

  CrackBucket* make_standalone(CrackBucket *lb, int value) {
    while (lb->next()) {
      CrackBucket *rb = stochastic_split_chain(lb, rng);
      // fprintf(stderr, "stochastic_split_chain %d %d\n", lb->size(), rb->size());
      insert_root(rb);
      lb = (value < rb->data(0)) ? lb : rb;
    }
    return lb;
  }

  bool erase(int const &value) {
    fprintf(stderr, "ERASE %d\n", value);
    // art_debug = 1;
    if (n_buckets) {
      Node *n = find_bucket(value);
      assert(isLeaf(n));
      uintptr_t v = getData((uintptr_t) n);
      if (isPointer(v)) {
        CrackBucket *lb = (CrackBucket*) v;
        int lo = lb->data(0);
        lb = make_standalone(lb, value);
        assert(lo == lb->data(0));
        if (value == lo) {

        }
      } else {
        assert(0);
        assert(getData(v) == (uintptr_t) value);
        return erase_root(value);
      }
    }
    assert(0);
    return erase_root(value);
  }

  /* TODO: lazy lower_bound */
  int lower_bound(int const &value) {
    static int nth = 0; nth++;
    // art_debug = 1;
    // fprintf(stderr, "lower_bound %d\n", value);

    if (n_buckets) {
      Node *n = find_bucket(value);
      assert(isLeaf(n));
      uintptr_t v = getData((uintptr_t) n);
      assert(v);
      if (isPointer(v)) {
        CrackBucket *lb = make_standalone((CrackBucket*) v, value);
        assert(lb);
        int pos = lb->lower_pos(value, rng);
        if (pos < lb->size()) {
          int ret = lb->data(pos);
          if (lb->n_cracks() > 30 && lb->touch() > 10) transition_to_art(lb);
          return ret;
        }
      }
    }

    uint64_t value64 = value;
    uint8_t key[8];
    loadKey(value64, key);
    // fprintf(stderr, "a");
    Node* leaf = ::lower_bound(tree,key,8,0,8);
    // fprintf(stderr, "b %p", leaf);
    return isLeaf(leaf) ? getLeafValue(leaf) : 0;
  }
};

#endif
