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

#include "random"
#include "art.h"

#define CRACK_AT (LARGE_SIZE >> 5)
#define DECRACK_AT (LARGE_SIZE >> 6)

int n_large; // Number of leaf large buckets.
int n_small; // Number of leaf small buckets.
int n_index; // Number of buckets chains.

template <typename T>
class Bucket {
 protected:
  int N;        // Number of data elements in this bucket pointed by D.
 public:
  int large_type;
  int size() { return N; }
};

template <typename T>
class SmallBucket : public Bucket<T> {
  T D[SMALL_SIZE];
  int n_touch;

 public:
  SmallBucket(T *arr, int n) {
    this->N = n;
    memcpy(D, arr, sizeof(T) * n);
    sort(D, D + this->N);
    this->large_type = 0;
    n_small++;
  }

  ~SmallBucket() {
    n_small--;
  }

  int slack() { return SMALL_SIZE - this->N; }

  int n_touched() { return ++n_touch; }

  T data(int i) {
    // fprintf(stderr, "i = %d / %d\n", i, N);
    assert(i >= 0 && i < this->N); return D[i];
  }

  int lower_pos(T v) {
    for (int i = 0; i < this->N; i++) {
      if (D[i] >= v) return i;
    }
    return this->N;
  }

  int erase(T v) {
    assert(this->N);
    int i = 0;
    for (; i < this->N; i++) {
      if (D[i] == v) break;
    }
    if (i == this->N) return -1;
    this->N--;
    for (int j = i; j < this->N; j++) {
      D[j] = D[j + 1];
    }
    return i;
  }

  bool insert(T &v) {
    if (this->N == SMALL_SIZE) {
      if (D[this->N - 1] < v) return false;
      T t = D[this->N - 1];
      for (int i = this->N - 1; i > 0; i--) {
        if (D[i - 1] > v) {
          D[i] = D[i - 1];
        } else {
          D[i] = v;
          v = t;
          return false;
        }
      }
      assert(0);
    }
    int i = this->N++;
    for (; i > 0; i--) {
      if (D[i - 1] > v) {
        D[i] = D[i - 1];
      } else break;
    }
    assert(i);
    D[i] = v;
    return true;
  }
};

// Dynamically resize COMB bucket sizes.
// If number of cracks > 32, it splits to two smaller buckets.
template <typename T>
class LargeBucket : public Bucket<T> {
  static const unsigned short MAX_CRACK = 64;

  int I;                // last indexed position
  unsigned char nC;     // the number of cracker indices
  unsigned long long S; // sorted bits
  int C[MAX_CRACK-1];   // the cracker indices
  T V[MAX_CRACK-1];     // the cracker value
  T D[LARGE_SIZE];      // the data elements
  int n_erase;  // number of erase operations performed to this bucket.
  int n_touch;
  LargeBucket* next_b;   // buckets can be chained like a linked list of buckets
                        // the value of next is -1 if there is no next chain
                        // otherwise the index of the bucket [0, num_of_buckets)
  LargeBucket* tail_b;   // pointer to the last bucket in the chain.

  T* partition(T *F, T *L, T const &v) {
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
  T* rough_middle_partition(T *L, T *R, Random &rng, int gap) {
    int ntry = 10;
    for (T *i=L, *j=R, *p; ntry--; ){
      std::iter_swap(j-1, i + rng.nextInt(j-i));
      std::iter_swap(j-1, p = partition(i, j-1, *(j-1)));
      if (p-L <= gap) i = p;
      else if (R-p <= gap) j = p;
      else return p;
    }
    T *M = L+((R-L)>>1);
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
      T tmp = D[j];            // store the pending tuple
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
  int get_piece_by_value(T v, int &L, int &R, Random &rng) {
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

  LargeBucket(): next_b(nullptr), tail_b(nullptr) {
    this->N = 0;
    clear_indexes();
    n_erase = n_touch = 0;
    this->large_type = 1;
    n_large++;
  }

  LargeBucket(T *arr, int n): next_b(nullptr), tail_b(nullptr) {
    memcpy(D, arr, sizeof(T) * n);
    this->N = n;
    n_erase = n_touch = 0;
    this->large_type = 1;
    n_large++;
  }

  ~LargeBucket() {
    n_large--;
  }

  int n_cracks() { return nC; }
  int slack() const { return LARGE_SIZE - this->N; }
  void set_next(LargeBucket *b){ next_b = b; }
  void set_tail(LargeBucket *b){ tail_b = b; }
  void clear_indexes(){ S = nC = I = 0; }
  int n_touched() { return ++n_touch; }
  int n_erased() { return ++n_erase; }
  T remove_first() { assert(this->N > 0); T ret = D[0]; D[0] = D[--this->N]; return ret; }
  int capacity() const { return LARGE_SIZE; }
  LargeBucket* next() const { return next_b; }
  LargeBucket* tail() const { return tail_b; }
  T& data(int i) { assert(i >= 0 && i < this->N); return D[i]; }
  void sort() { std::sort(D, D + this->N); };
  void set_data(int i, T v) { assert(i >= 0 && i < this->N); D[i] = v; }
  T* data() { return D; }

  void rec_split(int L, int R, Random &rng, std::function<void(T*, int)> callback) {
    if (R - L <= SMALL_SIZE) {
      callback(D + L, R - L);
    } else {
      int M = rough_middle_partition(D + L + 1, D + R, rng, (R - L) / 3) - D;
      rec_split(L, M, rng, callback);
      rec_split(M, R, rng, callback);
    }
  }

  void split(Random &rng, std::function<void(T*, int)> callback) {
    rec_split(0, this->N, rng, callback);
  }

  void add_chain(LargeBucket *next) {
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

  void insert(T const &v) {
    if (this->N < LARGE_SIZE) {
      D[this->N++] = v;
      return;
    }
    assert(!tail_b || !tail_b->next());
    if (!(!next_b || tail_b)) fprintf(stderr, "%p %p\n", next_b, tail_b);
    assert(!next_b || tail_b);
    if (!next_b || tail_b->N == tail_b->capacity()) {
      // fprintf(stderr, "c");
      add_chain(new LargeBucket());
      // fprintf(stderr, "d");
    }
    assert(tail_b && tail_b->N < tail_b->capacity());
    tail_b->D[tail_b->N++] = v;
  }

  T remove_random_data(Random &rng) {
    int j = rng.nextInt(this->N);
    T ret = D[j];
    D[j] = D[--this->N];
    return ret;
  }

  void swap_random_data_with(T &R, Random &rng) {
    assert(this->N > 0);
    swap(R, D[rng.nextInt(this->N)]);
  }

  void append(T value) {
    assert(this->N < LARGE_SIZE);
    D[this->N++] = value;
  }

  void each(std::function<void(T)> callback) {
    for (int i = 0; i < this->N; i++) {
      callback(D[i]);
    }
  }

  int bulk_insert(T const *v, int length) {
    assert(this->N == 0);
    memcpy(D, v, sizeof(T) * length);
    return this->N = length;
  }

  // partition this bucket based on value v, destroying all cracker indices
  int partition(T const &v) {
    clear_indexes();
    assert(this->N > 0 && this->N <= LARGE_SIZE);
    return std::partition(D, D + this->N, [&](T x){ return (x < v); }) - D;
  }

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified LargeBucket "to", destroying all cracker indices
  void moveToFromIdx(LargeBucket *to, int fromIdx) {
    clear_indexes(); to->clear_indexes();    // destroy both buckets' cracker indices
    assert(this->N > fromIdx);            // make sure there is something to move
    assert(to->N + this->N - fromIdx <= LARGE_SIZE);    // make sure the receiver has enough space
    memmove(to->D + to->N, D+fromIdx, (this->N - fromIdx) * sizeof(T));
    to->N += this->N - fromIdx;
    this->N = fromIdx;
  }

  int lower_pos(T value, Random &rng) {
    int i, L, R;
    return crack(value, i, L, R, true, rng);
  }

  int crack(T v, int &i, int &L, int &R, bool sort_piece, Random &rng) {
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

  int erase(T v, Random &rng) {
    assert(!next_b);
    int i, L, R, at = crack(v, i, L, R, false, rng);
    // if (at >= R || !eq(D[at],v)) return false;  // the element to be erased is not found!
    if (at >= R) {
      fprintf(stderr, "R = %d, N = %d, v = %d, DR = %d\n", R, this->N, v, D[R]);
      for (int i = 0; i < this->N; i++) {
        fprintf(stderr, "D[%d] = %d\n", i, D[i]);
      }
      assert(R == this->N);
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

template <typename T>
T data(Bucket<T> *b) {
  return b->large_type ? ((LargeBucket<T>*) b)->data(0) : ((SmallBucket<T>*) b)->data(0);
}

template <typename T>
class Comb {
  Random rng;  // The random number generator.
  Node* tree;

  // add a LargeBucket (bidx) to the root chain 'ridx'
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
  T get_random_pivot(B *b, Random &rng) {  // pick the pivot near the median
    assert(b->next());
    T R[11]; // Randomly pick 11 elements from b.
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

  void mark_hi(T* D, int N, T const &P, int *hi, int &nhi){
    for (int i=0; i < N; i++){
      hi[nhi] = i;
      nhi += (D[i] >= P);
    }
  }

  void mark_lo(T* D, int N, T const &P, int *lo, int &nlo){
    for (int i=0; i < N; i++){
      lo[nlo] = i;
      nlo += (D[i] < P);
    }
  }

  void fusion(T *Lp, T *Rp, int *hi, int *lo, int &nhi, int &nlo){
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  // fusion
  template <typename B>
  B* stochastic_split_chain(B *b, Random &rng) {
    // fprintf(stderr, "bsz = %d\n", b->size());
    const T &bp = b->remove_first();
    const T &p = get_random_pivot(b, rng);

    B *left_chain = nullptr;
    B *right_chain = nullptr;
    B *Lb = nullptr, *Rb = nullptr;
    int hi[LARGE_SIZE], lo[LARGE_SIZE];
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
    T first = left_chain->data(0);
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

  Comb() { tree = NULL; }

  void insert_root(Bucket<T> *b) {
    uint8_t key[8];
    uint64_t value = data(b);
    assert(value >= 0);
    loadKey(value, key);
    // fprintf(stderr, "add root %d\n", data(b));
    n_index++;
    ::insert(tree, &tree, key, 0, (uintptr_t) b, 8);
  }

  void insert_root_value(uintptr_t value) {
    uint8_t key[8];
    assert(value >= 0);
    assert(value < (1ULL << 40));
    // fprintf(stderr, "insert root value %llu\n", (unsigned long long) value);
    n_index++;
    loadKey(value, key);
    if (!::lookup(tree, key, 8, 0, 8)) {
      ::insert(tree, &tree, key, 0, (value << 1) | 1, 8);
    } else {
      // fprintf(stderr, "D"); // Duplicate.
    }
  }

  void load(T const *arr, int n) {
    LargeBucket<T> *root = new LargeBucket<T>();
    root->append(0); // Dummy leftmost bucket.
    assert(root);
    assert(!( ((uintptr_t) root) & 3));
    insert_root(root);

    int i = 0;
    while (i + LARGE_SIZE <= n) {
      LargeBucket<T> *b = new LargeBucket<T>();
      b->bulk_insert(arr + i, LARGE_SIZE);
      i += LARGE_SIZE;
      if (root) {
        ((LargeBucket<T>*) root)->add_chain(b);
      } else {
        root = b;
      }
    }

    while (i < n) {
      insert(arr[i++]);
    }
  }

  Node* find_bucket(uint64_t value64, Node **next) {
    uint8_t key[8];
    loadKey(value64, key);
    Node *ret = lower_bound_prev(tree, key, 8, 0, 8, next);
    assert(ret);
    assert(isLeaf(ret));
    return ret;
  }

  Node* find_bucket_upper(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    return ::lower_bound(tree,key,8,0,8);
  }

  void insert(T value) {
    // fprintf(stderr, "ins %d\n", value);
    Node *n = find_bucket(value, nullptr);
    assert(isLeaf(n));
    uintptr_t v = getData((uintptr_t) n);
    if (isPointer(v)) {
      Bucket<T> *b = (Bucket<T>*) v;
      if (b->large_type) {
        ((LargeBucket<T>*) v)->insert(value);
      } else if (!((SmallBucket<T>*) b)->insert(value)) {
        insert_root_value(value);
      }
    } else {
      insert_root_value(value);
    }
  }

  template <typename B>
  bool transition_to_art(B *b) {
    bool ok = erase_root(b->data(0));
    assert(ok);
    for (int i = 0; i < b->size(); i++) {
      // assert(b->data(i) != 888859321);
      // fprintf(stderr, "%d ", b->data(i));
      insert_root_value(b->data(i));
    }
    delete b;
    return (n_large + n_small == 0);
  }

  bool erase_root(uint64_t value64) {
    uint8_t key[8];
    loadKey(value64, key);
    // fprintf(stderr, "erase root %llu\n", value64);
    ::erase(tree,&tree,key,8,0,8);
    n_index--;
    return true;
  }

  LargeBucket<T>* make_standalone(LargeBucket<T> *lb, T value) {
    while (lb->next()) {
      LargeBucket<T> *rb = stochastic_split_chain(lb, rng);
      // fprintf(stderr, "stochastic_split_chain %d\n", rb->data(0));
      insert_root(rb);
      lb = (value < rb->data(0)) ? lb : rb;
    }
    return lb;
  }

  bool erase(T const &value) {
    // fprintf(stderr, "ERASE %d\n", value);
    // art_debug = 1;
    if (n_large + n_small) {
      Node *n = find_bucket(value, nullptr);
      assert(isLeaf(n));
      uintptr_t v = getData((uintptr_t) n);
      if (isPointer(v)) {
        Bucket<T> *b = (Bucket<T>*) v;
        if (b->large_type) {
          LargeBucket<T> *lb = make_standalone((LargeBucket<T>*) b, value);
          if (LARGE_TOUCH == 0) {
            transition_to_art(lb);
            return erase_root(value);
          } else if (lb->data(0) == value || lb->n_touched() + lb->n_erased() > LARGE_TOUCH) {
            b = to_small_bucket(lb, value);
          } else {
            return lb->erase(value, rng);
          }
        }
        SmallBucket<T> *sb = (SmallBucket<T>*) b;
        // fprintf(stderr, "smallize\n");
        if (!sb) return false;
        // fprintf(stderr, "trans\n");
        if (sb->data(0) == value || sb->size() < 10 || sb->n_touched() > SMALL_TOUCH) {
          // fprintf(stderr, ".");
          transition_to_art(sb);
        } else {
          int idx = sb->erase(value);
          assert(idx > 0);
          return true;
        }
      } else {
        // fprintf(stderr, "%llu %d\n", (unsigned long long) getData(v), value);
        assert(getData(v) == (uintptr_t) value);
      }
    }
    // fprintf(stderr, "err %d\n", value);
    return erase_root(value);
  }

  SmallBucket<T>* to_small_bucket(Bucket<T> *b, T value) {
    assert(b);
    if (b->large_type) {
      LargeBucket<T> *lb = (LargeBucket<T>*) b;
      erase_root(lb->data(0));
      SmallBucket<T> *target = nullptr;
      vector<Bucket<T>*> arr;
      lb->split(rng, [&](T *D, int n) {
        SmallBucket<T> *sb = new SmallBucket<T>(D, n);
        assert(sb->size());
        arr.push_back(sb);
        if (sb->data(0) <= value && (!target || sb->data(0) > target->data(0))) {
          target = sb;
        }
      });
      for (Bucket<T>* sb : arr) insert_root(sb);
      delete lb;
      b = target;
    }
    return (SmallBucket<T>*) b;
  }

  /* TODO: lazy lower_bound */
  T lower_bound(T const &value) {
    static int nth = 0; nth++;
    // art_debug = 1;
    // fprintf(stderr, "lower_bound %d\n", value);

    Node *next = nullptr;
    if (n_large + n_small) {
      Node *n = find_bucket(value, &next);
      assert(isLeaf(n));
      uintptr_t v = getData((uintptr_t) n);
      assert(v);
      if (isPointer(v)) {
        Bucket<T> *b = (Bucket<T>*) v;
        if (b->large_type) {
          LargeBucket<T> *lb = (LargeBucket<T>*) b;
          bool has_next = lb->next();
          lb = make_standalone(lb, value);
          int pos = lb->lower_pos(value, rng);
          if (pos < lb->size()) {
            T ret = lb->data(pos);
            if (LARGE_TOUCH == 0) transition_to_art(lb);
            else if (lb->n_touched() > LARGE_TOUCH) to_small_bucket(lb, value);
            return ret;
          }
          if (has_next) {
            // the root array may have changed.
            uint8_t key[8];
            loadKey(value, key);
            next = ::lower_bound(tree,key,8,0,8);
          }
        } else {
          SmallBucket<T> *sb = (SmallBucket<T>*) b;
          if (sb) {
            int pos = sb->lower_pos(value);
            if (pos < sb->size()) {
              T ret = sb->data(pos);
              if (sb->n_touched() > SMALL_TOUCH) transition_to_art(sb);
              return ret;
            }
          }
        }
      } else if ((T) getData(v) == value) {
        return value;
      }
    } else {
      uint8_t key[8];
      loadKey(value, key);
      next = ::lower_bound(tree,key,8,0,8);
    }

    return isLeaf(next) ? getLeafValue(next) : 0;
  }

  void statistics(std::function<void(int, int, long long, int, int, int, int, int, int, int, int, int, int, int)> cb) {
    long long n_bytes = 0;
    int N = -1, n_slack_art = 0, n_slack_leaves = 0, n_chain = 0, n_internal = 0,
        n_leaf = 0, art_n4 = 0, art_n16 = 0, art_n48 = 0, art_n256 = 0;

    art_visit(tree, [&](Node *n) {
      if (isLeaf(n)) {
        n_leaf++;
        uintptr_t v = getData((uintptr_t) n);
        assert(v);
        if (isPointer(v)) {
          Bucket<T> *b = (Bucket<T>*) v;
          if (b->large_type) {
            LargeBucket<T> *lb = (LargeBucket<T>*) b;
            n_slack_leaves += lb->slack();
            n_chain--;
            while (lb) {
              n_bytes += sizeof(LargeBucket<T>);
              N += lb->size();
              lb = lb->next();
              n_chain++;
            }
          } else {
            n_bytes += sizeof(SmallBucket<T>);
            SmallBucket<T> *sb = (SmallBucket<T>*) b;
            n_slack_leaves += sb->slack();
            N += sb->size();
          }
        } else {
          n_bytes += sizeof(uintptr_t);
          N++;
        }
      } else {
        n_internal++;
        switch (n->type) {
          case NodeType4: {
             Node4* node = static_cast<Node4*>(n);
             n_slack_art += 4 - node->count;
             art_n4++;
             n_bytes += sizeof(NodeType4);
             break;
          }
          case NodeType16: {
             Node16* node=static_cast<Node16*>(n);
             n_slack_art += 16 - node->count;
             art_n16++;
             n_bytes += sizeof(NodeType16);
             break;
          }
          case NodeType48: {
             Node48* node=static_cast<Node48*>(n);
             n_slack_art += 48 - node->count;
             art_n48++;
             n_bytes += sizeof(NodeType48);
             break;
          }
          case NodeType256: {
             Node256* node=static_cast<Node256*>(n);
             n_slack_art += 256 - node->count;
             art_n256++;
             n_bytes += sizeof(NodeType256);
             break;
          }
        }
      }
    });
    cb(N, n_index, n_bytes, n_slack_art, n_slack_leaves, n_internal, n_leaf, n_small, n_large, n_chain, art_n4, art_n16, art_n48, art_n256);
  }
};

#endif

/*
NOTEs:
- cracker_index can be any sorted data structure (e.g. ART, BTree, etc..)
- the transition to cracker_index only happens during delete operation. It's not beneficial to transit for insert (see experiment).

function insert(v) {
  ptr = search v in the cracker_index and get the pointer to leaf.
  if (ptr is a pointer to a bucket chain) {
    append v to bucket chain pointed by ptr.
  } else {
    inserts v to cracker_index or
    create a new bucket chain (in case of large inserts) and append.
  }
}

function erase(v) {
  ptr = search v in the cracker_index and get the pointer to leaf.
  if (ptr is a pointer to a bucket chain) {
    perform DDR to the bucket chain pointed by ptr until the target bucket has no chain.
    bucket = the first and only bucket pointed in the bucket chain.
    if (the bucket has been touched >= X times) {
      move all elements in the bucket to cracker_index.
      erase v from cracker_index.
    } else {
      erase v from the bucket.
    }
  } else {
    erase v from cracker_index.
  }
}

function query(v) {
  ptr = search v in the cracker_index and get the pointer to leaf.
  if (ptr is a pointer to a bucket chain) {
    perform DDR to the bucket chain pointed by ptr until the target bucket has no chain.
    bucket = the first and only bucket pointed in the bucket chain.
    pos = crack the bucket on v
    return pointer to (bucket, pos).
  } else {
    return pointer to ptr.
  }
}

NOTEs:
- the return value of query() can be a pointer to (cracked) bucket or pointer to cracker_index.
- we can traverse the next element in the (cracked) bucket or in the cracker_index.

*/
