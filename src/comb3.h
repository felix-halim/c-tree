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

#include <sys/mman.h>

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

#define CRACK_AT (cap >> 5)
#define DECRACK_AT (cap >> 6)
#ifndef INTERNAL_BSIZE
  #define INTERNAL_BSIZE 32
#endif

template <typename T, typename CMP  = std::less<T> >
static bool eq(T const &a, T const &b, CMP const &cmp){ return !cmp(a,b) && !cmp(b,a); }


template <typename T, typename CMP>
class Bucket {
 protected:
  bool leaf;
  int N;        // Number of data elements in this bucket pointed by D.
  Bucket *par;  // Pointer to the parent bucket (must be an InternalBucket)

 public:

  virtual ~Bucket() {};
  int size() const { return N; }
  bool is_leaf() const { return leaf; }
  Bucket* parent() const { return par; }

  void set_parent(Bucket* p) { par = p; }
};


template <typename T, typename CMP>
class InternalBucket : public Bucket<T, CMP> {
  T D[INTERNAL_BSIZE];            // Data values.
  Bucket<T, CMP>* C[INTERNAL_BSIZE + 1];  // Pointer to children (can be internals or leaves).

 public:

  InternalBucket(Bucket<T, CMP> *parent, Bucket<T, CMP> *left_child) {
    this->leaf = false;
    this->par = parent;
    this->N = 0;
    C[0] = left_child;
    left_child->set_parent(this);
  }

  void set_data(int i, T value) { assert(i >= 0 && i < this->N); D[i] = value; };

  int capacity() const { return INTERNAL_BSIZE; }
  Bucket<T, CMP>* child(int i) const { assert(i >= 0 && i <= this->N); return C[i]; }
  int child_pos(Bucket<T, CMP> *c) const {
    for (int i = 0; i <= this->N; i++) {
      if (C[i] == c) return i;
    }
    assert(0);
    return 0;
  }
  T& data(int i) { assert(i >= 0 && i < this->N); return D[i]; }
  const T* data_pointer(int i) const { assert(i >= 0 && i < this->N); return &D[i]; }
  T promote_last() { return D[--this->N]; }
  bool is_full() const { return this->N == INTERNAL_BSIZE; }
  virtual bool is_leaf() const { return false; }

  int lower_pos(T value, CMP &cmp) {
    int pos = 0;
    while (pos < this->N && cmp(D[pos], value)) pos++;
    return pos;
  }

  virtual void debug() const {
    fprintf(stderr, "N = %d  [ ", this->N);
    for (int i = 0; i < this->N; i++) fprintf(stderr, "%d ", D[i]);
    fprintf(stderr, "]\n");
  }

  void insert(T value, Bucket<T, CMP> *nb, int left, CMP &cmp) {
    assert(!this->is_full());
    // fprintf(stderr, "ii %d\n", b);
    int i = this->N - 1;
    while (i >= 0 && cmp(value, D[i])) {
      D[i + 1] = D[i];
      C[i + 2] = C[i + 1];
      i--;
    }
    D[i + 1] = value;
    if (left == -1) {
      C[i + 2] = C[i + 1];
      C[i + 1] = nb;
    } else {
      C[i + 2] = nb;
    }
    this->N++;
    nb->set_parent(this);
  }

  void erase(int pos, int stride) {
    this->N--;
    while (pos < this->N) {
      D[pos] = D[pos + 1];
      C[pos + stride] = C[pos + stride + 1];
      pos++;
    }
    if (!stride) C[pos] = C[pos + 1];
  }

  InternalBucket* middle_split() {
    int H = this->N / 2;
    InternalBucket* ib = new InternalBucket(this->par, C[H]);
    for (int i = H, j = 0; i < this->N; i++) {
      ib->D[j++] = D[i];
      ib->C[j] = C[i + 1];
      ib->C[j]->set_parent(ib);
    }
    ib->N = this->N - H;
    this->N = H;
    return ib;
  }
};


// Dynamically resize COMB bucket sizes.
// If number of cracks > 32, it splits to two smaller buckets.
template <typename T, typename CMP, int MAX_BSIZE>
class LeafBucket : public Bucket<T, CMP> {
  static const unsigned short MAX_CRACK = 64;

  int I;                // last indexed position
  unsigned char nC;     // the number of cracker indices
  unsigned long long S; // sorted bits
  int C[MAX_CRACK-1];   // the cracker indices
  T V[MAX_CRACK-1];     // the cracker value
  T *D;                 // the data elements
  int cap;              // the maximum number of elements in D.
  LeafBucket* next_b;   // buckets can be chained like a linked list of buckets
                        // the value of next is -1 if there is no next chain
                        // otherwise the index of the bucket [0, num_of_buckets)
  LeafBucket* tail_b;   // pointer to the last bucket in the chain.

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

  void flush_pending_inserts(CMP &cmp) {
    assert(I <= this->N);                     // Indexed index should be less than the number of elements
    if (!nC){ I = this->N; S = 0; return; }   // no index yet, all the elements are considered "inserted"
    assert(next_b == NULL);               // Indexes only makes sense when there is no chain

    // IMPROVE: bulk insert? (Currently using Merge Completely)
    int minC = nC;
    for (int j = I; I < this->N; j= ++I) {        // insert all pending elements (from I to N)
      int i = nC - 1;
      T tmp = D[j];            // store the pending tuple
      for (; i>=0 && cmp(tmp,V[i]); i--){  // insert by shuffling through the cracker indices C
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

  T* partition(T *F, T *L, T const &v, CMP &cmp) {
    while (true) {
      while (true)
        if (F == L) return F;
        else if (cmp(*F,v)) ++F;
        else break;
      --L;
      while (true)
        if (F == L) return F;
        else if (!bool(cmp(*L,v))) --L;
        else break;
      std::iter_swap(F, L);
      ++F;
    }
    assert(false);
  }

  // partitions roughly in the middle satisfying the DECRACK_AT
  T* rough_middle_partition(T *L, T *R, CMP &cmp, Random &rng) {
    assert(R-L >= CRACK_AT);
    int ntry = 10;
    for (T *i=L, *j=R, *p; ntry--; ){
      std::iter_swap(j-1, i + rng.nextInt(j-i));
      std::iter_swap(j-1, p = partition(i, j-1, *(j-1), cmp));
      if (p-L <= DECRACK_AT) i = p;
      else if (R-p <= DECRACK_AT) j = p;
      else return p;
    }
    T *M = L+((R-L)>>1);
    std::nth_element(L,M,R, cmp);
    return M;
  }

  // returns a piece [L,R) that contains v
  // it will reorganize the elements so that the R-L range gets smaller overtime
  int get_piece_by_value(T const &v, int &L, int &R, CMP &cmp, Random &rng) {
    flush_pending_inserts(cmp);
    int i = 0;
    while (i<nC && !cmp(v, V[i])) i++;      // find the cracker indices that covers v
    L = i==0? 0 : C[i-1];            // the left crack boundary
    R = i==nC? this->N : C[i];            // the right crack boundary
    while (R-L > CRACK_AT){            // narrow down the piece using DDR
      int M = rough_middle_partition(D+L+(i?1:0), D+R, cmp, rng) - D;
      add_cracker_index(i, M);
//        fprintf(stderr,"CRACKING %d %d, [%d %d]\n",M,D[M],L,R);
      if (cmp(v,D[M])) R=M; else L=M, i++;  // adjust the cracker index i
    }
    assert(i>=0 && i<=nC);
//      assert(check(D[0],false,D[0],false, cmp));
    return i;
  }

public:
  LeafBucket(Bucket<T, CMP> *p, int c): D(new T[c]), cap(c), next_b(NULL) {
    this->leaf = true;
    this->par = p;
    this->N = 0;
    clear_indexes();
  }

  ~LeafBucket() { delete[] D; }

  int n_cracks() const { return nC; }
  int slack() const { return cap - this->N; }
  void set_next(LeafBucket *b){ next_b = b; }
  void set_tail(LeafBucket *b){ tail_b = b; }
  void clear_indexes(){ S = nC = I = 0; }
  T randomValue(Random &rng) const { return D[rng.nextInt(this->N)]; }

  int capacity() const { return cap; }
  Bucket<T, CMP>* next() const { return next_b; }
  Bucket<T, CMP>* tail() const { return tail_b; }
  T& data(int i) { assert(i >= 0 && i < this->N); return D[i]; }
  const T* data_pointer(int i) const { assert(i >= 0 && i < this->N); return &D[i]; }
  bool is_full() const { return this->N == cap; }
  virtual bool is_leaf() const { return true; };

  void add_chain(LeafBucket *next) {
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
    // if (v == 1846371397) fprintf(stderr, "ins D[%d] = %d, to bucket = %p, cap = %d\n", N, v, this, cap);
    if (this->N < cap) {
      D[this->N++] = v;
      return;
    }
    assert(!tail_b || !tail_b->next());
    if (!tail_b || tail_b->N == cap) {
      add_chain(new LeafBucket(this->par, tail_b ? MAX_BSIZE : cap));
    }
    assert(tail_b && tail_b->N < tail_b->cap);
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
    assert(!is_full());
    D[this->N++] = value;
  }

  vector<pair<T, LeafBucket*>> split(CMP &cmp) {
    assert(cap >= 512);
    flush_pending_inserts(cmp);

    int nsplits = 16;
    assert(cap % nsplits == 0);
    cap /= nsplits;

    assert(nC > 1);
    assert(nC < 64 - nsplits);
    assert(this->N > cap);

    if (C[0] > cap) {
      assert(!piece_is_sorted(0));
      nth_element(D, D + cap, D + C[0]);
      // fprintf(stderr, "add cracker front %d\n", cap);
      add_cracker_index(0, cap);
      // debug("asdf", 0, 0);
      // assert(check(0, 0, 0, 0, cmp));
    }

    int newC = 0;
    T *DD = nullptr;
    vector<pair<T, LeafBucket*>> ret;
    for (int next_pos = cap, i = 0; ; next_pos += cap) {
      while (i + 1 < nC && C[i + 1] <= next_pos) i++;
      assert(C[i] <= next_pos);

      if (C[i] < next_pos) {
        assert(i + 1 == nC || C[i + 1] > next_pos);
        int R = ((i + 1) == nC) ? this->N : C[i + 1];
        if (next_pos < R) {
          piece_set_sorted(i, false);
          // assert(!piece_is_sorted(i));
          nth_element(D + C[i] + 1, D + next_pos, D + R);
          add_cracker_index(++i, next_pos);
        }
      }

      if (DD == nullptr) {
        DD = new T[cap];
        memcpy(DD, D, sizeof(T) * cap);
        newC = i;
      } else {
        int pos = next_pos - cap;
        int n_move = min(cap, this->N - pos);
        // fprintf(stderr, "n_move = %d, cap = %d, D[%d] = %d\n", n_move, cap, pos, D[pos]);
        assert(0 <= n_move && n_move <= cap);

        LeafBucket *b = new LeafBucket(this->par, cap);
        b->I = b->N = n_move;
        for (int j = 0; j < n_move; j++) b->D[j] = D[pos + j];

        assert(i < nC);
        int j = i + 1;
        for (; j < nC && C[j] < next_pos; j++) {
          b->C[b->nC] = C[j] - pos;
          b->V[b->nC] = V[j];
          b->piece_set_sorted(b->nC, piece_is_sorted(j));
          b->nC++;
        }
        // b->piece_set_sorted(b->nC, piece_is_sorted(j));
        // assert(b->check(0, 0, 0, 0, cmp));

        if (b->N < b->cap / 2) {
          // Chain it
          ret.back().second->set_next(b);
        } else {
          ret.push_back(make_pair(D[pos], b));
        }
      }
      if (next_pos >= this->N) break;
    }

    delete[] D;
    D = DD;
    nC = newC;
    this->N = I = cap;
    // if (!piece_is_sorted(0)) {
    //   std::sort(D, D + cap, cmp);
    //   piece_set_sorted(0, true);
    // }
    // debug("asdf", 0, 0);
    // b->debug("asdf", 0, 0);

    // assert(check(0, 0, 0, 0, cmp));

    return ret;
  }

  int bulk_insert(T const *v, int length) {
    assert(this->N == 0);
    memcpy(D, v, sizeof(T) * length);
    return this->N = length;
  }

  void mark_hi(T const &P, CMP &cmp, int *hi, int &nhi){
    for (int i=0; i < this->N; i++){
      hi[nhi] = i;
      nhi += !cmp(D[i], P);
    }
  }

  void mark_lo(T const &P, CMP &cmp, int *lo, int &nlo){
    for (int i=0; i < this->N; i++){
      lo[nlo] = i;
      nlo += cmp(D[i], P);
    }
  }

  // partition this bucket based on value v, destroying all cracker indices
  int partition(T const &v, CMP &cmp){
    clear_indexes();
    assert(this->N > 0 && this->N <= cap);
    return partition(D, D + this->N, v, cmp) - D;
  }

  void fusion(LeafBucket *that, int *hi, int *lo, int &nhi, int &nlo){
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    T *Lp = D, *Rp = that->D;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  virtual void debug() const {
    debug("debug", 0, 0, true);
  }

  bool debug(const char *msg, int i, int j, bool verbose) const {
    for (int k=0; k<nC; k++)
      fprintf(stderr,"C[%d/%d] = %d, %d (sorted = %d)\n",
        k,nC,C[k],(int)D[C[k]],piece_is_sorted(k));
    fprintf(stderr,"%s : i=%d/N=%d, j=%d/nC=%d, D[i,i+1] = %d, %d; I=%d, N=%d, next=%p\n",
      msg, i,this->N, j,nC, (int)D[i],(int)D[i+1], I,this->N,next_b);
    if (verbose) {
      fprintf(stderr, "[[");
      for (int k = 0; k < this->N; k++) {
        fprintf(stderr, "%d ", D[k]);
      }
      fprintf(stderr, "]]\n");
    }
    return false;
  }

  // call this function to check the consistency of this LeafBucket structure
  bool check(T lo, bool useLo, T hi, bool useHi, CMP &cmp) const {
    if (useLo) for (int i=0; i < this->N; i++) if (cmp(D[i],lo)){
      fprintf(stderr,"D[%d] = %d, lo = %d, bucket = %p\n",i, D[i], lo, this);
      return debug("useLo failed", i,0);
    }
    if (useHi) for (int i=0; i < this->N; i++) if (!cmp(D[i],hi)){
      fprintf(stderr,"D[%d] = %d, hi = %d, bucket = %p\n",i, D[i], hi, this);
      return debug("useHi failed", i,0);
    }
    for (int i = 0; i < nC; i++) {
      if (!eq(V[i], D[C[i]], cmp)) {
        fprintf(stderr, "not equal %d, %d != %d\n", i, V[i], D[C[i]]);
        return false;
      }
    }
    for (int i=0,j=0; i<I; i++){                    // check cracker indices
      if (j<nC && C[j]==i) assert(eq(V[j],D[i],cmp)), lo = D[i], j++;
      if (piece_is_sorted(j) && (j==nC? (i+1<I) : (i<C[j])) && cmp(D[i+1], D[i]))
        return debug("sortedness violation", i,j);
      if (j>0 && cmp(D[i], D[C[j-1]])) return debug("lower bound fail", i,j);
      if (j<nC && cmp(D[C[j]], D[i])) return debug("upper bound fail", i,j);
      // if (j<nC && !(cmp(D[i], D[C[j]]))) return debug("upper bound fail", i,j);
    }
    return true;
  }

  int index_of(T const &v, CMP &cmp) const {
    for (int i=0; i < this->N; i++) if (eq(D[i],v,cmp)) return i;
    return -1;
  }

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified LeafBucket "to", destroying all cracker indices
  void moveToFromIdx(LeafBucket *to, int fromIdx) {
    clear_indexes(); to->clear_indexes();    // destroy both buckets' cracker indices
    assert(this->N > fromIdx);            // make sure there is something to move
    assert(to->N + this->N - fromIdx <= cap);    // make sure the receiver has enough space
    memmove(to->D + to->N, D+fromIdx, (this->N - fromIdx) * sizeof(T));
    to->N += this->N - fromIdx;
    this->N = fromIdx;
  }

  int lower_pos(T value, CMP &cmp, Random &rng) {
    int i, L, R;
    return crack(value, i, L, R, true, cmp, rng);
  }

  int crack(T const &v, int &i, int &L, int &R, bool sort_piece, CMP &cmp, Random &rng){
    assert(!next());            // it doesn't make sense crack a chained bucket!
    i = get_piece_by_value(v,L,R,cmp,rng);    // find the piece [L,R) containing v
    assert(L>=0 && L<=R && R <= this->N);        // range check
    if (!piece_is_sorted(i)){
      if (sort_piece){            // sort the piece if requested
        std::sort(D+L,D+R,cmp);
        piece_set_sorted(i,true);
      } else {
        for (int at=L; at<R; at++)
          if (eq(D[at],v,cmp)) return at;
        return R;
      }
    }
    for (int pos = L; pos < R; pos++) if (!cmp(D[pos], v)) return pos;
    return R;
    // return std::lower_bound(D+L, D+R, v, cmp) - D;    // find the element v using binary search
  }

  T remove_largest(CMP &cmp, Random &rng) {
    flush_pending_inserts(cmp);
    while (true) {
      int i = nC;
      int L = (i == 0) ? 0 : C[i - 1];            // the left crack boundary
      int R = (i ==nC) ? this->N : C[i];            // the right crack boundary
      assert(L < R);
      if (R - L <= CRACK_AT) {
        if (nC && R-L+1 <= DECRACK_AT) {
          remove_cracker_index((i>0)?(--i):i);
        }
        int j = L;
        for (i = L + 1; i < R; i++)
          if (cmp(D[j], D[i])) j = i;
        swap(D[j], D[R - 1]);
        // fprintf(stderr, "REM LARGEST I = %d, N = %d, R = %d, largest = %d\n", I, this->N, R, D[R - 1]);
        assert(this->N == R);
        return D[I = this->N = R - 1];
      }
      get_piece_by_value(D[L + rng.nextInt(R - L)], L, R, cmp, rng);
    }
  }

  bool erase(T const &v, CMP &cmp, Random &rng) {
    assert(!next_b);
    int i, L, R, at = crack(v, i, L, R, false, cmp, rng);
    if (at >= R || !eq(D[at],v,cmp)) return false;  // the element to be erased is not found!

    // decrack this cracker piece (it becomes too small) or
    // if the deleted element index is a cracker index
    if (nC && R-L+1 <= DECRACK_AT){
      remove_cracker_index((i>0)?(--i):i);
    } else if (i > 0 && at == L){
      at = L+1;
      for (int j=L+2; j<R; j++)    // find a replacement element for the cracker index
        if (cmp(D[j], D[at])) at = j;  // that is the smallest element in the piece (L,R)
      std::swap(D[L], D[at]);
      piece_set_sorted(i,false);
      V[i-1] = D[L];
    }

    assert(at<R && eq(D[at],v,cmp));    // the element v must be found!

    // IMPROVE: use pending delete? antimatter?
    piece_set_unsorted_onwards(i);      // unset the sorted bit i onwards
    assert(i==nC || cmp(D[at],V[i]));
    for (int j=i; j<nC; j++){        // shuffle out the deleted element
      R = C[j]--;
      D[at] = D[R-1];
      D[R-1] = D[R];
      at = R;
    }
    D[at] = D[--this->N];    // the deleted element has been shuffled out from the bucket
    I--;        // adjust the pending index

//      assert(check(D[0],false,D[0],false,cmp));
    return true;
  }
};


template <typename T, typename CMP  = std::less<T>, int MAX_BSIZE = 8192>
class Comb {
  typedef LeafBucket<T, CMP, MAX_BSIZE> leaf_bucket_t;

  CMP cmp;     // The comparator functor.
  Random rng;  // The random number generator.
  Bucket<T, CMP> *root;

public:
  class iterator {
    void advance() {
      assert(bucket);
      assert(idx < bucket->size());
      if (++idx == bucket->size()) {
        if (bucket->is_leaf() && ((leaf_bucket_t*) bucket)->next()) {
          bucket = ((leaf_bucket_t*) bucket)->next();
          assert(bucket->size());
          idx = 0;
        } else {
          while (true) {
            InternalBucket<T, CMP> *parent = (InternalBucket<T, CMP>*) bucket->parent();
            if (!parent) { bucket = nullptr; idx = 0; break; }
            idx = parent->child_pos(bucket);
            bucket = parent;
            if (idx < parent->size()) break;
          }
        }
      }
    }

   public:
    typedef Comb<T, CMP, MAX_BSIZE> crack_type;
    crack_type *crack;
    Bucket<T, CMP> *bucket;
    int idx;

    // STL specific typedefs.
    typedef iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef int difference_type;
    typedef std::forward_iterator_tag iterator_category;

    iterator(crack_type *ct, Bucket<T, CMP> *b, int i): crack(ct), bucket(b), idx(i) {}

    self_type operator++() {
      self_type i = *this;
      advance();
      return i;
    }

    self_type operator++(int junk) {
      advance();
      return *this;
    }

    const reference operator*() {
      assert(bucket);
      assert(idx < bucket->size());
      return bucket->is_leaf()
        ? ((leaf_bucket_t*) bucket)->data(idx)
        : ((InternalBucket<T, CMP>*) bucket)->data(idx);
    }

    const pointer operator->() {
      assert(bucket);
      assert(idx < bucket->size());
      return bucket->data_pointer(idx);
    }

    bool operator==(const self_type &that) {
      return crack == that.crack && bucket == that.bucket && idx == that.idx;
    }

    bool operator!=(const self_type &that) {
      return !(*this == that);
    }

    // bool debug() {
    //   fprintf(stderr,"iter: v=%d, ridx=%d/%d, bidx=%d/%d, idx=%d/%d\n",
    //     crack->get_bucket_element(bidx,idx),
    //     ridx,crack->root_size(),
    //     bidx, crack->num_of_buckets(),
    //     idx, crack->bucket_size(bidx));
    //   return false;
    // }
  };

  Comb() {
    root = new leaf_bucket_t(nullptr, MAX_BSIZE);
  }
  
  void load(T const *arr, int n) {
    // fprintf(stderr, "batch %d\n", n);
    int i = 0;
    if (MAX_BSIZE <= n) {
      ((leaf_bucket_t*) root)->bulk_insert(arr, MAX_BSIZE);
      i = MAX_BSIZE;
    }
    while (i + MAX_BSIZE <= n) {
      leaf_bucket_t *b = new leaf_bucket_t(nullptr, MAX_BSIZE);
      b->bulk_insert(arr + i, MAX_BSIZE);
      i += MAX_BSIZE;
      // fprintf(stderr, "chain %d, %d\n", i, b->size());
      ((leaf_bucket_t*) root)->add_chain(b);
    }
    // fprintf(stderr, "done %d %d\n", i, n);
    while (i < n) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "inserted %d elements\n", size());
  }

  int size(Bucket<T, CMP> *b = nullptr) {
    if (!b) b = root;
    int ret = 0;
    if (b->is_leaf()) {
      while (b) {
        ret += b->size();
        b = ((leaf_bucket_t*) b)->next();
      }
      return ret;
    }
    for (int i = 0; i <= b->size(); i++) {
      ret += size(((InternalBucket<T, CMP>*) b)->child(i));
    }
    return ret;
  }

  void insert(T const &value) {
    // fprintf(stderr, "ins %d\n", value);
    assert(root);
    for (Bucket<T, CMP> *b = root; ; ) {
      if (b->is_leaf()) {
        ((leaf_bucket_t*) b)->insert(value);
        break;
      }
      int pos = ((InternalBucket<T, CMP>*) b)->lower_pos(value, cmp);
      b = ((InternalBucket<T, CMP>*) b)->child(pos);
    }
  }

  // Returns <bucket, pos> if found in internal node, otherwise returns <bucket, splitted> for leaf node.
  pair<Bucket<T, CMP>*, int> find_bucket(T value, bool include_internal) {
    Bucket<T, CMP> *b = root;
    int splitted = 0;
    // fprintf(stderr, "find_bucket %p\n", b);
    while (true) {
      if (b->is_leaf()) {
        // fprintf(stderr, "find_bucket2 %p\n", b);
        if (!split_chain((leaf_bucket_t*) b)) {
          return make_pair(b, splitted);
        }
        // fprintf(stderr, "find_bucket3 %p\n", b);
        splitted = 1;
        b = b->parent();
        assert(b);
      } else {
        int pos = ((InternalBucket<T, CMP>*) b)->lower_pos(value, cmp);
        if (include_internal && pos < b->size() && eq(((InternalBucket<T, CMP>*) b)->data(pos), value, cmp)) {
          // fprintf(stderr, "find_bucket4 %p\n", b);
          return make_pair(b, pos); // Found in the internal bucket.
        }
        b = ((InternalBucket<T, CMP>*) b)->child(pos);    // Search the child.
      }
    }
  }

  bool split_bucket(leaf_bucket_t *b) {
    // To avoid having too many cracker indexes in a bucket.
    if (b->n_cracks() > 40 && b->capacity() == MAX_BSIZE) {
      // assert(check());
      // fprintf(stderr, ".%lu", R.size());

      vector<pair<T, leaf_bucket_t*>> nb = b->split(cmp);
      for (int i = 1; i < (int) nb.size(); i++) {
        assert(cmp(nb[i-1].first, nb[i].first));
      }
      for (auto it = nb.rbegin(); it != nb.rend(); it++) {
        insert_internal(b, *it);
      }
      // assert(check());
      return true;
    }
    return false;
  }

  bool erase(T const &value) {
    // assert(check());
    // fprintf(stderr, "ERASED %d\n", value);
    // debug();

    pair<Bucket<T, CMP>*, int> p = find_bucket(value, true);
    if (p.first->is_leaf()) {
      // fprintf(stderr, "ERASE1 %d\n", value);
      // fprintf(stderr, "a");
      bool ret = ((leaf_bucket_t*) p.first)->erase(value, cmp, rng);
      // fprintf(stderr, "b");
      // split_bucket((leaf_bucket_t*) p.first);
      return ret;
    }

    // Found in an internal node, delete the largest node <= value.
    pair<Bucket<T, CMP>*, int> upper = find_bucket(value, false);

    // It is possible that finding upper invalidated p's references.
    if (upper.second) p = find_bucket(value, true); // Refresh.

    // fprintf(stderr, "ERASE2 %d\n", value);

    // upper.first is the leaf bucket containing the value.
    // upper.second signify whether a leaf_split happened when finding the bucket.
    Bucket<T, CMP> *b = upper.first;
    T next_largest = 0;
    assert(b->is_leaf());
    if (b->size()) {
      // fprintf(stderr, "ERASE3 %d\n", value);
      next_largest = ((leaf_bucket_t*) b)->remove_largest(cmp, rng);
    } else {
      // Bucket b is empty, search ancestors.
      while (true) {
        // fprintf(stderr, "ERASE4 %d\n", b);
        assert(b != p.first);
        InternalBucket<T, CMP> *parent = (InternalBucket<T, CMP>*) b->parent();
        assert(parent);
        delete b;
        if (parent->size()) {
          if (parent == p.first) {
            // Found in the same internal node as p.first, just delete it.
            ((InternalBucket<T, CMP>*) p.first)->erase(p.second, 0);
            // fprintf(stderr, "yay\n");
            return true;
          }
          next_largest = parent->promote_last();
          break;
        }
        b = parent;
      }
    }
    // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
    assert(!p.first->is_leaf());
    assert(((InternalBucket<T, CMP> *) p.first)->data(p.second) == value);
    assert(next_largest <= value);
    ((InternalBucket<T, CMP> *) p.first)->set_data(p.second, next_largest);
    // assert(check());
    // fprintf(stderr, "ERASE5 %d\n", value);
    return true;
  }

  // add a LeafBucket (bidx) to the root chain 'ridx'
  bool add_to_chain(leaf_bucket_t *&chain, leaf_bucket_t *b) {
    if (!chain) { // the root chain is empty.
      chain = b;  // b is the head of the chain.
      chain->set_next(nullptr);
      chain->set_tail(nullptr);
    } else {
      assert(!chain->tail() || !((leaf_bucket_t*) chain->tail())->next());
      if (chain->slack()) {
        b->moveToFromIdx(chain, b->size() - std::min(b->size(), chain->slack()));
      } else if (chain->tail()) {
        if (((leaf_bucket_t*) chain->tail())->slack()) {
          b->moveToFromIdx((leaf_bucket_t*) chain->tail(), b->size() - std::min(b->size(), ((leaf_bucket_t*) chain->tail())->slack()));
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

  T get_random_pivot(leaf_bucket_t *B, CMP &cmp, Random &rng) {  // pick the pivot near the median
    // leaf_bucket_t *C = (leaf_bucket_t*) B->tail();

    assert(B->size() > 11);

    // fprintf(stderr, "Picking random 11 elements, sizes = %d + %d\n", B->size(), T->size());
    // while (B->size() < 11) B->append(C->leaf_promote_last());

    T R[11]; // Randomly pick 11 elements from B.
    for (int i = 0; i < 11; i++) R[i] = B->remove_random_data(rng);

    // Replace R with the next buckets in the chain using reservoir sampling.
    leaf_bucket_t *Nb = (leaf_bucket_t*) B->next();
    for (int i = 1; Nb; i++) {
      int j = rng.nextInt(i);
      if (j < 11) Nb->swap_random_data_with(R[j], rng);
      Nb = (leaf_bucket_t*) Nb->next();
    }

    std::nth_element(R, R + 5, R + 11);
    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      if (i != 5) B->append(R[i]);
    }
    // fprintf(stderr, "pivot = %d\n", pivot);
    return R[5];
  }

  // fusion
  pair<T, leaf_bucket_t*> stochastic_split_chain(leaf_bucket_t *b, CMP &cmp, Random &rng) {
    const T &p = get_random_pivot(b, cmp, rng);

    // fprintf(stderr, "stochastic_split_chain %p\n", b);
    leaf_bucket_t *left_chain = nullptr;
    leaf_bucket_t *right_chain = nullptr;
    leaf_bucket_t *Lb = nullptr, *Rb = nullptr;
    int hi[MAX_BSIZE], lo[MAX_BSIZE];
    int nhi = 0, nlo = 0;

    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        Lb->fusion(Rb, hi, lo, nhi, nlo);
        if (!nhi) { add_to_chain(left_chain, Lb); Lb = nullptr; }
        if (!nlo) { add_to_chain(right_chain, Rb); Rb = nullptr; }
      } else if (!Lb) {
        if (!b) break;
        Lb = b;
        b = (leaf_bucket_t*) b->next();
      } else if (!nhi) {
        assert(Lb);
        Lb->mark_hi(p, cmp, hi, nhi);
        if (!nhi) { add_to_chain(left_chain, Lb); Lb = nullptr; }
      } else if (!Rb) {
        if (!b) break;
        Rb = b;
        b = (leaf_bucket_t*) b->next();
      } else if (!nlo) {
        assert(Rb);
        Rb->mark_lo(p, cmp, lo, nlo);
        if (!nlo) { add_to_chain(right_chain, Rb); Rb = nullptr; }
      } else {
        assert(0);
      }
    }

    if (Rb) { assert(!Lb); Lb = Rb; }
    if (Lb) {
      if (Lb->size()) {
        int i = Lb->partition(p, cmp);
        if (i == 0) {
          add_to_chain(right_chain, Lb);
        } else if (i == Lb->size()) {
          add_to_chain(left_chain, Lb);
        } else {
          Rb = new leaf_bucket_t(Lb->parent(), Lb->capacity());
          Lb->moveToFromIdx(Rb, i);
          add_to_chain(left_chain, Lb);
          add_to_chain(right_chain, Rb);
        }
      } else {
        delete Lb;
      }
    }

//    assert(check());
    return make_pair(p, right_chain);
  }

  void insert_internal(leaf_bucket_t *leafb, pair<T, Bucket<T, CMP>*> right_chain) {
    InternalBucket<T, CMP> *parent = (InternalBucket<T, CMP>*) leafb->parent();
    // fprintf(stderr, "parent = %p, right_chain = %p\n", parent, right_chain.second);
    while (parent && right_chain.second) {
      if (parent->is_full()) {
        // fprintf(stderr, "parful\n");
        // Optional optimization: transfer_one_to_left_or_right();
        InternalBucket<T, CMP> *inb = parent->middle_split();
        T midValue = parent->promote_last();
        if (!cmp(right_chain.first, midValue)) {
          inb->insert(right_chain.first, right_chain.second, 0, cmp);
        } else {
          parent->insert(right_chain.first, right_chain.second, 0, cmp);
        }
        right_chain.first = midValue;
        right_chain.second = inb;
        parent = (InternalBucket<T, CMP>*) right_chain.second->parent();
      } else {
        // fprintf(stderr, "internal\n");
        parent->insert(right_chain.first, right_chain.second, 0, cmp);
        right_chain.second = nullptr;
      }
    }
    // assert(check());
    // fprintf(stderr, "split_chain3 = %p\n", leafb);
    if (right_chain.second) {
      // fprintf(stderr, "OLD ROOT %p\n", root);
      assert(!parent);
      assert(!root->parent());
      root = new InternalBucket<T, CMP>(nullptr, root);
      ((InternalBucket<T, CMP>*) root)->insert(right_chain.first, right_chain.second, 0, cmp);
      // fprintf(stderr, "NEW ROOT %d\n", size(root));
    }
  }

  bool split_chain(leaf_bucket_t *leafb) {
    // assert(check());
    if (!leafb->next()) return false;

    // fprintf(stderr, "split_chain %p\n", leafb);
    // fprintf(stderr, "split_chain %d\n", size(leafb));
    pair<T, leaf_bucket_t*> right_chain = stochastic_split_chain(leafb, cmp, rng);
    // fprintf(stderr, "split_chain2 %d %d\n", size(leafb), size(right_chain.second));
    // assert(check());
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);

    insert_internal(leafb, right_chain);

    // fprintf(stderr, "done split %p\n\n\n", leafb);
    // debug();
    // assert(check());
    return true;
  }

  iterator begin() {
    Bucket<T, CMP> *bucket = root;
    while (!bucket->is_leaf()) {
      bucket = ((InternalBucket<T, CMP>*) bucket)->child(0);
    }
    return iterator(this, bucket, 0);
  }

  iterator end() {
    return iterator(this, nullptr, 0);
  }

  /* TODO: lazy lower_bound */
  iterator lower_bound(T const &value, bool sort_piece = true) {
    static int nth = 0; nth++;
    // assert(check());
    // TODO: optimize leaf slack

    // fprintf(stderr, "lower_bound %d, %d\n", value, size());

    pair<Bucket<T, CMP>*, int> p = find_bucket(value, true);

    if (!p.first->is_leaf()) {
      // Found in internal bucket.
      // fprintf(stderr, "lower_bound2 %d\n", value);
      return iterator(this, p.first, p.second);
    }

    int pos = ((leaf_bucket_t*) p.first)->lower_pos(value, cmp, rng);
    if (pos < p.first->size()) {
      // fprintf(stderr, "lower_bound3 %d, %d, bucket = %p\n", value, ((leaf_bucket_t*) p.first)->data(pos), p.first);
      // ((leaf_bucket_t*) p.first)->debug();
      return iterator(this, p.first, pos);
    }

    InternalBucket<T, CMP> *ib = (InternalBucket<T, CMP>*) p.first->parent();
    while (ib) {
      // fprintf(stderr, "lower_bound4 %d\n", value);
      pos = ib->lower_pos(value, cmp);
      if (pos < ib->size()) {
        return iterator(this, ib, pos);
      }
      ib = (InternalBucket<T, CMP>*) ib->parent();
    }

    return end();
  }

  // erase from [v1, v2)
  void erase(T const &v1, T const &v2) {
    assert(!cmp(v2, v1));
    break_chain(find_root(v2), v2);
    break_chain(find_root(v1), v1);
    // TODO: destroy B from i1 to i2-1
  }

  int count(T const &v1, T const &v2){
    assert(!cmp(v2,v1));
    iterator it1 = lower_bound(v1);
    iterator it2 = lower_bound(v2);
    return count(it1, it2);
  }

  int count(iterator it1, iterator it2){
    // assert(!it1.bucket->next());
    // assert(!it2.bucket->next());
    if (it1.bucket == it2.bucket) return it2.idx - it1.idx;
    int ret = it1.bucket->size() - it1.idx;
    // for (it1.root_iter++; it1.root_iter != it2.root_iter; it1.root_iter++) {
    //   ret += size(root_value(it1.root_iter));
    // }
    return ret + it2.idx;
  }

  bool check() {
    for (iterator it = begin(), e = end(); it != e; it++) {
      /*
      if (b->next()) {
        if (!b->tail()) {
          fprintf(stderr, "Has next but no tail\n");
          return false;
        }
        if (b->next() == b->tail()) {
          fprintf(stderr, "Next and tail must point to different bucket\n");
          return false;
        }
      } else {
        if (b->tail()) {
          fprintf(stderr, "No next but has tail\n");
          return false;
        }
      }

      if (it != R.end()) upper = *it;

      for (; b; b = b->next()) {
        if (!b->check(lower, i > 0, upper, it != R.end(), cmp)) {
  //          for (int j=0; j<root_size(); j++)
  //            fprintf(stderr,"root[%d] = %d\n",j,R[j]);
          fprintf(stderr,"Fail ridx = %d / %lu\n", i, R.size());
          return false;
        }
      }
      */
    }
    return true;
  }

  int exists(T const &v, Bucket<T, CMP> *b = nullptr, bool print = false) {
    if (!b) b = root;
    int ret = 0;
    if (b->is_leaf()) {
      leaf_bucket_t *lb = (leaf_bucket_t*) b;
      while (lb) {
        ret += lb->size();
        lb = (leaf_bucket_t*) lb->next();
      }
    } else {
      InternalBucket<T, CMP> *ib = (InternalBucket<T, CMP>*) b;
      for (int j = 0; j < ib->size(); j++) {
        if (ib->data(j) == v) {
          if (print) fprintf(stderr,"v=%d, Exists at ridx = %p, at=%d/%d\n", v, ib, j, b->size());
        }
        ret += exists(v, ib->child(j), print);
      }
      ret += exists(v, ib->child(ib->size()), print);
      ret += ib->size();
    }
    return ret;
  }
};

/*

  void eager_insert(T const &v, bool use_pos = true){
    int ridx = find_ridx(v);
    assert(Pb[ridx]!=-1);
    if (ridx == 0 && (B[Pb[0]].size() == 0 || cmp(v, R[0]))) R[0] = v;
    assert(Pe[ridx]!=-1);
    int cnt = 0;
    while (!B[Pe[ridx]].free()){
      split_chain(ridx);
      if (!cmp(v, R[ridx+1])) ridx++;  // readjust root idx
      cnt++;
    }
    B[Pe[ridx]].insert(v);
  }

  bool nth(int idx, T &res){
    static_assert(USE_POS, "nth must set USE_POS");
    int lo=1, hi=root_size()-1, ridx=0;
    while (lo<=hi){
      int mid = lo + ((hi-lo)>>1);
      if (bit_sum(mid-1) > idx){
        hi = mid-1;
      } else {
        ridx = mid;
        lo = mid+1;
      }
    }
    if (B[Pb[ridx]].next() != -1){
      split_chain(ridx);  // randomly split the chain into two
      bit_init();      // TODO: improve!
      return nth(idx, res);
    }
    int bnth = idx - bit_sum(ridx-1);
    if (bnth >= B[Pb[ridx]].size()) return false;
    res = B[Pb[ridx]].nth(bnth, cmp, rng);
    return true;
  }

  int rank(T const &v){
    assert(USE_POS);
    int ridx = find_ridx(v);
    if (ridx == 0 && B[Pb[ridx]].size()==0) return 0;
    //int cnt = 
    break_chain(ridx,v);
    int i,L,R,idx = B[Pb[ridx]].crack(v,i,L,R,true,cmp,rng);
    if (ridx == 0) return idx;
    return bit_sum(ridx - 1) + idx;
  }

  int get_piece_by_index(int idx, int &L, int &R, CMP &cmp, Random &rng) {
    assert(idx>=0 && idx<N);        // index must fall within the bucket's size
    flush_pending_inserts(cmp);
    int i = 0;
    while (i<nC && idx >= C[i]) i++;    // find the cracker indices that covers the idx
    L = i==0? 0 : C[i-1];          // the left crack
    R = i==nC? N : C[i];          // the right crack
    while (R-L > CRACK_AT){          // narrow down the piece using DDR
      int M = rough_middle_partition(D+L+(i?1:0), D+R, cmp, rng) - D;
      add_cracker_index(i,M);
      if (idx < M) R=M; else L=M, i++;  // adjust the cracker index
    }
    assert(i>=0 && i<=nC);
//      assert(check(D[0],false,D[0],false, cmp));
    return i;
  }

*/

#endif
