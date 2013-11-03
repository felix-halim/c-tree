#ifndef _COMB_H_
#define _COMB_H_

#include <stdio.h>
#include <string.h>

#include <functional>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

using namespace std;

#include <sys/mman.h>

#include "google/btree_set.h"
#include "google/btree_map.h"

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

template <typename T, typename CMP  = std::less<T> >
static bool eq(T const &a, T const &b, CMP const &cmp){ return !cmp(a,b) && !cmp(b,a); }

// Dynamically resize COMB bucket sizes.
// If number of cracks > 32, it splits to two smaller buckets.
template <typename T, typename CMP  = std::less<T>>
class Bucket {
  static const unsigned short MAX_CRACK = 64;

  int N;                // the number of data elements
  int I;                // last indexed position
  unsigned char nC;     // the number of cracker indices
  unsigned long long S;       // sorted bits
  int C[MAX_CRACK-1];   // the cracker indices
  T V[MAX_CRACK-1];     // the cracker value
  T *D;                 // the data elements
  int cap;              // the maximum number of elements in D.
  Bucket* next_b;       // buckets can be chained like a linked list of buckets
                        // the value of next is -1 if there is no next chain
                        // otherwise the index of the bucket [0, num_of_buckets)

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
    assert(I <= N);                     // Indexed index should be less than the number of elements
    if (!nC){ I = N; S = 0; return; }   // no index yet, all the elements are considered "inserted"
    assert(next_b == NULL);               // Indexes only makes sense when there is no chain

    // IMPROVE: bulk insert? (Currently using Merge Completely)
    int minC = nC;
    for (int j = I; I < N; j= ++I) {        // insert all pending elements (from I to N)
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
    R = i==nC? N : C[i];            // the right crack boundary
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
  Bucket(int c): N(0), D(new T[c]), cap(c), next_b(NULL) { clear_indexes(); }
  ~Bucket() { delete[] D; }

  int size() const { return N; }
  int slack() const { return cap - N; }
  int is_full() const { return N == cap; }
  int capacity() const { return cap; }
  int n_cracks() const { return nC; }
  bool empty() { return N == 0; }
  Bucket* next() const { return next_b; }
  void set_next(Bucket *b){ next_b = b; }
  void clear_indexes(){ S = nC = I = 0; }
  T randomValue(Random &rng) const { return D[rng.nextInt(N)]; }
  T data(int i) const { assert(i >= 0 && i < N); return D[i]; }
  T* getp(int i) { assert(i>=0 && i<N); return &D[i]; }
  void insert(T const &v) {
    // if (v == 1846371397) fprintf(stderr, "ins D[%d] = %d, to bucket = %p, cap = %d\n", N, v, this, cap);
    assert(N < cap); D[N++] = v; }

  vector<pair<T, Bucket*>> split(CMP &cmp) {
    assert(cap >= 512);
    flush_pending_inserts(cmp);

    int nsplits = 16;
    assert(cap % nsplits == 0);
    cap /= nsplits;

    assert(nC > 1);
    assert(nC < 64 - nsplits);
    assert(N > cap);

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
    vector<pair<T, Bucket*>> ret;
    for (int next_pos = cap, i = 0; next_pos <= N; next_pos += cap) {
      while (i + 1 < nC && C[i + 1] <= next_pos) i++;
      assert(C[i] <= next_pos);

      if (C[i] < next_pos) {
        assert(i + 1 == nC || C[i + 1] > next_pos);
        int R = ((i + 1) == nC) ? N : C[i + 1];
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
        int n_move = min(cap, N - pos);
        // fprintf(stderr, "n_move = %d, cap = %d, D[%d] = %d\n", n_move, cap, pos, D[pos]);
        assert(0 <= n_move && n_move <= cap);

        Bucket *b = new Bucket(cap);
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
    }

    delete[] D;
    D = DD;
    nC = newC;
    N = I = cap;
    // if (!piece_is_sorted(0)) {
    //   std::sort(D, D + cap, cmp);
    //   piece_set_sorted(0, true);
    // }
    // debug("asdf", 0, 0);
    // b->debug("asdf", 0, 0);

    // assert(check(0, 0, 0, 0, cmp));

    return ret;
  }

  int bulk_insert(T const *v, int length){
    assert(N == 0);
    memcpy(D, v, sizeof(T) * length);
    return N = length;
  }

  void mark_hi(T const &P, CMP &cmp, int *hi, int &nhi){
    for (int i=0; i<N; i++){
      hi[nhi] = i;
      nhi += !cmp(D[i], P);
    }
  }

  void mark_lo(T const &P, CMP &cmp, int *lo, int &nlo){
    for (int i=0; i<N; i++){
      lo[nlo] = i;
      nlo += cmp(D[i], P);
    }
  }

  // partition this bucket based on value v, destroying all cracker indices
  int partition(T const &v, CMP &cmp){
    clear_indexes();
    assert(N > 0 && N <= cap);
    return partition(D, D + N, v, cmp) - D;
  }

  void fusion(Bucket *that, int *hi, int *lo, int &nhi, int &nlo){
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    T *Lp = D, *Rp = that->D;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  bool debug(const char *msg, int i, int j) const {
    for (int k=0; k<nC; k++)
      fprintf(stderr,"C[%d/%d] = %d, %d (sorted = %d)\n",
        k,nC,C[k],(int)D[C[k]],piece_is_sorted(k));
    fprintf(stderr,"%s : i=%d/N=%d, j=%d/nC=%d, D[i,i+1] = %d, %d; I=%d, N=%d, next=%p\n",
      msg, i,N, j,nC, (int)D[i],(int)D[i+1], I,N,next_b);
    return false;
  }

  // call this function to check the consistency of this Bucket structure
  bool check(T lo, bool useLo, T hi, bool useHi, CMP &cmp) const {
    if (useLo) for (int i=0; i<N; i++) if (cmp(D[i],lo)){
      fprintf(stderr,"D[%d] = %d, lo = %d, bucket = %p\n",i, D[i], lo, this);
      return debug("useLo failed", i,0);
    }
    if (useHi) for (int i=0; i<N; i++) if (!cmp(D[i],hi)){
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
    for (int i=0; i<N; i++) if (eq(D[i],v,cmp)) return i;
    return -1;
  }

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified Bucket "to", destroying all cracker indices
  void moveToFromIdx(Bucket *to, int fromIdx){
    clear_indexes(); to->clear_indexes();    // destroy both buckets' cracker indices
    assert(N > fromIdx);            // make sure there is something to move
    assert(to->N + N-fromIdx <= cap);    // make sure the receiver has enough space
    memmove(to->D + to->N, D+fromIdx, (N-fromIdx) * sizeof(T));
    to->N += N-fromIdx;
    N = fromIdx;
  }

  int crack(T const &v, int &i, int &L, int &R, bool sort_piece, CMP &cmp, Random &rng){
    assert(!next());            // it doesn't make sense crack a chained bucket!
    i = get_piece_by_value(v,L,R,cmp,rng);    // find the piece [L,R) containing v
    assert(L>=0 && L<=R && R<=N);        // range check
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
    // for (int pos = L; pos < R; pos++) if (!cmp(D[pos], v)) return pos;
    // return R;
    return std::lower_bound(D+L, D+R, v, cmp) - D;    // find the element v using binary search
  }

  bool erase(T const &v, CMP &cmp, Random &rng){
    int i,L,R,at=crack(v,i,L,R,false,cmp,rng);
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
    D[at] = D[--N];    // the deleted element has been shuffled out from the bucket
    I--;        // adjust the pending index

//      assert(check(D[0],false,D[0],false,cmp));
    return true;
  }
};


template <typename T, typename CMP  = std::less<T>, int MAX_BSIZE = 4096>
class Comb {
  class RangeError {};    // an exception class

  typedef Bucket<T, CMP> bucket_type;    // A COMB bucket.
  typedef pair<bucket_type*, bucket_type*> bucket_chain;       // Pointer to first and last bucket in the chain.
  // typedef btree::btree_map<T, bucket_chain> root_map;
  typedef pair<T, bucket_chain> root_entry;
  typedef std::map<T, bucket_chain> root_map;
  typedef typename root_map::const_iterator root_iterator;

  CMP cmp;     // The comparator functor.
  Random rng;  // The random number generator.
  root_map R;  // The root array, using Google's BTree.

  void set_root(T key, bucket_chain value) {
    assert(!R.count(key));
    R[key] = value;
  }

  // add a Bucket (bidx) to the root chain 'ridx'
  void _add(bucket_chain &c, bucket_type *b) {
    if (!c.first){                   // the root chain is empty.
      c.first = c.second = b;        // add the new bucket directly.
    } else {
      assert(!c.second->next());
      // assert(!B[Pe[ridx]].free());
      if (c.second->slack()) {
        b->moveToFromIdx(c.second, b->size() - std::min(b->size(), c.second->slack()));
      }
      if (b->size()) {
        c.second->set_next(b);
        c.second = b;
      } else {
        delete b;
      }
    }
    c.second->set_next(NULL);
  }

  T get_random_pivot(bucket_chain c) {  // pick the pivot near the median
    T rtmp[11]; int ntmp = 0, ni = 0;
    for (bucket_type *b = c.first; b; b = b->next(), ni++){
      b->clear_indexes();
      if (ntmp < 11){        // pick a random object from a stream!
        rtmp[ntmp++] = b->randomValue(rng);
      } else if (rng.nextInt(ni) < 11) {
        rtmp[rng.nextInt(11)] = b->randomValue(rng);
      }
    }
    while (ntmp < 11) {
      rtmp[ntmp++] = c.first->randomValue(rng);
    }
    std::nth_element(rtmp, rtmp+5, rtmp+11, cmp);
    return rtmp[5];      // the chosen pivot is here
  }

  // fusion
  pair<root_entry, root_entry> split_chain(bucket_chain c1) {
    const T &p = get_random_pivot(c1);        // split based on random pivot

    bucket_type *b = c1.first;
    bucket_chain c2 = c1 = make_pair(nullptr, nullptr);
    bucket_type *Lb = NULL, *Rb = NULL;
    int hi[MAX_BSIZE], lo[MAX_BSIZE];
    int nhi = 0, nlo = 0;

    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        Lb->fusion(Rb, hi, lo, nhi, nlo);
        if (nhi == 0){ _add(c1, Lb); Lb = NULL; }
        if (nlo == 0){ _add(c2, Rb); Rb = NULL; }
      } else if (!Lb){
        if (!b) break;
        Lb = b;
        b = b->next();
      } else if (nhi == 0) {
        assert(Lb);
        Lb->mark_hi(p, cmp, hi, nhi);
        if (nhi == 0){ _add(c1, Lb); Lb = NULL; }
      } else if (!Rb){
        if (!b) break;
        Rb = b;
        b = b->next();
      } else if (nlo == 0) {
        assert(Rb);
        Rb->mark_lo(p, cmp, lo, nlo);
        if (nlo == 0){ _add(c2, Rb); Rb = NULL; }
      } else {
        assert(0);
      }
    }

    if (Rb) { assert(!Lb); Lb = Rb; }
    if (Lb) {
      if (Lb->size()) {
        int i = Lb->partition(p, cmp);
        if (i == 0){
          _add(c2, Lb);
        } else if (i == Lb->size()) {
          _add(c1, Lb);
        } else {
          Rb = new bucket_type(Lb->capacity());
          Lb->moveToFromIdx(Rb, i);
          _add(c1, Lb);
          _add(c2, Rb);
        }
      } else {
        delete Lb;
      }
    }

    assert(!c1.second->next());
    assert(!c2.second->next());
//    assert(check());

    assert(valid_chain(c1));
    assert(valid_chain(c2));
    root_entry left_entry = make_pair(c1.first->data(0), c1);
    root_entry right_entry = make_pair(p, c2);
    return make_pair(left_entry, right_entry);
  }

  bool valid_chain(bucket_chain c) {
    return c.first && c.second && ((c.first->next() && c.first != c.second) || (!c.first->next() && c.first == c.second));
  }

  root_iterator break_chain(root_iterator it, T const &v) {
    // assert(check());
    while (it->second.first->next()) {
      // fprintf(stderr, "break_chain %lu\n", R.size());

      bool is_begin = it == R.begin();
      T left_key = it->first;
      bucket_chain left_chain = it->second;
      R.erase(it);

      pair<root_entry, root_entry> cs = split_chain(left_chain);  // randomly split the chain into two

      root_entry right_entry = cs.second;

      if (cmp(cs.first.first, left_key)) {
        // fprintf(stderr, "%d %d\n", cs.first.first, left_key);
        assert(is_begin); // Can only happen for the leftmost bucket.
        left_key = cs.first.first;
      }

      assert(cmp(left_key, right_entry.first));

      set_root(left_key, cs.first.second);
      it = find_root(left_key);
      assert(it->first == left_key);

      set_root(right_entry.first, right_entry.second);

      if (!(cmp(v, right_entry.first))) {
        // fprintf(stderr, "READJUST\n");
        // it = find_root(v);
        it++; // readjust root idx
        assert(it->first == right_entry.first);
        assert(it->second.first == right_entry.second.first);
        assert(it->second.second == right_entry.second.second);
      }
    }
    // assert(check());
    return it;
  }

  int size(bucket_chain c) {
    bucket_type *b = c.first;
    int ret = 0;
    while (b) {
      ret += b->size();
      b = b->next();
    }
    return ret;
  }

  int slack(bucket_chain c) {
    bucket_type *b = c.first;
    int ret = 0;
    while (b) {
      ret += b->slack();
      b = b->next();
    }
    return ret;
  }

public:
  class iterator {
   public:
    typedef Comb<T, CMP, MAX_BSIZE> crack_type;
    crack_type *crack;
    root_iterator root_iter;
    bucket_type *bucket;
    int idx;

    // STL specific typedefs.
    typedef iterator self_type;
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef int difference_type;
    typedef std::forward_iterator_tag iterator_category;


    iterator(crack_type *ct, root_iterator rit, bucket_type *b, int i): crack(ct), root_iter(rit), bucket(b), idx(i) {}
    // self_type operator++() { self_type i = *this; ptr_++; return i; }
    // self_type operator++(int junk) { ptr_++; return *this; }
    // const reference operator*() { return *ptr_; }
    // const pointer operator->() { return ptr_; }
    // bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
    // bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }

    bool operator==(iterator &that) {
      if (has_next() != that.has_next()) return false;
      return root_iter == that.root_iter && bucket == that.bucket && idx == that.idx;
    }

    bool operator!=(iterator &that) { return !(*this == that); }

    bool erase() { return false; }

    // bool debug() {
    //   fprintf(stderr,"iter: v=%d, ridx=%d/%d, bidx=%d/%d, idx=%d/%d\n",
    //     crack->get_bucket_element(bidx,idx),
    //     ridx,crack->root_size(),
    //     bidx, crack->num_of_buckets(),
    //     idx, crack->bucket_size(bidx));
    //   return false;
    // }

    bool has_next() {
      if (root_iter == crack->root_end()) return false;
      if (idx < bucket->size()) return true;  // The current bucket still have unread element.
      if (bucket->next()) {                   // Seek the next bucket in the bucket chain.
        bucket = bucket->next();
        assert(!bucket->empty());
        idx = 0;
        return true;
      }
      root_iter++;                                   // Seek the next root
      if (root_iter != crack->root_end()) {
        bucket = root_iter->second.first;
        idx = 0;
        return true;
      }
      return false; // no next element
    }

    bool next(T &res) {
      if (!has_next()) return false;
      res = bucket->data(idx++);
      return true;
    }

    bool prev(T &res) {
      // TODO fix this!
      if (idx == 0) return false;
      res = crack->data(--idx);
      return true;
    }

    // T* next() {
    //   if (!has_next()) return NULL;
    //   return crack->get_bucket_element_ptr(bidx,idx++);
    // }
  };

  Comb() {
    bucket_type *leaf = new bucket_type(MAX_BSIZE);
    set_root(0, make_pair(leaf, leaf));
  }
  
  void load(T const *arr, int n) {
    assert(!R.empty());
    root_iterator it = R.begin();
    bucket_chain c = it->second;
    assert(c.first == c.second);
    assert(c.second->empty());
    R.erase(it);

    int i = 0, CAP = c.second->capacity();
    while (i + CAP <= n) {
      i += c.second->bulk_insert(arr + i, CAP);
      bucket_type *nb = new bucket_type(CAP);
      c.second->set_next(nb);
      c.second = nb;
    }
    while (i < n) {
      c.second->insert(arr[i++]);
    }
    set_root(c.first->data(0), c);
  }

  root_iterator find_root(T const &v) const {
    root_iterator it = R.lower_bound(v);
    if (it == R.end()) {
      if (it != R.begin()) it--;
    } else if (v < it->first) {
      if (it != R.begin()) it--;
    }
    return it;
  }

  root_iterator root_end() {
    return R.end();
  }

  bucket_chain insert(T const &v, bucket_chain p) {
    // fprintf(stderr, "iins %d, size/cap = %d / %d\n", v, p.second->size(), p.second->capacity());
    if (p.second->is_full()) {
      // fprintf(stderr, "isfull %p %p, %p %p\n", p.first, p.second, p.first->next(), p.second->next());
      // assert(check());
      int new_cap = min(MAX_BSIZE, p.second->capacity() * (p.first == p.second ? 1 : 2));
      bucket_type *nb = new bucket_type(new_cap);
      p.second->set_next(nb);
      p.second = nb;
      // fprintf(stderr, "isfull2 %p %p, %p %p\n", p.first, p.second, p.first->next(), p.second->next());
    }
    // assert(check());
    // fprintf(stderr, "inserted to bucket = %p\n", p.second);
    p.second->insert(v);
    // assert(check());
    return p;
  }

  void insert(T const &v) {
    // assert(check());
    root_iterator it = find_root(v);
    // fprintf(stderr, "insert %d, root = %d\n", v, it->first);

    // if (v == 1846371397) {
    //   fprintf(stderr, "has next = %p, FL = %p %p\n", it->second.first->next(), it->second.first, it->second.second);
    //   it->second.first->debug("a", 0,0);
    //   it->second.second->debug("b", 0,0);
    // }

    bucket_chain c1 = it->second;
    bucket_chain c2 = insert(v, c1);
    if (c1 != c2) {
      // fprintf(stderr, "DELELELE, %p %p, c2 next = %p\n", c2.first, c2.second, c2.first->next());
      root_entry t = *it;
      R.erase(it);
      // fprintf(stderr, "INS KEY = %d, %p %p, c2 next = %p\n", t.first, c2.first, c2.second, c2.first->next());
      R[t.first] = c2;
    }
    // assert(check());
    // fprintf(stderr, "insert2 %d, root = %d\n", v, it->first);
  }

  int size() {
    int ret = 0;
    for (auto it : R) ret += size(it.second);
    return ret;
  }

  int slack() {
    int ret = 0;
    for (auto it : R) ret += slack(it.second);
    return ret;
  }

  bool erase(T const &v) {
    // 
    // if (v == 1489490129) { fprintf(stderr, "erase %d\n", v); assert(check()); }
    root_iterator it = find_root(v);
    const bucket_chain &c = it->second;
    if (c.first->empty() && c.first == c.second) return false;

    it = break_chain(it, v);
    bucket_type *b = it->second.first;
    bool ret = b->erase(v, cmp, rng);

    // TODO: merge with the left bucket if too small.

    if (b->empty()) {
      R.erase(it);
      delete b;
    } else {
      split_bucket(it, b);
    }

    // assert(check());
    // fprintf(stderr, "erase2 %d\n", v);

    return ret;
  }

  bool split_bucket(root_iterator it, bucket_type *b) {
    // To avoid having too many cracker indexes in a bucket.
    if (b->n_cracks() > 40 && b->capacity() == MAX_BSIZE) {
      // assert(check());
      fprintf(stderr, ".%lu", R.size());
      vector<pair<T, bucket_type*>> nb = b->split(cmp);
      // assert(check());

      if (it->first >= nb[0].first) {
        // Must be the leftmost root entry.
        // fprintf(stderr, "REPLACE LEFTMOST ERASE %d\n", it->first);
        R.erase(it);
        set_root(b->data(0), make_pair(b, b));
        it = find_root(b->data(0));
        // assert(it->first == b->data(0));
        // fprintf(stderr, "REPLACE LEFTMOST ADD %d\n", it->first);
      }
      // assert(check());

      for (auto it : nb) {
        // fprintf(stderr, "SET ROOT %d\n", it.first);
        // it.second->debug("he", 0, 0);
        bucket_chain c(it.second, it.second);
        while (c.second->next()) c.second = c.second->next();
        set_root(it.first, c);
        // assert(check());
      }
      // assert(check());
      return true;
    }
    return false;
  }

  /* TODO: lazy lower_bound */
  iterator lower_bound(T const &v, bool sort_piece = true) {
    static int nth = 0; nth++;
    // if (nth % 300 == 0) { fprintf(stderr, "lower_bound = %d\n", v); assert(check()); }
    root_iterator it = break_chain(find_root(v), v);
    // assert(check());
    bucket_type *b = it->second.first;
    int i, LL, RR, idx = b->crack(v, i, LL, RR, sort_piece, cmp, rng);
    // assert(check());

    if (0 && split_bucket(it, b)) {
      it++;
      if (v >= it->first) {
        b = it->second.first;
      } else {
        it--;
      }
      idx = b->crack(v, i, LL, RR, sort_piece, cmp, rng);
    }

    // assert(check());

    if (idx == b->size()) {
      // assert(check());
      it++;
      if (it != R.end()) {
        it = break_chain(it, v);
        b = it->second.first;
        assert(b->size() > 0);
        idx = b->crack(v, i, LL, RR, sort_piece, cmp, rng);
        // fprintf(stderr, "idx = %d, v = %d\n", idx, v);
        assert(idx == 0);
        // B[Pb[ridx]].nth(0,cmp,rng);          // just to crack it
      }
    }
    // assert(check());
    // if (ith > 1e9) assert(check());
    // fprintf(stderr, "lower_bound2 = %d\n", v);
    return iterator(this, it, b, idx);
  }

  // erase from [v1, v2)
  void erase(T const &v1, T const &v2){
    if (cmp(v2, v1)) throw RangeError();
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
    assert(!(it1.bucket->next()));
    assert(!it2.bucket->next());
    if (it1.root_iter == it2.root_iter){
      assert(it1.bucket == it2.bucket);
      return it2.idx - it1.idx;
    }
    int ret = it1.bucket->size() - it1.idx;
    for (it1.root_iter++; it1.root_iter != it2.root_iter; it1.root_iter++) {
      ret += size(it1.root_iter->second);
    }
    return ret + it2.idx;
  }

  void erase(iterator it1, iterator it2){ throw RangeError(); }

  bool check() {
    int i = 0;
    for (auto it = R.begin(); it != R.end(); i++) {
      bucket_type *b = it->second.first;
      if (!b->next() && it->second.first != it->second.second) {
        fprintf(stderr, "Head mismatch %p %p, %p\n", it->second.first, it->second.second, it->second.first->next());
        return false;
      }
      T lower = it->first;
      it++;
      T upper = (T){0};
      if (it != R.end()) upper = it->first;

      for (; b; b = b->next()) {
        if (!b->check(lower, i > 0, upper, it != R.end(), cmp)) {
//          for (int j=0; j<root_size(); j++)
//            fprintf(stderr,"root[%d] = %d\n",j,R[j]);
          fprintf(stderr,"Fail ridx = %d / %lu\n", i, R.size());
          return false;
        }
      }
    }
    return true;
  }

  int exists(T const &v, bool print=false) {
    int cnt = 0, i = 0;
    for (auto it = R.begin(); it != R.end(); i++, it++) {
      bucket_type *b = it->second.first;
      for (; b; b = b->next()) {
        int at = b->index_of(v, cmp);
        if (at != -1){
          if (print) fprintf(stderr,"v=%d, Exists at ridx = %d/%lu, key = %d, at=%d/%d, next=%p\n",
            v, i,R.size(), it->first, at,b->size(),b->next());
          b->debug("ignore",0,0);
          cnt++;
        }
      }
    }
    return cnt;
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
