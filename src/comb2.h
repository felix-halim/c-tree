#ifndef _COMB2_H_
#define _COMB2_H_

#include <stdio.h>
#include <string.h>

#include <functional>
#include <algorithm>
#include "random.h"

namespace comb {

using namespace std;

template<typename T,
  int MAX_SIZE    = 2048,           // Bucket capacity.
  int CRACK_AT    = 128,            // Creates a cracker index if no index exists in the span of CRACK_AT length.
  int DECRACK_AT  = 32>             // Destroys cracker index that spans less than DECRACK_AT length.

class Bucket {
  // Maximum number of cracker indices in a bucket.
  const static int MAX_INDEX = (MAX_SIZE + DECRACK_AT - 1) / DECRACK_AT - 1;

  // static_assert(MAX_SIZE < 65536, "MAX_SIZE must be less than 65536");
  static_assert(MAX_INDEX <= 64, "MAX_INDEX be at most 64");

  int size_;              // Number of data elements in D.
  int last_indexed_pos;   // Last indexed position in D.
  int next_bucket_id;     // Pointer to the next bucket ID.
  unsigned char nC;       // Number of cracker indices.
  unsigned long long S1;  // Sorted bits of each cracker index.
  int C[MAX_INDEX];       // Cracker index positions.
  T V[MAX_INDEX];         // Cached cracker values.
  T D[MAX_SIZE];          // Data elements.

 public:

  int size() const { return size_; }
  int free() const { return MAX_SIZE - size(); }
  int next() const { return next_bucket_id; }
  T* data() { return D; }
  T get(int i) const { assert(i >= 0 && i < size()); return D[i]; }
  void next(int id) { next_bucket_id = id; }
  void insert(T const &v) { assert(size() < MAX_SIZE); D[size_++] = v; }
  void init() { size_ = 0; next_bucket_id = -1; clear_indexes(); }
  void clear_indexes() { S1 = /* S2 = S3 = S4 = */ nC = last_indexed_pos = 0; }

  void piece_set_sorted(int i, bool sorted);
  bool piece_is_sorted(int i) const;
  void piece_set_unsorted_onwards(int i);

  void add_cracker_index(int at, int M);
  void remove_cracker_index(int at);

  // partition this bucket based on value v, destroying all cracker indices
  template<typename CMP>
  int partition(T const &v, CMP &cmp) {
    clear_indexes();
    assert(size() > 0 && size() <= MAX_SIZE);
    return partition(D, D + size(), v, cmp) - D;
  }

  template<typename CMP>
  T* partition(T *F, T *L, T const &v, CMP &cmp) {
    while (true){
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

  // move this bucket data in range [fromIdx, end) and append
  // it to the specified Bucket "to", destroying all cracker indices
  void moveToFromIdx(Bucket &to, int fromIdx){
    clear_indexes(); to.clear_indexes();    // destroy both buckets' cracker indices
    assert(size() > fromIdx);            // make sure there is something to move
    assert(to.size() + size() - fromIdx <= MAX_SIZE);    // make sure the receiver has enough space
    memmove(to.D + to.size(), D + fromIdx, (size() - fromIdx) * sizeof(T));
    to.size_ += size() - fromIdx;
    size_ = fromIdx;
  }

  template<typename CMP>
  int get_piece_by_index(int idx, int &L, int &R, CMP &cmp, Random &rng);

  template<typename CMP>
  T nth(int idx, CMP &cmp, Random &rng) {
    assert(next() == -1);      // it doesn't make sense crack a chained bucket!
    int L, R, i = get_piece_by_index(idx, L, R, cmp, rng);
    assert(L >= 0 && L <= R && R <= size());
    if (!piece_is_sorted(i)){    // sort the piece if it isn't
      std::sort(D+L,D+R,cmp);
      piece_set_sorted(i,true);
    }
    return D[idx];          // this is the correct nth element
  }

  T randomValue(Random &rng) const { return D[rng.nextInt(size())]; }

  template<typename CMP>
  void flush_pending_inserts(CMP &cmp);

  // Reorganizes the elements and returns a piece [L,R) that contains v.
  template<typename CMP>
  int get_piece_by_value(T const &v, int &L, int &R, CMP &cmp, Random &rng);

  template<typename CMP>
  int crack(T const &v, int &i, int &L, int &R, bool sort_piece, CMP &cmp, Random &rng);

  template<typename CMP>
  bool erase(T const &v, CMP &cmp, Random &rng);

  bool debug(const char *msg, int i, int j) const;

  // call this function to check the consistency of this Bucket structure
  template<typename CMP>
  bool check(T lo, bool useLo, T hi, bool useHi, CMP &cmp) const;

};


template<typename T, typename CMP>
static bool eq(T const &a, T const &b, CMP const &cmp) { return !cmp(a,b) && !cmp(b,a); }


template <typename T,
  typename CMP    = std::less<T>,
  bool USE_POS    = false,
  int MAX_SIZE    = 2048,
  int CRACK_AT    = 128,
  int DECRACK_AT  = 32>

class Comb {
  static_assert(DECRACK_AT * 2 <= CRACK_AT, "insufficient gap for crack decrack");
  typedef Bucket<T, MAX_SIZE, CRACK_AT, DECRACK_AT> Bucket_t;

  Random rng;   // The random number generator.
  CMP cmp;      // The comparator functor.
  T *R;         // The root data.
  Bucket_t *B;   // The buckets.
  int *Pb;      // Root pointers to first and last bucket.
  int *Pe;
  // int *bit;     // The cumulative number of elements in the root chains.
  int *F;       // Free bucket indices.
                // F[0] = num of free buckets
                // F[1] = nR = number of roots
                // F[2] = nB = number of buckets

  // Get a "free" bucket and initialize.
  int new_bidx();

  // Set free the 'bidx' Bucket.
  void free_bidx(int bidx);

  // Adds a Bucket (bidx) to the root chain 'ridx'.
  void _add(int ridx, int bidx);

  int size(int ridx);

  T get_random_pivot(int ridx);
  void split_chain(int ridx);

  int break_chain(int &ridx, T const &v);

  // find the correct root index (binary search)
  int find_ridx(T const &v) const;

  void destroy_root(int ridx);


 public:

  class iter {
   public:
    typedef Comb<T,CMP,USE_POS,MAX_SIZE,CRACK_AT,DECRACK_AT> crack_type;
    crack_type *crack;
    int ridx, bidx, idx;

    // STL specific typedefs.
    // typedef iter self_type;
    // typedef T value_type;
    // typedef T& reference;
    // typedef T* pointer;
    // typedef int difference_type;
    // typedef std::forward_iter_tag iter_category;


    iter(crack_type *ct, int r, int b, int i): crack(ct), ridx(r), bidx(b), idx(i) {}
    // self_type operator++() { self_type i = *this; ptr_++; return i; }
    // self_type operator++(int junk) { ptr_++; return *this; }
    // const reference operator*() { return *ptr_; }
    // const pointer operator->() { return ptr_; }
    // bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
    // bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }

    bool operator==(iter &that) {
      if (has_next() != that.has_next()) return false;
      return ridx==that.ridx && bidx==that.bidx && idx==that.idx;
    }

    bool operator!=(iter &that){ return !(*this == that); }
    bool erase(){ return false; }

    bool debug(){
      fprintf(stderr,"iter: v=%d, ridx=%d/%d, bidx=%d/%d, idx=%d/%d\n",
        crack->get_bucket_element(bidx,idx),
        ridx,crack->root_size(),
        bidx, crack->num_of_buckets(),
        idx, crack->bucket_size(bidx));
      return false;
    }

    bool has_next(){
      if (ridx >= crack->root_size()) return false;
      if (idx < crack->bucket_size(bidx)) return true;  // the current bucket still have unread element
      if (crack->bucket_next(bidx) != -1){        // seek the next bucket in the bucket chain
        bidx = crack->bucket_next(bidx);
        assert(crack->bucket_size(bidx) > 0);
        idx = 0;
        return true;
      }
      ridx++;                        // seek the next root
      if (ridx < crack->root_size()){
        bidx = crack->get_root_bidx(ridx);
        idx = 0;
        return true;
      }
      return false; // no next element
    }

    bool next(T &res){
      if (!has_next()) return false;
      res = crack->get_bucket_element(bidx,idx++);
      return true;
    }

    bool prev(T &res){
      // TODO fix this!
      if (idx == 0) return false;
      res = crack->get_bucket_element(bidx,--idx);
      return true;
    }

    T* next(){
      if (!has_next()) return NULL;
      return crack->get_bucket_element_ptr(bidx,idx++);
    }
  };


  Comb() {
    int nB = 1;
    allocate(nB);
    set_num_of_buckets(nB);
    clear();
  }
  
  ~Comb(){
    delete[] R;
    // if (USE_POS) delete[] bit;
    delete[] Pb;
    delete[] Pe;
    delete[] F;
    delete[] B;
  }

  int root_size() const { return F[1]; }
  int num_of_buckets() const { return F[2]; }
  void set_root_size(int nR){ F[1] = nR; }
  void set_num_of_buckets(int nB){ F[2] = nB; }
  int bucket_size(int bidx) const { return B[bidx].size(); }
  int bucket_next(int bidx) const { return B[bidx].next(); }
  int get_root_bidx(int ridx){
    assert(ridx>=0 && ridx<root_size());
    return Pb[ridx];
  }

  void allocate(int nB) {
    R = new T[nB];    // the root elements
    // bit = USE_POS? new int[nB] : NULL;  // the bit for position
    Pb = new int[nB]; // pointer to the first bucket
    Pe = new int[nB]; // pointer to the last bucket (for insert/append)
    F = new int[nB+3];  // F[0] is the starting index  of the free buckets
    B = new Bucket_t[nB]; // the bucket itself containing the data
  }

  void clear(){
    int nB = num_of_buckets();
    for (int i=3; i<nB+3; i++) F[i] = i-3;
    set_root_size(1);          // root has 1 element now
    F[0] = 2;
    Pb[0] = Pe[0] = new_bidx();      // allocate the first bucket to the root
  }

  T get_bucket_element(int bidx, int idx){
    assert(bidx >= 0 && bidx < num_of_buckets());
    assert(idx < MAX_SIZE);
    assert(idx >= 0 && idx < B[bidx].size());
    return B[bidx].get(idx);
  }

  template<typename X>
  void doubles(X *&oldX, int oldN, int newN) {
    X *newX = new X[newN];
    for (int i = 0; i < oldN; i++)
      newX[i] = oldX[i];
    delete[] oldX;
    oldX = newX;
  }


  // Doubles the number of buckets if needed.
  void doubles_num_of_buckets() {
    int nB = num_of_buckets();
    int nB2 = nB * 2;

    doubles(R, nB, nB2);
    // if (USE_POS) doubles(bit, nB, nB2);
    doubles(Pb, nB, nB2);
    doubles(Pe, nB, nB2);
    doubles(F, nB+3, nB2+3);
    doubles(B, nB, nB2);

    // Set free indices.
    for (int i=nB+3; i<nB2+3; i++) F[i] = i-3;

    set_num_of_buckets(nB2);
  }

  void insert(T const &v, bool use_pos = true) {
    int ridx = find_ridx(v);
    assert(Pb[ridx]!=-1);
    if (ridx == 0 && (B[Pb[0]].size() == 0 || cmp(v, R[0]))) R[0] = v;
    assert(Pe[ridx]!=-1);
    if (!B[Pe[ridx]].free()){
      int bidx = new_bidx();
      B[Pe[ridx]].next(bidx);
      Pe[ridx] = bidx;
    }
    assert(ridx+1>=root_size() || cmp(v, R[ridx+1]));
    B[Pe[ridx]].insert(v);
    // if (USE_POS && use_pos){
    //   assert(ridx>=0 && ridx<root_size());
    //   bit_update(ridx,1);
    // }
  }

  /* TODO: lazy lower_bound */
  iter lower_bound(T const &v, bool sort_piece=true) {
    int ridx = find_ridx(v);
    if (ridx==0 && B[Pb[ridx]].size()==0) return iter(this, 1,0,0);
    assert(ridx+1 >= root_size() || cmp(v,R[ridx+1]));
    assert(ridx==0 || !cmp(v,R[ridx]));
    int cnt = break_chain(ridx, v);
    // if (cnt) fprintf(stderr, "b");
    assert(ridx+1 >= root_size() || cmp(v,R[ridx+1]));
    assert(ridx==0 || !cmp(v,R[ridx]));
    assert(B[Pb[ridx]].next() == -1);
    int i,L,R,idx = B[Pb[ridx]].crack(v,i,L,R,sort_piece,cmp,rng);
    if (idx == B[Pb[ridx]].size()){
      ridx++;
      if (ridx < root_size()){
        cnt = 0;
        while (B[Pb[ridx]].next() != -1){
          split_chain(ridx);
          cnt++;
        }
        // if (cnt && USE_POS) bit_init();
        assert(B[Pb[ridx]].size() > 0);
        B[Pb[ridx]].nth(0, cmp, rng);          // just to crack it
        // fprintf(stderr, "hoho\n");
        return iter(this, ridx, Pb[ridx], 0);
      }
      // fprintf(stderr, "here \n");
      return iter(this, ridx, 0, 0);
    }
    return iter(this, ridx, Pb[ridx], idx);
  }

  bool erase(T const &v) {
    int ridx = find_ridx(v);
    if (ridx == 0 && B[Pb[ridx]].size()==0) return false;
    //int cnt =
    break_chain(ridx,v);
    bool ret = B[Pb[ridx]].erase(v,cmp,rng);
    if (ridx > 0){  // move to left if too small
      int p = Pb[ridx-1], c = Pb[ridx];
      assert(B[c].next() == -1);
      if (B[p].next()==-1 && B[c].size() && B[c].size() + B[p].size() < MAX_SIZE - 100)
        B[c].moveToFromIdx(B[p], 0);
    }
    if (B[Pb[ridx]].size()==0 && ridx){
      assert(B[Pb[ridx]].next() == -1);
      destroy_root(ridx);
    } else {
      if (B[Pb[ridx]].size() > 0 && eq(v, R[ridx], cmp)){
        R[ridx] = B[Pb[ridx]].nth(0, cmp, rng);
        // fix_root_cache();
      }
      // if (USE_POS && ret){
      //   assert(ridx >=0 && ridx < root_size());
      //   bit_update(ridx,-1);
      // }
    }
    return ret;
  }

  bool check() {
    for (int i=0,cnt=0; i<root_size(); i++){
      T upper = (T){0};
      if (i+1<root_size()) upper = R[i+1];
      if (USE_POS) {
        cnt += size(i);
        // if (cnt != bit_sum(i))
        //   fprintf(stderr,"i = %d, cnt = %d, %d\n",i,cnt,bit_sum(i));
        // assert(cnt == bit_sum(i));
      }
      for (int idx = Pb[i]; idx != -1; idx = B[idx].next()){
        if (!B[idx].check(R[i], true, upper, i+1 < root_size(), cmp)) {
//          for (int j=0; j<root_size(); j++)
//            fprintf(stderr,"root[%d] = %d\n",j,R[j]);
          fprintf(stderr,"Fail ridx = %d/%d, %d/%d\n",
            i, root_size(), idx, num_of_buckets());
          return false;
        }
      }
    }
    return true;
  }
};

#define BUCKET_TPL(type) \
  template<typename T, int MAX_SIZE, int CRACK_AT, int DECRACK_AT> \
  type Bucket<T, MAX_SIZE, CRACK_AT, DECRACK_AT>

#define BUCKET_TPLC(type) \
  template<typename T, int MAX_SIZE, int CRACK_AT, int DECRACK_AT> \
  template<typename CMP> \
  type Bucket<T, MAX_SIZE, CRACK_AT, DECRACK_AT>


BUCKET_TPL(void)::piece_set_sorted(int i, bool sorted) {
  assert(i >= 0 && i <= nC && nC <= MAX_INDEX);
  if (sorted) {
    S1 |= 1ULL << i;
  } else {
    S1 &= ~(1ULL << i);
  }
}

BUCKET_TPL(bool)::piece_is_sorted(int i) const {
  assert(i >= 0 && i <= nC && nC <= MAX_INDEX);
  return S1 & (1ULL << i);
}

BUCKET_TPL(void)::piece_set_unsorted_onwards(int i) {
  assert(i >= 0 && i <= nC && nC <= MAX_INDEX);
  S1 &= (1ULL << i) - 1;          // destroy sorted bit std::vector from i onwards
}

// Insert bit 0 at position {at} in {S}
static void insert_bit_at(unsigned long long &S, int at) {
  S = ((S << 1) & ~((((1ULL << at) - 1) << 1) | 1)) | (S & ((1ULL << at) - 1));
}

// Remove bit at position {at} in {S}
static void remove_bit_at(unsigned long long &S, int at) {
  S = ((S & ~((((1ULL << at) - 1) << 1) | 1)) >> 1) | (S & ((1ULL << at) - 1));
}

BUCKET_TPL(void)::add_cracker_index(int at, int M) {
  assert(at >= 0 && at <= nC && nC < MAX_INDEX);
  for (int i = nC-1; i >= at; i--) {
    C[i + 1] = C[i];
    V[i + 1] = V[i];
  }
  C[at] = M;
  V[at] = D[M];
  nC++;
  assert(at == 0 || C[at - 1] < C[at]);
  assert(at + 1 == nC || C[at] < C[at + 1]);
  insert_bit_at(S1, at);
}

BUCKET_TPL(void)::remove_cracker_index(int at) {
  assert(at >= 0 && at < nC && nC > 0);
  for (int i = at + 1; i < nC; i++) {
    C[i - 1] = C[i];
    V[i - 1] = V[i];
  }
  nC--;
  remove_bit_at(S1, at);
}

BUCKET_TPLC(void)::flush_pending_inserts(CMP &cmp) {
  // Indexed index should be less than the number of elements.
  assert(last_indexed_pos <= size());

  // No index yet, all the elements are considered "inserted".
  if (!nC) { last_indexed_pos = size(); S1 = 0; return; } // S2 = S3 = S4 =

  // Flushing only makes sense when the bucket is not chained.
  assert(next_bucket_id == -1);

  // IMPROVE: bulk insert? (Currently using Merge Completely).
  int minC = nC;
  // Inserts all pending elements (from last_indexed_pos to size).
  for (int j = last_indexed_pos; last_indexed_pos < size(); j = ++last_indexed_pos) {
    int i = nC - 1;
    T tmp = D[j];            // Store the pending tuple.
    for (; i >= 0 && cmp(tmp, V[i]); i--) {  // insert by shuffling through the cracker indices C
      int &L = C[i];         // Left boundary of this cracker piece.
      D[j] = D[L + 1];       // Replace the pending with the next to cracker boundary.
      D[L + 1] = D[L];       // Shift the cracker boundary to the right.
      j = L++;               // reposition the cracker piece separator.
    }
    D[j] = tmp;              // The pending tuple is now merged in.
    minC = std::min(minC, i + 1);      // Keep track the lowest piece that is touched.
  }

  piece_set_unsorted_onwards(minC);
}


template<typename T, typename CMP>
static T* doUpperCrack(T *lo, T *hi, T *pivot, CMP &cmp) {
  while (lo <= hi) {
    while (lo <= hi && !cmp(*pivot, *lo)) lo++;
    while (lo <= hi &&  cmp(*pivot, *hi)) hi--;
    if (lo < hi) std::iter_swap(lo, hi);
  }
  return lo;
}

template<typename T, typename CMP>
static T* doLowerCrack(T *lo, T *hi, T *pivot, CMP &cmp) {
  while (lo <= hi) {
    while (lo <= hi &&  cmp(*lo, *pivot)) lo++;
    while (lo <= hi && !cmp(*hi, *pivot)) hi--;
    if (lo < hi) std::iter_swap(lo, hi);
  }
  return lo;
}

template<typename T, typename CMP>
static T* partition(T *lo, T *hi, T *pivot, bool upperCrack, CMP &cmp) {
  assert(lo <= pivot && pivot <= hi);
  if (upperCrack) {
    std::iter_swap(lo, pivot);
    // fprintf(stderr, "upperCrack %d, %d\n", *lo, hi - lo);
    T* ret = doUpperCrack(lo + 1, hi, lo, cmp);
    std::iter_swap(ret, lo);
    return ret;
  }
  std::iter_swap(hi, pivot);
  T* ret = doLowerCrack(lo, hi - 1, hi, cmp);
  std::iter_swap(ret, hi);
  return ret;
}

// Partitions roughly in the middle satisfying the DECRACK_AT.
template<typename T, typename CMP>
static T* rough_middle_partition2(T *L, T *R, int min_gap, CMP &cmp, Random &rng) {
  T *lo = L;
  T *hi = R - 1;
  T *H = L + (R - L) / 2;
  assert(R - L >= min_gap * 2);
  while (lo < hi) {
    T* pivot = lo + rng.nextInt(hi - lo + 1);
    // fprintf(stderr, "lo = %d(%d), hi = %d(%d), pivot = %d(%d)\n", lo - L, *lo, hi - L, *hi, pivot - L, *pivot);
    pivot = partition(lo, hi, pivot, false, cmp);

    if (pivot - L >= min_gap) return pivot;
    if (R - pivot >= min_gap) return pivot;

    // fprintf(stderr, "pivot = %d(%d)\n", pivot - L, *pivot);
    if (pivot == lo) {
      pivot = partition(lo, hi, pivot, true, cmp);
      if (pivot > H) break;
      if (pivot <= H) lo = pivot;
      else hi = pivot - 1;
      continue;
    } else {
      assert(pivot > lo);
    }
    if (pivot < H) lo = pivot + 1;
    else if (pivot > H) hi = pivot - 1;
    else {
      assert(false);
      assert(hi == pivot);
      std::iter_swap(H, pivot);
      break;
    }
  }
  // int ret = H - L;
  // while (L < H) assert(*L++ <= *H);
  // R--;
  // while (H < R) assert(*H <= *R--);
  return H;
}

// Partitions roughly in the middle satisfying the DECRACK_AT.
template<typename T, typename CMP>
static T* rough_middle_partition(T *L, T *R, int min_gap, CMP &cmp, Random &rng) {
  T *arr = L;
  int idx[10000];
  int N = R - L, H = N / 2;
  int lo = 0, hi = N, ret = -1;
  assert(N >= min_gap * 2);
  while (true) {
    int n = hi - lo;
    // fprintf(stderr, "%d %d, min_gap = %d, gap left = %d, %d\n", lo, hi, min_gap, lo, N - hi);
    assert(n > 0);
    T* p = arr + lo + rng.nextInt(n);
    T pval = *p;
    if (lo >= min_gap && N - hi >= min_gap) {
      // Normal Crack
      ret = partition(arr + lo, arr + hi, p, false, cmp) - arr;
      // fprintf(stderr, "normal crac = %d, pos = %d < %d < %d\n", *p, lo, ret, hi);

      // for (int i = lo; i < hi; i++) {
      //   fprintf(stderr, "%d = %d\n", i, arr[i]);
      // }
      break;
    }
    // Fusion!
    int ilo = 0, ihi = n;
    for (int i = lo, _n = lo + n; i < _n; i++) {
      idx[ilo] = idx[ihi] = i;
      ilo += !cmp(arr[i], *p);
      ihi +=  cmp(arr[i], *p);
    }
    int nswap = std::min(ilo, ihi - n);
    for (int i = 0; i < nswap; i++) {
      int j = ihi - i - 1;
      // fprintf(stderr, "swapping %d (%d) <>  %d (%d), pivot = %d\n", arr[idx[i]], idx[i], arr[idx[j]], idx[j], pval);
      if (idx[i] >= idx[j]) { nswap = i; break; }
      swap(arr[idx[i]], arr[idx[j]]);
    }
    int j = ihi - nswap;
    assert(j >= 0);
    int next = std::min(idx[nswap], idx[j]);
    // fprintf(stderr, "n = %d, nswap = %d, %d %d, ilo = %d, ihi = %d, next = %d, pivot = %d\n",
    //   n, nswap, idx[nswap], idx[n + nswap], ilo, ihi - n, next, p - arr);
    if (nswap == 0) {
      std::iter_swap(p, arr + next);
      ret = next;
      break;
    }
    // for (int i = 0; i < N; i++) fprintf(stderr, "[%d]%d ", i, arr[i]);
    // fprintf(stderr, "\n" );
    if (next <= H) {
      lo = next;
    } else {
      hi = next;
    }
  }
  assert(ret != -1);
  // fprintf(stderr, "ret = %d\n", ret);
  // for (int i = 0; i < N; i++) fprintf(stderr, "[%d]%d ", i, arr[i]);

  // for (int i = 0; i <= ret; i++) {
  //   fprintf(stderr, "\ni = %d, L = %d <= %d", i, arr[i], arr[ret]);
  //   assert(arr[i] <= arr[ret]);
  // }
  // for (int i = ret; i < N; i++) {
  //   fprintf(stderr, "\ni = %d, R = %d >= %d", i, arr[i], arr[ret]);
  //   assert(arr[i] >= arr[ret]);
  // }
  // fprintf(stderr, "\n\n");
  return L + ret;
}

BUCKET_TPLC(int)::get_piece_by_value(T const &v, int &L, int &R, CMP &cmp, Random &rng) {
  flush_pending_inserts(cmp);
  // assert(check(D[0],false,D[0],false, cmp));
  int i = 0;
  while (i < nC && !cmp(v, V[i])) i++;  // Find the cracker indices that covers v
  L = (i == 0) ? 0 : C[i - 1];          // Left crack boundary.
  R = (i == nC) ? size() : C[i];        // Right crack boundary.
  while (R - L > CRACK_AT) {            // Narrow down the piece using DDR.
    int M = rough_middle_partition(D + L + (i ? 1 : 0), D + R, DECRACK_AT, cmp, rng) - D;
    // assert(abs((L + R) / 2 - M) <= 1);
    add_cracker_index(i, M);
    // fprintf(stderr,"CRACKING %d %d, [%d %d]\n", M, D[M], L, R);
    if (cmp(v, D[M])) R = M; else L = M, i++;  // Adjust the cracker index i.
  }
  assert(i >= 0 && i <= nC);
  // assert(check(D[0],false,D[0],false, cmp));
  return i;
}

BUCKET_TPLC(int)::get_piece_by_index(int idx, int &L, int &R, CMP &cmp, Random &rng) {
  assert(idx >= 0 && idx < size());        // index must fall within the bucket's size
  flush_pending_inserts(cmp);
  int i = 0;
  while (i<nC && idx >= C[i]) i++;    // find the cracker indices that covers the idx
  L = i==0? 0 : C[i-1];          // the left crack
  R = i==nC? size() : C[i];          // the right crack
  while (R-L > CRACK_AT){          // narrow down the piece using DDR
    int M = rough_middle_partition(D + L + (i ? 1 : 0), D + R, DECRACK_AT, cmp, rng) - D;
    // assert(abs((L + R) / 2 - M) <= 1);
    add_cracker_index(i,M);
    if (idx < M) R=M; else L=M, i++;  // adjust the cracker index
  }
  assert(i>=0 && i<=nC);
  // assert(check(D[0], false, D[0], false, cmp));
  return i;
}

BUCKET_TPLC(int)::crack(T const &v, int &i, int &L, int &R, bool sort_piece, CMP &cmp, Random &rng) {
  assert(next() == -1);                     // It doesn't make sense crack a chained bucket!
  // assert(check(D[0], false, D[0], false, cmp));
  i = get_piece_by_value(v,L,R,cmp,rng);    // Find the piece [L,R) containing v.
  assert(L >= 0 && L <= R && R <= size());
  // assert(check(D[0], false, D[0], false, cmp));
  if (!piece_is_sorted(i)) {
    if (sort_piece) {                       // Sort the piece if requested.
      std::sort(D + L, D + R, cmp);
      piece_set_sorted(i, true);
    } else {
      for (int at = L; at < R; at++)
        if (eq(D[at], v, cmp)) return at;
      return R;
    }
  }
  // assert(check(D[0], false, D[0], false, cmp));
  return std::lower_bound(D + L, D + R, v, cmp) - D;
}


BUCKET_TPLC(bool)::erase(T const &v, CMP &cmp, Random &rng) {
  int i, L, R, at = crack(v, i, L, R, false, cmp, rng);
  if (at >= R || !eq(D[at], v, cmp)) return false;  // The element to be erased is not found!

  // Decrack this cracker piece (it becomes too small) or
  // if the deleted element index is a cracker index
  if (nC && R - L + 1 <= DECRACK_AT) {
    remove_cracker_index((i > 0) ? --i : i);
  } else if (i > 0 && at == L) {
    at = L + 1;
    for (int j = L + 2; j < R; j++)    // find a replacement element for the cracker index
      if (cmp(D[j], D[at])) at = j;  // that is the smallest element in the piece (L,R)
    std::swap(D[L], D[at]);
    piece_set_sorted(i, false);
    V[i - 1] = D[L];
  }

  assert(at < R && eq(D[at], v, cmp));    // the element v must be found!

  // IMPROVE: use pending delete? antimatter?
  piece_set_unsorted_onwards(i);      // unset the sorted bit i onwards
  assert(i == nC || cmp(D[at], V[i]));
  for (int j = i; j < nC; j++) {        // shuffle out the deleted element
    R = C[j]--;
    D[at] = D[R - 1];
    D[R - 1] = D[R];
    at = R;
  }
  D[at] = D[--size_];   // The deleted element has been shuffled out from the bucket.
  last_indexed_pos--;   // Adjust the pending index.

  // assert(check(D[0],false,D[0],false,cmp));
  return true;
}

BUCKET_TPL(bool)::debug(const char *msg, int i, int j) const {
  for (int k=0; k<nC; k++)
    fprintf(stderr,"C[%d/%d] = %d, %d (sorted = %d)\n",
      k,nC,C[k],(int)D[C[k]],piece_is_sorted(k));
  fprintf(stderr,"%s : i=%d/N=%d, j=%d/nC=%d, D[i,i+1] = %d, %d; I=%d, N=%d, next=%d\n",
    msg, i,size(), j,nC, (int)D[i],(int)D[i+1], last_indexed_pos,size(),next_bucket_id);
  return false;
}

BUCKET_TPLC(bool)::check(T lo, bool useLo, T hi, bool useHi, CMP &cmp) const {
  if (useLo) for (int i = 0; i < size(); i++) if (cmp(D[i], lo)) {
    fprintf(stderr,"D[%d] = %d, lo = %d\n", i, D[i], lo);
    return debug("useLo failed", i, 0);
  }
  if (useHi) for (int i = 0; i < size(); i++) if (!cmp(D[i], hi)) {
    fprintf(stderr,"D[%d] = %d, hi = %d\n", i, D[i], hi);
    return debug("useHi failed", i, 0);
  }
  for (int i=0,j=0; i<last_indexed_pos; i++){                    // check cracker indices
    if (j<nC && C[j]==i) assert(eq(V[j],D[i],cmp)), lo = D[i], j++;
    if (piece_is_sorted(j) && (j==nC? (i+1<last_indexed_pos) : (i<C[j])) && cmp(D[i+1], D[i]))
      return debug("sortedness violation", i,j);
    if (j>0 && cmp(D[i], D[C[j-1]])) return debug("lower bound fail", i,j);
    if (j<nC && cmp(D[C[j]], D[i])) return debug("upper bound fail", i,j);
    // if (j<nC && !(cmp(D[i], D[C[j]]))) return debug("upper bound fail", i,j);
  }
  return true;
}

#undef BUCKET_TPL
#undef BUCKET_TPLC



/********************************
 *           COMB               *
 ********************************/

#define COMB_TPL(type) \
  template<typename T, typename CMP, bool USE_POS, int MAX_SIZE, int CRACK_AT, int DECRACK_AT> \
  type Comb<T, CMP, USE_POS, MAX_SIZE, CRACK_AT, DECRACK_AT>


COMB_TPL(int)::new_bidx() {
  ++F[0];
  if (F[0] >= num_of_buckets()+3) {
    // fprintf(stderr, "F[0] = %d, num_of_buckets = %d\n", F[0], num_of_buckets()+3);
    doubles_num_of_buckets();
  }
  assert(F[0] < num_of_buckets()+3);
  int bidx = F[F[0]];
  B[bidx].init();
  return bidx;
}

COMB_TPL(void)::free_bidx(int bidx) {
  assert(F[0] > 0);
  F[F[0]--] = bidx;
  // TODO: shrink requires quite a number of data movement
  // if (F[0] < num_of_buckets()/2-2) resize_buckets(false);
}

COMB_TPL(void)::_add(int ridx, int bidx) {
  if (Pb[ridx] == -1){                        // the root chain is empty
    Pb[ridx] = Pe[ridx] = bidx;                    // add the new bucket directly
  } else {
    assert(B[Pe[ridx]].next() == -1);
    // assert(!B[Pe[ridx]].free());
    if (B[Pe[ridx]].free())
      B[bidx].moveToFromIdx(B[Pe[ridx]],
        B[bidx].size() - std::min(B[bidx].size(), B[Pe[ridx]].free()));
    if (B[bidx].size()){
      B[Pe[ridx]].next(bidx);
      Pe[ridx] = bidx;
    } else {
      free_bidx(bidx);
    }
  }
  B[Pe[ridx]].next(-1);
}

COMB_TPL(int)::size(int ridx) {
  assert(ridx >=0 && ridx < root_size());
  int ret = 0;
  for (int i=Pb[ridx]; i!=-1; i=B[i].next())
    ret += B[i].size();
  return ret;
}

COMB_TPL(T)::get_random_pivot(int ridx) {  // pick the pivot near the median
  T rtmp[11]; int ntmp = 0;
  for (int i=Pb[ridx],ni=0; i!=-1; i=B[i].next(),ni++){
    B[i].clear_indexes();
    if (ntmp < 11){        // pick a random object from a stream!
      rtmp[ntmp++] = B[i].randomValue(rng);
    } else if (rng.nextFloat()*ni < 11){
      rtmp[rng.nextInt(11)] = B[i].randomValue(rng);
    }
  }
  while (ntmp < 11)
    rtmp[ntmp++] = B[Pb[ridx]].randomValue(rng);
  std::nth_element(rtmp, rtmp+5, rtmp+11, cmp);
  return rtmp[5];      // the chosen pivot is here
}


template<typename T, typename CMP>
void mark_hi(T* D, int N, T const &P, CMP &cmp, int *hi, int &nhi){
  for (int i = 0; i < N; i++) {
    hi[nhi] = i;
    nhi += !cmp(D[i], P);
  }
}

template<typename T, typename CMP>
void mark_lo(T* D, int N, T const &P, CMP &cmp, int *lo, int &nlo){
  for (int i = 0; i < N; i++) {
    lo[nlo] = i;
    nlo += cmp(D[i], P);
  }
}

template<typename T>
static void fusion(T *Lp, T *Rp, int *hi, int *lo, int &nhi, int &nlo) {
  int m = std::min(nhi, nlo); assert(m > 0);
  int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
  nhi -= m; nlo -= m;
  while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
}

COMB_TPL(void)::split_chain(int ridx) {
  assert(ridx >=0 && ridx < root_size());
  assert(root_size() < num_of_buckets());

  for (int i=root_size()-1; i>ridx; i--){              // make a room for the new root chain (ridx+1)
    R[i+1] = R[i];
    Pb[i+1] = Pb[i];
    Pe[i+1] = Pe[i];
  }
  set_root_size(root_size()+1);

  const T &p = R[ridx+1] = get_random_pivot(ridx);        // split based on random pivot
  assert(cmp(R[ridx], R[ridx+1]));                // root chain must have different value

  // fix_root_cache();

  int i = Pb[ridx];
  Pb[ridx] = Pb[ridx+1] = Pe[ridx] = Pe[ridx+1] = -1;
  int Lb=-1, Rb=-1, *hi=new int[MAX_SIZE], *lo=new int[MAX_SIZE], nhi=0, nlo=0;
  while (true){
    if (nhi && nlo){
      assert(Lb != -1 && Rb != -1);
      fusion(B[Lb].data(), B[Rb].data(), hi, lo, nhi, nlo);
      if (nhi==0){ _add(ridx, Lb); Lb = -1; }
      if (nlo==0){ _add(ridx+1, Rb); Rb = -1; }
    } else if (Lb == -1){
      if (i == -1) break;
      i = B[Lb = i].next();
    } else if (nhi == 0){
      assert(Lb!=-1);
      mark_hi(B[Lb].data(), B[Lb].size(), p, cmp, hi, nhi);
      if (nhi == 0){ _add(ridx, Lb); Lb = -1; }
    } else if (Rb == -1){
      if (i == -1) break;
      i = B[Rb = i].next();
    } else if (nlo == 0){
      assert(Rb!=-1);
      mark_lo(B[Rb].data(), B[Rb].size(), p, cmp, lo, nlo);
      if (nlo == 0){ _add(ridx+1, Rb); Rb = -1; }
    } else {
      assert(0);
    }
  }
  delete[] hi;
  delete[] lo;

  if (Rb!=-1){ assert(Lb==-1); Lb = Rb; }
  if (Lb!=-1){
    if (B[Lb].size()){
      int i = B[Lb].partition(p, cmp);
      if (i == 0){
        _add(ridx+1, Lb);
      } else if (i == B[Lb].size()) {
        _add(ridx, Lb);
      } else {
        Rb = new_bidx();
        B[Lb].moveToFromIdx(B[Rb], i);
        _add(ridx, Lb);
        _add(ridx+1, Rb);
      }
    } else {
      free_bidx(Lb);
    }
  }

  assert(Pe[ridx]>=0 && Pe[ridx] < num_of_buckets());
  assert(Pe[ridx+1]>=0 && Pe[ridx+1] < num_of_buckets());
  assert(B[Pe[ridx]].next() == -1);
  assert(B[Pe[ridx+1]].next() == -1);
  assert(B[Pb[ridx]].next() == -1 || B[Pb[ridx]].size() == MAX_SIZE);
  assert(B[Pb[ridx+1]].next() == -1 || B[Pb[ridx+1]].size() == MAX_SIZE);
//    assert(check());
}

COMB_TPL(int)::break_chain(int &ridx, T const &v) {
  int cnt = 0;
  if (!(ridx >=0 && ridx < root_size())) fprintf(stderr, "ridx = %d, %d\n", ridx, root_size());
  assert(ridx >=0 && ridx < root_size());
  assert(Pb[ridx]>=0 && Pb[ridx]<num_of_buckets());
  assert(ridx==0 || B[Pb[ridx]].size() > 0);
  while (B[Pb[ridx]].next() != -1){
//      fprintf(stderr,"break ridx=%d, %d\n",ridx,size(ridx));
    split_chain(ridx);  // randomly split the chain into two
//      fprintf(stderr,"borken\n");
    assert(ridx+1 < root_size());
    assert(Pb[ridx+1]>=0 && Pb[ridx+1]<num_of_buckets());
    assert(B[Pb[ridx+1]].size() > 0);
    if (!(cmp(v, R[ridx+1]))) ridx++;  // readjust root idx
    cnt++;
  }
  // if (USE_POS && cnt > 0) bit_init();
//    fprintf(stderr,"broken ridx=%d, %d\n",ridx,cnt);
  return cnt;
}

COMB_TPL(int)::find_ridx(T const &v) const {
  if (root_size() <= 1) return 0;
  int lo = 0, hi = root_size() - 1, res = 0;
  while (lo<=hi){
    int mid = lo + ((hi-lo)>>1);
    if (cmp(v, R[mid])) hi = mid-1;
    else res = mid, lo = mid+1;
  }
  assert(res==0 || !cmp(v,R[res]));
  return res;
}

COMB_TPL(void)::destroy_root(int ridx) {
//    fprintf(stderr,"DESTROY ROOT %d\n",ridx);
  assert(B[Pb[ridx]].next() == -1);
  assert(Pb[ridx]== Pe[ridx]);
  free_bidx(Pb[ridx]);
  for (int i=ridx+1; i<root_size(); i++){
    R[i-1] = R[i];
    Pb[i-1] = Pb[i];
    Pe[i-1] = Pe[i];
  }
  set_root_size(root_size()-1);
  // if (USE_POS) bit_init();
  // fix_root_cache();
}


#undef COMB_TPL

/*
  void piece_set_sorted(int i, bool sorted){
    assert(i>=0 && i<MAXC);
    if (sorted){
          if (i<64)   S1 |= 1ULL << i;
      else   if (i<128)  S2 |= 1ULL << (i-64);
      else   if (i<192)  S3 |= 1ULL << (i-128);
      else          S4 |= 1ULL << (i-192);
    } else {
          if (i<64)  S1 &= ~(1ULL << i);
      else  if (i<128)  S2 &= ~(1ULL << (i-64));
      else  if (i<192)  S3 &= ~(1ULL << (i-128));
      else        S4 &= ~(1ULL << (i-192));
    }
  }

  bool piece_is_sorted(int i) const {
    assert(i>=0 && i<MAXC);
    if (i<64)  return S1 & (1ULL << i);
    if (i<128)  return S2 & (1ULL << (i-64));
    if (i<192)  return S3 & (1ULL << (i-128));
          return S4 & (1ULL << (i-192));
  }

  void piece_set_unsorted_onwards(int i){
    if (i<64){
      S1 &= (1ULL<<i)-1;          // destroy sorted bit std::vector from i onwards
      S2 = S3 = S4 = 0;
    } else if (i<128){
      S2 &= (1ULL<<(i-64))-1;          // destroy sorted bit std::vector from i onwards
      S3 = S4 = 0;
    } else if (i<192){
      S3 &= (1ULL<<(i-128))-1;          // destroy sorted bit std::vector from i onwards
      S4 = 0;
    } else {
      S4 &= (1ULL<<(i-192))-1;          // destroy sorted bit std::vector from i onwards
    }
  }

  void add_cracker_index(int at, int M){
    assert(at>=0 && at<=nC && nC<MAXC-1);
    for (int i=nC-1; i>=at; i--){
      C[i+1] = C[i];
      V[i+1] = V[i];
    }
    C[at] = M;
    V[at] = D[M];
    nC++;
    // assert(nC < MAXC);
    assert(at==0 || C[at-1] < C[at]);
    assert(at+1==nC || C[at] < C[at+1]);

    if (at<64){
      S4 = (S4<<1) | (piece_is_sorted(191)?1:0);
      S3 = (S3<<1) | (piece_is_sorted(127)?1:0);
      S2 = (S2<<1) | (piece_is_sorted(63)?1:0);
      insert_bit_at(S1, at);
    } else if (at<128){
      S4 = (S4<<1) | (piece_is_sorted(191)?1:0);
      S3 = (S3<<1) | (piece_is_sorted(127)?1:0);
      insert_bit_at(S2, at-64);
    } else if (at<192){
      S4 = (S4<<1) | (piece_is_sorted(191)?1:0);
      insert_bit_at(S3, at-128);
    } else {
      insert_bit_at(S4, at-192);
    }
  }

  void remove_cracker_index(int at){
    assert(at>=0 && at<nC);
    for (int i=at+1; i<nC; i++){
      C[i-1] = C[i];
      V[i-1] = V[i];
    }
    nC--;
    // assert(nC>=0);
    if (at<64){
      remove_bit_at(S1,at);
      piece_set_sorted(63, piece_is_sorted(64)); S2 >>= 1;
      piece_set_sorted(127, piece_is_sorted(128)); S3 >>= 1;
      piece_set_sorted(191, piece_is_sorted(192)); S4 >>= 1;
    } else if (at<128){
      remove_bit_at(S2,at-64);
      piece_set_sorted(127, piece_is_sorted(128)); S3 >>= 1;
      piece_set_sorted(191, piece_is_sorted(192)); S4 >>= 1;
    } else if (at<192){
      remove_bit_at(S3,at-128);
      piece_set_sorted(191, piece_is_sorted(192)); S4 >>= 1;
    } else {
      remove_bit_at(S4,at-192);
    }
  }












public:
  T* getp(int i) { assert(i>=0 && i<N); return &D[i]; }
  void bulk_insert(T const *v, int sidx, int length){
    assert(N==0);
    for (int j = 0; j < length; j++) {
      D[N++] = v[sidx++];
    }
  }





  int index_of(T const &v, CMP &cmp) const {
    for (int i=0; i<N; i++) if (equal(D[i],v,cmp)) return i;
    return -1;
  }



};






  void bit_update(int idx, int val){
    assert(idx>=0 && idx < root_size());
    for (idx++; idx <= root_size(); ){
      bit[idx-1] += val;
      idx += (idx & -idx);
    }
  }

  int bit_sum(int idx){
    if (idx<0) return 0;
    assert(idx>=0 && idx<root_size());
    int sum = 0; idx++;
    while (idx > 0){
      sum += bit[idx-1];
      idx -= (idx & -idx);
    }
    return sum;
  }






  T* get_bucket_element_ptr(int bidx, int idx) {
    assert(bidx >= 0 && bidx < num_of_buckets());
    assert(idx < MAX_SIZE);
    assert(idx >= 0 && idx < B[bidx].size());
    return B[bidx].getp(idx);
  }

public:


  void load(T const *v, int n) {
    int i = 0;
    while (i + MAX_SIZE <= n) {
      assert(B[Pe[0]].free() == MAX_SIZE);
      B[Pe[0]].bulk_insert(v, i, MAX_SIZE);
      i += MAX_SIZE;
      int bidx = new_bidx();
      B[Pe[0]].next(bidx);
      Pe[0] = bidx;
    }
    while (i < n) {
      B[Pe[0]].insert(v[i++]);
    }
  }

  void bit_init(){
    assert(USE_POS);
    memset(bit,0,sizeof(int)*root_size());
    for (int i=0; i<root_size(); i++) bit_update(i,size(i));
  }

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
    assert(ridx+1>=root_size() || cmp(v, R[ridx+1]));
    B[Pe[ridx]].insert(v);
    if (USE_POS && use_pos){
      if (cnt) bit_init();
      else {
        assert(ridx>=0 && ridx<root_size());
        bit_update(ridx,1);
      }
    }
  }

  int size(){
    int ret = 0;
    for (int i=0; i<root_size(); i++)
      for (int idx = Pb[i]; idx!=-1; idx=B[idx].next())
        ret += B[idx].size();
    return ret;
  }

  int exists(T const &v, bool print=false){
    int cnt = 0;
    for (int i=0; i<root_size(); i++)
      for (int idx = Pb[i]; idx!=-1; idx=B[idx].next()){
        int at = B[idx].index_of(v,cmp);
        if (at != -1){
          if (print) fprintf(stderr,"v=%d, Exists at ridx = %d/%d, idx=%d/%d, at=%d/%d, next=%d\n",
            v, i,root_size(), idx,num_of_buckets(), at,B[idx].size(),B[idx].next());
          B[idx].debug("ignore",0,0);
          cnt++;
        }
      }
    return cnt;
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

  iterator begin(){
    return iterator(this, 1, 0, 0);
  }


  // erase from [v1, v2)
  void erase(T const &v1, T const &v2){
    if (cmp(v2, v1)) throw RangeError();
    int i1 = find_ridx(v1);
    int i2 = find_ridx(v2);
    break_chain(i2,v2);
    break_chain(i1,v1);

    // TODO: destroy B from i1 to i2-1
    if (USE_POS) bit_init();
  }

  int count(T const &v1, T const &v2){
    assert(!cmp(v2,v1));
    iterator it1 = lower_bound(v1);
    iterator it2 = lower_bound(v2);
    return count(it1,it2);
  }

  int rank(iterator it){
    if (it.ridx==0) return it.idx;
    return bit_sum(it.ridx - 1) + it.idx;
  }

  int count(iterator it1, iterator it2){
    if (USE_POS) return rank(it2) - rank(it1);
    assert(B[Pb[it1.ridx]].next() == -1);
    assert(B[Pb[it2.ridx]].next() == -1);
    if (it1.ridx == it2.ridx){
      assert(it1.bidx == it2.bidx);
      return it2.idx - it1.idx;
    }
    assert(it1.ridx < it2.ridx);
    int ret = B[it1.bidx].size() - it1.idx;
    for (int i=it1.ridx+1; i<it2.ridx; i++) ret += size(i);
    return ret + it2.idx;
  }

  void erase(iterator it1, iterator it2){ throw RangeError(); }

};
*/


}

#endif
