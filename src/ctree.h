#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <chrono>
#include <algorithm>
#include "random.h"

using namespace std;
using namespace chrono;

namespace ctree {

#define INTERNAL_BSIZE 64   // Must be power of two.
#define LEAF_BSIZE 2048     // Must be power of two.
#define MAX_INDEX 64
#define CRACK_AT 64
#define DECRACK_AT 32

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

int nLeaves, nInternals, nCap, nDes;
Random rng;

class Bucket {
 protected:
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N, cap;

 public:
  int get_cap() { return cap; }
  int size() { return N; }
  bool is_full() { return size() == cap; }
  bool is_leaf() { return pending_insert >= 0; }
  void optimize();
  int debug(int depth);
  bool check();
};

class LeafBucket : public Bucket {
 protected:
  int *D;
  int last_indexed_pos;   // Last indexed position in D.
  unsigned char nC;       // Number of cracker indices.
  unsigned long long S1;  // Sorted bits of each cracker index.
  int C[MAX_INDEX];       // Cracker index positions.
  int V[MAX_INDEX];       // Cached cracker values.
  LeafBucket *next, *tail;  // Store pending inserts in a linked list.

 public:
  ~LeafBucket();
  LeafBucket(int cap);
  void init(int cap);
  LeafBucket* next_bucket() { return next; }
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }

  void clear_indexes() { S1 = /* S2 = S3 = S4 = */ nC = last_indexed_pos = 0; }
  void piece_set_sorted(int i, bool sorted);
  bool piece_is_sorted(int i) const;
  void piece_set_unsorted_onwards(int i);
  void add_cracker_index(int at, int M);
  void remove_cracker_index(int at);

  int get_piece_by_index(int idx, int &L, int &R);

  // T nth(int idx) {
  //   assert(next() == -1);      // it doesn't make sense crack a chained bucket!
  //   int L, R, i = get_piece_by_index(idx, L, R);
  //   assert(L >= 0 && L <= R && R <= size());
  //   if (!piece_is_sorted(i)){    // sort the piece if it isn't
  //     std::sort(D+L,D+R);
  //     piece_set_sorted(i,true);
  //   }
  //   return D[idx];          // this is the correct nth element
  // }

  // T randomValue(Random &rng) const { return D[rng.nextInt(size())]; }

  void flush_pending_inserts();

  // Reorganizes the elements and returns a piece [L,R) that contains v.
  int get_piece_by_value(int &v, int &L, int &R);

  int crack(int &v, int &i, int &L, int &R, bool sort_piece);

  // template<typename T, typename CMP>
  // bool erase(int &v);

  // bool debug(const char *msg, int i, int j) const;

  // // call this function to check the consistency of this Bucket structure
  // template<typename T, typename CMP>
  // bool check(T lo, bool useLo, T hi, bool useHi) const;



  void leaf_insert(int v);
  void leaf_split(vector<pair<int, LeafBucket*>> &nbs);
  void leaf_optimize();
  int promote_last();
  int lower_pos(int value);
  pair<bool,int> leaf_lower_bound(int value);
  LeafBucket* detach_and_get_next();
  void add_chain(LeafBucket *b);
  void distribute_values(int pivot, LeafBucket* chain[2]);
  LeafBucket* transfer_to(LeafBucket *b, int pivot);
};

class InternalBucket : public LeafBucket {
 protected:
  Bucket *C[INTERNAL_BSIZE + 1];

 public:
  ~InternalBucket();
  InternalBucket();
  InternalBucket(int promotedValue, Bucket *left, Bucket *right);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  int child_pos(int value);
  Bucket*& child_bucket(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
};



vector<LeafBucket*> free_leaves[30];

LeafBucket* new_leaf(int cap) {
  for (int i = 2; ; i++) {
    if ((1 << i) == cap) {
      if (free_leaves[i].empty()) {
        free_leaves[i].push_back(new LeafBucket(cap));
      }
      LeafBucket* leaf = free_leaves[i].back();
      leaf->init(cap);
      free_leaves[i].pop_back();
      return leaf;
    }
  }
}

void delete_leaf(LeafBucket *b) {
  for (int i = 2; ; i++) {
    if ((1 << i) == b->get_cap()) {
      free_leaves[i].push_back(b);
      break;
    }
  }
}

void Bucket::optimize() {
  if (is_leaf()) {
    ((LeafBucket*) this)->leaf_optimize();
  } else {
    for (int i = 0; i <= N; i++) {
      ((InternalBucket*) this)->child(i)->optimize();
    }
  }
}

int Bucket::debug(int depth) {
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  LeafBucket *b = (LeafBucket *) this;
  int sz = 0;
  while (b) {
    sz += b->size();
    fprintf(stderr, "N = %d (%d%s), ", b->size(), b->pending_insert, b->is_leaf() ? ", leaf" : "");
    for (int i = 0; i < b->size(); i++) {
      fprintf(stderr, "%d ", b->data(i));
    }
    if ((b = b->next_bucket())) {
      fprintf(stderr, "------ ");
    }
  }
  fprintf(stderr, "\n");
  if (!is_leaf()) {
    for (int i = 0; i <= ((LeafBucket*) this)->size(); i++) {
      sz += ((InternalBucket*) this)->child(i)->debug(depth + 1);
    }
  }
  for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
  fprintf(stderr, "size = %d\n", sz);
  return sz;
}

bool Bucket::check() {
  assert(N >=0 && N <= cap);
  assert(is_leaf() || ((InternalBucket*) this)->child(N));
  for (int i = 0; i < N; i++) {
    if (is_leaf()) {
    } else if (i > 0) {
      assert(((InternalBucket*) this)->child(i));
      assert(((InternalBucket*) this)->data(i - 1) < ((InternalBucket*) this)->data(i));
    }
  }
  return true;
}


LeafBucket::LeafBucket(int cap) {
  this->cap = cap;
  D = new int[cap];
  init(cap);
}

void LeafBucket::init(int cap) {
  assert(this->cap == cap);
  N = 0;
  pending_insert = 0;
  next = tail = NULL;
  nCap += cap;
  nLeaves++;
}

LeafBucket::~LeafBucket() {
  nCap -= cap;
  nLeaves--;
  nDes++;
}

void LeafBucket::leaf_insert(int value) {
  assert(is_leaf());
  assert(N >= 0);
  pending_insert++;
  if (!is_full()) {
    D[N++] = value;
  } else {
    if (!tail) {
      assert(cap == LEAF_BSIZE);
      add_chain(new_leaf(LEAF_BSIZE));
    } else if (tail->is_full()) {
      add_chain(new_leaf(LEAF_BSIZE));
    }
    tail->D[tail->N++] = value;
    tail->pending_insert++;
  }
  // assert(check());
}




void LeafBucket::piece_set_sorted(int i, bool sorted) {
  assert(i >= 0 && i <= nC && nC <= MAX_INDEX);
  if (sorted) {
    S1 |= 1ULL << i;
  } else {
    S1 &= ~(1ULL << i);
  }
}

bool LeafBucket::piece_is_sorted(int i) const {
  assert(i >= 0 && i <= nC && nC <= MAX_INDEX);
  return S1 & (1ULL << i);
}

void LeafBucket::piece_set_unsorted_onwards(int i) {
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

void LeafBucket::add_cracker_index(int at, int M) {
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

void LeafBucket::remove_cracker_index(int at) {
  assert(at >= 0 && at < nC && nC > 0);
  for (int i = at + 1; i < nC; i++) {
    C[i - 1] = C[i];
    V[i - 1] = V[i];
  }
  nC--;
  remove_bit_at(S1, at);
}

void LeafBucket::flush_pending_inserts() {
  // Indexed index should be less than the number of elements.
  assert(last_indexed_pos <= size());

  // No index yet, all the elements are considered "inserted".
  if (!nC) { last_indexed_pos = size(); S1 = 0; return; } // S2 = S3 = S4 =

  // Flushing only makes sense when the bucket is not chained.
  assert(!next);

  // IMPROVE: bulk insert? (Currently using Merge Completely).
  int minC = nC;
  // Inserts all pending elements (from last_indexed_pos to size).
  for (int j = last_indexed_pos; last_indexed_pos < size(); j = ++last_indexed_pos) {
    int i = nC - 1;
    int tmp = D[j];            // Store the pending tuple.
    for (; i >= 0 && (tmp < V[i]); i--) {  // insert by shuffling through the cracker indices C
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


template<typename T>
static T* doUpperCrack(T *lo, T *hi, T *pivot) {
  while (lo <= hi) {
    while (lo <= hi && !(*pivot < *lo)) lo++;
    while (lo <= hi &&  (*pivot < *hi)) hi--;
    if (lo < hi) std::iter_swap(lo, hi);
  }
  return lo;
}

template<typename T>
static T* doLowerCrack(T *lo, T *hi, T *pivot) {
  while (lo <= hi) {
    while (lo <= hi &&  (*lo < *pivot)) lo++;
    while (lo <= hi && !(*hi < *pivot)) hi--;
    if (lo < hi) std::iter_swap(lo, hi);
  }
  return lo;
}

template<typename T>
static T* partition(T *lo, T *hi, T *pivot, bool upperCrack) {
  assert(lo <= pivot && pivot <= hi);
  if (upperCrack) {
    std::iter_swap(lo, pivot);
    // fprintf(stderr, "upperCrack %d, %d\n", *lo, hi - lo);
    T* ret = doUpperCrack(lo + 1, hi, lo);
    std::iter_swap(ret, lo);
    return ret;
  }
  std::iter_swap(hi, pivot);
  T* ret = doLowerCrack(lo, hi - 1, hi);
  std::iter_swap(ret, hi);
  return ret;
}

// Partitions roughly in the middle satisfying the DECRACK_AT.
template<typename T, typename CMP>
static T* rough_middle_partition2(T *L, T *R, int min_gap) {
  T *lo = L;
  T *hi = R - 1;
  T *H = L + (R - L) / 2;
  assert(R - L >= min_gap * 2);
  while (lo < hi) {
    T* pivot = lo + rng.nextInt(hi - lo + 1);
    // fprintf(stderr, "lo = %d(%d), hi = %d(%d), pivot = %d(%d)\n", lo - L, *lo, hi - L, *hi, pivot - L, *pivot);
    pivot = partition(lo, hi, pivot, false);

    if (pivot - L >= min_gap) return pivot;
    if (R - pivot >= min_gap) return pivot;

    // fprintf(stderr, "pivot = %d(%d)\n", pivot - L, *pivot);
    if (pivot == lo) {
      pivot = partition(lo, hi, pivot, true);
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
template<typename T>
static T* rough_middle_partition(T *L, T *R, int min_gap) {
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
    // T pval = *p;
    if (lo >= min_gap && N - hi >= min_gap) {
      // Normal Crack
      ret = partition(arr + lo, arr + hi, p, false) - arr;
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
      ilo += !(arr[i] < *p);
      ihi +=  (arr[i] < *p);
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

int LeafBucket::get_piece_by_value(int &v, int &L, int &R) {
  flush_pending_inserts();
  // assert(check(D[0],false,D[0],false));
  int i = 0;
  while (i < nC && (v >= V[i])) i++;  // Find the cracker indices that covers v
  L = (i == 0) ? 0 : C[i - 1];          // Left crack boundary.
  R = (i == nC) ? size() : C[i];        // Right crack boundary.
  while (R - L > CRACK_AT) {            // Narrow down the piece using DDR.
    int M = rough_middle_partition(D + L + (i ? 1 : 0), D + R, DECRACK_AT) - D;
    // assert(abs((L + R) / 2 - M) <= 1);
    add_cracker_index(i, M);
    // fprintf(stderr,"CRACKING %d %d, [%d %d]\n", M, D[M], L, R);
    if ((v < D[M])) R = M; else L = M, i++;  // Adjust the cracker index i.
  }
  assert(i >= 0 && i <= nC);
  // assert(check(D[0],false,D[0],false));
  return i;
}

// BUCKET_TPLC(int)::get_piece_by_index(int idx, int &L, int &R) {
//   assert(idx >= 0 && idx < size());        // index must fall within the bucket's size
//   flush_pending_inserts(cmp);
//   int i = 0;
//   while (i<nC && idx >= C[i]) i++;    // find the cracker indices that covers the idx
//   L = i==0? 0 : C[i-1];          // the left crack
//   R = i==nC? size() : C[i];          // the right crack
//   while (R-L > CRACK_AT){          // narrow down the piece using DDR
//     int M = rough_middle_partition(D + L + (i ? 1 : 0), D + R, DECRACK_AT) - D;
//     // assert(abs((L + R) / 2 - M) <= 1);
//     add_cracker_index(i,M);
//     if (idx < M) R=M; else L=M, i++;  // adjust the cracker index
//   }
//   assert(i>=0 && i<=nC);
//   // assert(check(D[0], false, D[0], false));
//   return i;
// }

int LeafBucket::crack(int &v, int &i, int &L, int &R, bool sort_piece) {
  assert(!next);                     // It doesn't make sense crack a chained bucket!
  // assert(check(D[0], false, D[0], false));
  i = get_piece_by_value(v, L, R);    // Find the piece [L,R) containing v.
  assert(L >= 0 && L <= R && R <= size());
  // assert(check(D[0], false, D[0], false));
  if (!piece_is_sorted(i)) {
    if (sort_piece) {                       // Sort the piece if requested.
      std::sort(D + L, D + R);
      piece_set_sorted(i, true);
    } else {
      for (int at = L; at < R; at++)
        if (D[at] == v) return at;
      return R;
    }
  }
  // assert(check(D[0], false, D[0], false));
  int pos = L;
  while (pos < R && D[pos] < v) pos++;
  // debug(0); fprintf(stderr, "pos = %d, v = %d\n", pos, v);
  return pos;
  return std::lower_bound(D + L, D + R, v) - D;
}


// BUCKET_TPLC(bool)::erase(int &v) {
//   int i, L, R, at = crack(v, i, L, R, false);
//   if (at >= R || !eq(D[at], v)) return false;  // The element to be erased is not found!

//   // Decrack this cracker piece (it becomes too small) or
//   // if the deleted element index is a cracker index
//   if (nC && R - L + 1 <= DECRACK_AT) {
//     remove_cracker_index((i > 0) ? --i : i);
//   } else if (i > 0 && at == L) {
//     at = L + 1;
//     for (int j = L + 2; j < R; j++)    // find a replacement element for the cracker index
//       if (cmp(D[j], D[at])) at = j;  // that is the smallest element in the piece (L,R)
//     std::swap(D[L], D[at]);
//     piece_set_sorted(i, false);
//     V[i - 1] = D[L];
//   }

//   assert(at < R && eq(D[at], v));    // the element v must be found!

//   // IMPROVE: use pending delete? antimatter?
//   piece_set_unsorted_onwards(i);      // unset the sorted bit i onwards
//   assert(i == nC || cmp(D[at], V[i]));
//   for (int j = i; j < nC; j++) {        // shuffle out the deleted element
//     R = C[j]--;
//     D[at] = D[R - 1];
//     D[R - 1] = D[R];
//     at = R;
//   }
//   D[at] = D[--size_];   // The deleted element has been shuffled out from the bucket.
//   last_indexed_pos--;   // Adjust the pending index.

//   // assert(check(D[0],false,D[0],false));
//   return true;
// }

// BUCKET_TPL(bool)::debug(const char *msg, int i, int j) const {
//   for (int k=0; k<nC; k++)
//     fprintf(stderr,"C[%d/%d] = %d, %d (sorted = %d)\n",
//       k,nC,C[k],(int)D[C[k]],piece_is_sorted(k));
//   fprintf(stderr,"%s : i=%d/N=%d, j=%d/nC=%d, D[i,i+1] = %d, %d; I=%d, N=%d, next=%d\n",
//     msg, i,size(), j,nC, (int)D[i],(int)D[i+1], last_indexed_pos,size(),next_bucket_id);
//   return false;
// }

// BUCKET_TPLC(bool)::check(T lo, bool useLo, T hi, bool useHi) const {
//   if (useLo) for (int i = 0; i < size(); i++) if (cmp(D[i], lo)) {
//     fprintf(stderr,"D[%d] = %d, lo = %d\n", i, D[i], lo);
//     return debug("useLo failed", i, 0);
//   }
//   if (useHi) for (int i = 0; i < size(); i++) if (!cmp(D[i], hi)) {
//     fprintf(stderr,"D[%d] = %d, hi = %d\n", i, D[i], hi);
//     return debug("useHi failed", i, 0);
//   }
//   for (int i=0,j=0; i<last_indexed_pos; i++){                    // check cracker indices
//     if (j<nC && C[j]==i) assert(eq(V[j],D[i])), lo = D[i], j++;
//     if (piece_is_sorted(j) && (j==nC? (i+1<last_indexed_pos) : (i<C[j])) && cmp(D[i+1], D[i]))
//       return debug("sortedness violation", i,j);
//     if (j>0 && cmp(D[i], D[C[j-1]])) return debug("lower bound fail", i,j);
//     if (j<nC && cmp(D[C[j]], D[i])) return debug("upper bound fail", i,j);
//     // if (j<nC && !(cmp(D[i], D[C[j]]))) return debug("upper bound fail", i,j);
//   }
//   return true;
// }




void mark_hi(int *D, int N, int P, int *hi, int &nhi) {
  for (int i = 0; i < N; i++) {
    hi[nhi] = i;
    nhi += D[i] >= P;
  }
}

void mark_lo(int *D, int N, int P, int *lo, int &nlo) {
  for (int i = 0; i < N; i++) {
    lo[nlo] = i;
    nlo += D[i] < P;
  }
}

void fusion(int *Lp, int *Rp, int *hi, int *lo, int &nhi, int &nlo) {
  int m = std::min(nhi, nlo); assert(m > 0);
  int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
  nhi -= m; nlo -= m;
  while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
}

LeafBucket* LeafBucket::detach_and_get_next() {
  LeafBucket* ret = next;
  next = tail = NULL;
  pending_insert = 1;
  return ret;
}

void LeafBucket::add_chain(LeafBucket *b) {
  if (!next || !tail) next = tail = b;
  else {
    tail->next = b;
    tail = b;
  }
}

void _add(LeafBucket *&b, LeafBucket *nb) {
  if (!b) b = nb;
  else b->add_chain(nb);
}

void LeafBucket::distribute_values(int pivot, LeafBucket *chain[2]) {
  // fprintf(stderr, "has left\n");
  while (N) {
    int i = !(D[--N] < pivot);
    // fprintf(stderr, "proc left %d, i = %d\n", N, i);
    chain[i]->leaf_insert(D[N]);
  }
}

LeafBucket* LeafBucket::transfer_to(LeafBucket *b, int pivot) {
  for (int i = 0; i < N; i++) {
    if (D[i] >= pivot) {
      b->leaf_insert(D[i]);
      D[i--] = D[--N];
    }
  }
  return b;
}

void LeafBucket::leaf_split(vector<pair<int, LeafBucket*>> &ret) {
  ret.clear();
  assert(next);
  assert(cap == LEAF_BSIZE); // The first bucket must be the smallest capacity.

  if (!next->next) {
    LeafBucket *b = detach_and_get_next(); b->detach_and_get_next();

    if (N + b->N <= LEAF_BSIZE) {
      for (int i = 0; i < b->N; i++)
        D[N++] = b->D[i];
    } else {
      // assert(b->cap == cap);

      // Ensure both have at least 5 elements.
      assert(N >= 5 || b->N >= 5);
      while (N < 5) D[N++] = b->D[--b->N];
      assert(N >= 5 || b->N >= 5);
      while (b->N < 5) b->D[b->N++] = D[--N];
      assert(N >= 5 && b->N >= 5);

      int R[5];
      Random rng(140384); // TODO: use randomized seed.
      for (int i = 0; i < 3; i++) {
        assert(N > 0);
        int j = rng.nextInt(N);
        R[i] = D[j];
        D[j] = D[--N];
      }
      for (int i = 0; i < 2; i++) {
        assert(b->N > 0);
        int j = rng.nextInt(b->N);
        R[i + 3] = b->D[j];
        b->D[j] = b->D[--b->N];
      }

      for (int i = 0; i < 5; i++) {
        // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      }

      std::nth_element(R, R + 2, R + 5);
      int pivot = R[2];
      D[N++] = R[0];
      D[N++] = R[1];
      D[N++] = R[3];
      b->D[b->N++] = R[4];

      LeafBucket *nb = transfer_to(new_leaf(LEAF_BSIZE), pivot);
      b->transfer_to(nb, pivot);
      for (int i = 0; i < b->N; i++) {
        leaf_insert(b->D[i]);
      }
      ret.push_back(make_pair(pivot, nb));
    }

    delete_leaf(b);

  } else {
    // fprintf(stderr, "split N = %d\n", N);
    // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
    assert(N + tail->N >= 11);
    while (N < 11) {
      D[N++] = tail->D[--tail->N];
    }
    int R[11];
    Random rng(140384); // TODO: use randomized seed.
    for (int i = 0; i < 11; i++) {
      assert(N > 0);
      int j = rng.nextInt(N);
      R[i] = D[j];
      D[j] = D[--N];
    }
    assert(N >= 0);

    LeafBucket *Nb = this;

    // Replace elements with gradually decreasing probability.
    for (int i = 1; Nb->next; i++) {
      Nb = Nb->next;
      assert(i > 0);
      int j = rng.nextInt(i);
      if (j < 11) {
        assert(Nb->N > 0);
        int k = rng.nextInt(Nb->N);
        // fprintf(stderr, "swap %d  <>  %d,   %d %d\n", R[j], Nb->D[k], j, k);
        swap(R[j], Nb->D[k]);
      }
    }
    // fprintf(stderr, "split2 N = %d\n", q.size());

    for (int i = 0; i < 11; i++) {
      // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
    }

    std::nth_element(R, R + 5, R + 11);
    int pivot = R[5];
    R[5] = R[10];
    for (int i = 0; i < 10; i++) {
      D[N++] = R[i];
    }
    D[N++] = Nb->D[--Nb->N];
    // fprintf(stderr, "queue size = %lu\n", q.size());
    // fprintf(stderr, "split3 N = %d, pivot = %d\n", next->N, pivot);

    // debug(10);

    LeafBucket *chain[2] { this, new_leaf(LEAF_BSIZE) };

    // Split the first bucket (this bucket).
    for (int i = 0; i < N; i++) {
      // fprintf(stderr, "hihi %d < %d\n", i, N);
      if (D[i] >= pivot) {
        chain[1]->leaf_insert(D[i]);
        D[i--] = D[--N];
      }
    }

    Nb = detach_and_get_next();

    LeafBucket *Lb = NULL, *Rb = NULL;
    // TODO: optimize locality.
    int hi[LEAF_BSIZE], nhi = 0;
    int lo[LEAF_BSIZE], nlo = 0;
    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        fusion(Lb->D, Rb->D, hi, lo, nhi, nlo);
        if (!nhi) { _add(chain[0], Lb); Lb = NULL; }
        if (!nlo) { _add(chain[1], Rb); Rb = NULL; }
      } else if (!Lb) {
        if (!Nb) break;
        Lb = Nb;
        Nb = Nb->detach_and_get_next();
        if (!Lb->is_full()) break;
      } else if (!nhi) {
        assert(Lb);
        mark_hi(Lb->D, Lb->N, pivot, hi, nhi);
        if (!nhi){ _add(chain[0], Lb); Lb = NULL; }
      } else if (!Rb) {
        if (!Nb) break;
        Rb = Nb;
        Nb = Nb->detach_and_get_next();
        if (!Rb->is_full()) break;
      } else if (!nlo) {
        assert(Rb);
        mark_lo(Rb->D, Rb->N, pivot, lo, nlo);
        if (!nlo){ _add(chain[1], Rb); Rb = NULL; }
      } else {
        assert(0);
      }
    }
    assert(!Nb);

    // fprintf(stderr, "splited\n");
    if (Lb) Lb->distribute_values(pivot, chain), delete_leaf(Lb);
    if (Rb) Rb->distribute_values(pivot, chain), delete_leaf(Rb);
    ret.push_back(make_pair(pivot, chain[1]));
  }
}

int LeafBucket::promote_last() {
  nth_element(D, D + N - 1, D + N);
  return D[--N];
}

void LeafBucket::leaf_optimize() {
  assert(pending_insert >= 0);
  sort(D, D + N);
  pending_insert = 0;
}

int LeafBucket::lower_pos(int value) {
  // return std::lower_bound(D, D + N, value) - D;
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

pair<bool,int> LeafBucket::leaf_lower_bound(int value) {
  int i, L, R, pos = crack(value, i, L, R, true);
  return (pos == N) ? make_pair(false, 0) : make_pair(true, D[pos]);
}

InternalBucket::InternalBucket() : LeafBucket(INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
}

InternalBucket::InternalBucket(int promotedValue, Bucket *left, Bucket *right) : LeafBucket(INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
  D[0] = promotedValue;
  C[0] = left;
  C[1] = right;
  N = 1;
}

InternalBucket::~InternalBucket() {
  nInternals--;
}

InternalBucket* InternalBucket::internal_split() {
  InternalBucket *nb = new InternalBucket();
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
    nb->C[j] = C[i + 1];
  }
  N /= 2;
  nb->N = N;
  nb->C[0] = C[N];
  // assert(check());
  return nb;
}


void InternalBucket::internal_insert(int value, Bucket *b) {
  // assert(check());
  assert(!is_full());
  int i = N - 1;
  assert(i >= 0);
  while (i >= 0 && D[i] > value) {
    D[i + 1] = D[i];
    C[i + 2] = C[i + 1];
    i--;
  }
  D[i + 1] = value;
  C[i + 2] = b;
  N++;
  // assert(check());
}

int InternalBucket::internal_promote_last() {
 return D[--N];
}

Bucket*& InternalBucket::child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}


class CTree {
  Bucket *root;

 public:

  const char *version = "Crack 2048 bucket";

  CTree() {
    root = new_leaf(LEAF_BSIZE);
  }

  void debug() {
    root->debug(0);
    // fprintf(stderr, "\n");
  }

  void insert(int value) {
    // fprintf(stderr, "ins %d\n", value);
    Bucket *b = root;
    while (!b->is_leaf()) b = ((InternalBucket*) b)->child_bucket(value);
    ((LeafBucket*) b)->leaf_insert(value);
    // root->debug(0);
  }

  bool erase(int value) {
    return true;
  }

  pair<InternalBucket*, int> parents[10];
  vector<pair<int, LeafBucket*>> nbs;
  double t1 = 0, t2 = 0, t3 = 0;

  pair<bool, int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    int i = 0;
    Bucket *b = root;
    // if (!root->is_leaf()) parents[i++] = make_pair((InternalBucket*) root, ((InternalBucket*) root)->child_pos(value));
      // fprintf(stderr, "leaf_split iii = %d\n", i);
    while (true) {
      // fprintf(stderr, "leaf_split ii = %d\n", i);
      if (!b->is_leaf()) {
        // t1 += time_it([&] {
          assert(i < 10);
          parents[i].first = (InternalBucket*) b;
          parents[i].second = ((InternalBucket*) b)->lower_pos(value);
          // fprintf(stderr, "child pos %d <= %d\n", parents[i].second, parents[i].first->size());
          b = parents[i].first->child(parents[i].second);
          i++;
        // });

        if (parents[i-1].second < parents[i-1].first->size() && parents[i-1].first->data(parents[i-1].second) == value)
          return make_pair(true, value);

      } else if (((LeafBucket*) b)->next_bucket()) {
        // t2 += time_it([&] {
          ((InternalBucket*) b)->leaf_split(nbs);
          while (!nbs.empty()) {
            if (i == 0) {
              b = root = new InternalBucket(nbs[0].first, b, nbs[0].second);
              for (int j = 1; j < (int) nbs.size(); j++) {
                ((InternalBucket*) b)->internal_insert(nbs[j].first, nbs[j].second);
              }
              break;
            }
            b = parents[--i].first;
            while (!b->is_full() && !nbs.empty()) {
              ((InternalBucket*) b)->internal_insert(nbs.back().first, nbs.back().second);
              nbs.pop_back();
            }
            if (nbs.empty()) break;
            InternalBucket *inb = ((InternalBucket*) b)->internal_split();
            int promotedValueInternal = ((InternalBucket*) b)->internal_promote_last();
            for (int j = 0; j < (int) nbs.size(); j++) {
              if (nbs[j].first >= promotedValueInternal) {
                inb->internal_insert(nbs[j].first, nbs[j].second);
              } else {
                ((InternalBucket*) b)->internal_insert(nbs[j].first, nbs[j].second);
              }
            }
            nbs.clear();
            nbs.push_back(make_pair(promotedValueInternal, inb));
          }
        // });
      } else {
        break;
      }
      // debug();
    }

    pair<bool, int> ret;
    // t3 += time_it([&] {
      ret = ((LeafBucket*) b)->leaf_lower_bound(value);
      while (!ret.first && i > 0) {
        auto &p = parents[--i];
        if (p.second < p.first->size()) {
          // fprintf(stderr, "here %d\n", p.first->data(p.second));
          ret = make_pair(true, p.first->data(p.second));
          break;
        }
      }
    // });
    // fprintf(stderr, "awww %d %d\n", ret.first, ret.second);
      // debug();
    return ret;
  }
};

}

#endif
