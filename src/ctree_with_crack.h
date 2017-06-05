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
  int *D;
  int pending_insert; // -1 for InternalBucket, >= 0 for LeafBucket.
  int N, cap;
  Bucket *parent;

 public:
  int get_cap() const { return cap; }
  int size() const { return N; }
  bool is_full() const { return size() == cap; }
  bool is_leaf() const { return pending_insert >= 0; }
  int data(int i) { assert(i >= 0 && i < N); return D[i]; }
  void set_data(int pos, int value) { D[pos] = value; }
  Bucket* get_parent() const { return parent; }
  void set_parent(Bucket *parent) { this->parent = parent; }
  void optimize();
  int debug(int depth);
  bool check(int lo) const;
  bool check(int lo, int hi) const;
};

class LeafBucket : public Bucket {
 protected:
  unsigned char nC;       // Number of cracker indices.
  unsigned long long S1;  // Sorted bits of each cracker index.
  int C[MAX_INDEX];       // Cracker index positions.
  int V[MAX_INDEX];       // Cached cracker values.
  LeafBucket *next, *tail;  // Store pending inserts in a linked list.

 public:
  ~LeafBucket();
  LeafBucket(Bucket *parent, int cap);
  void init(int cap);
  LeafBucket* next_bucket() { return next; }

  void clear_indexes() { S1 = /* S2 = S3 = S4 = */ nC = 0; pending_insert = 0; }
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
  bool leaf_erase(int &v);
  pair<bool,int> leaf_erase_largest();

  bool leaf_debug(const char *msg, int i, int j) const;

  // call this function to check the consistency of this Bucket structure
  bool leaf_check() const;
  bool leaf_check(int lo, bool useLo, int hi, bool useHi) const;

  void leaf_debug();
  void leaf_insert(int v);
  void leaf_split(int &promotedValue, LeafBucket *&nb);
  void leaf_optimize();
  int leaf_promote_last();
  int leaf_lower_pos(int value);
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
  InternalBucket(Bucket *parent, Bucket *left_child);
  Bucket*& child(int i) { assert(i >= 0 && i <= N); return C[i]; }
  int child_pos(int value);
  Bucket*& child_bucket(int value);
  int internal_lower_pos(int value);
  InternalBucket* internal_split();
  void internal_insert(int value, Bucket *b);
  int internal_promote_last();
  void internal_erase(int pos);
};


vector<LeafBucket*> free_leaves[30];

LeafBucket* new_leaf(Bucket *parent, int cap) {
  for (int i = 2; ; i++) {
    if ((1 << i) == cap) {
      if (free_leaves[i].empty()) {
        free_leaves[i].push_back(new LeafBucket(parent, cap));
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
    ((LeafBucket*) b)->leaf_debug();
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

bool Bucket::check(int lo, int hi) const {
  if (is_leaf()) return ((LeafBucket*) this)->leaf_check(lo, true, hi, true);
  return check(lo);
}

bool Bucket::check(int lo) const {
  assert(N >=0 && N <= cap);
  assert(is_leaf() || ((InternalBucket*) this)->child(N));
  if (is_leaf()) return ((LeafBucket*) this)->leaf_check();
  InternalBucket *ib = (InternalBucket*) this;
  if (N) {
    if (ib->child(0)->get_parent() != ib) return false;
    ib->child(0)->check(lo, D[0]);
  }
  for (int i = 0; i < N; i++) {
    assert(ib->child(i));
    if (i > 0) assert(ib->data(i - 1) < ib->data(i));
    if (ib->child(i + 1)->get_parent() != ib) return false;
    if (!ib->child(i + 1)->check(D[i], (i + 1 < N) ? D[i + 1] : 2147483647))
      return false;
  }
  return true;
}


LeafBucket::LeafBucket(Bucket *parent, int cap) {
  this->parent = parent;
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

void LeafBucket::leaf_debug() {
  fprintf(stderr, "N = %d (p=%d, c=%d, LEAF), ", N, pending_insert, nC);
  for (int i = 0; i < N; i++) {
    fprintf(stderr, "%d ", D[i]);
  }
  for (int i = 0; i < nC; i++) {
    fprintf(stderr, "C[%d] = (pos=%d, val=%d), ", i, C[i], V[i]);
  }
}

void LeafBucket::leaf_insert(int value) {
  // assert(leaf_check());
  // assert(is_leaf());
  // assert(N >= 0);
  if (!is_full()) {
    D[N++] = value;
    pending_insert++;
  } else {
    if (!tail || tail->is_full()) {
      // assert(cap == LEAF_BSIZE);
      // add_chain(new_leaf(parent, LEAF_BSIZE));
    // } else if (tail->is_full()) {
      add_chain(new_leaf(parent, LEAF_BSIZE));
    }
    tail->D[tail->N++] = value;
  }
  // assert(leaf_check());
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
  // fprintf(stderr, "%d <= %d\n", pending_insert, size());
  // debug(10);
  assert(pending_insert <= size());

  // No index yet, all the elements are considered "inserted".
  if (!nC) { pending_insert = 0; S1 = 0; return; } // S2 = S3 = S4 =

  // Flushing only makes sense when the bucket is not chained.
  assert(!next);

  // IMPROVE: bulk insert? (Currently using Merge Completely).
  int minC = nC;
  // Inserts all pending elements (from pending_insert to size).
  for (int j = size() - pending_insert; pending_insert; j = size() - --pending_insert) {
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
T* partition(T *F, T *L, T const &v){
  while (true){
    while (true)
      if (F == L) return F;
      else if ((*F < v)) ++F;
      else break;
    --L;
    while (true)
      if (F == L) return F;
      else if (!bool((*L < v))) --L;
      else break;
    std::iter_swap(F, L);
    ++F;
  }
  assert(false);
}

// partitions roughly in the middle satisfying the DECRACK_AT
template<typename T>
T* rough_middle_partition(T *L, T *R) {
  assert(R-L >= CRACK_AT);
  int ntry = 10;
  for (T *i=L, *j=R, *p; ntry--; ){
    std::iter_swap(j-1, i + rng.nextInt(j-i));
    std::iter_swap(j-1, p = partition(i, j-1, *(j-1)));
    if (p-L <= DECRACK_AT) i = p;
    else if (R-p <= DECRACK_AT) j = p;
    else return p;
  }
  T *M = L+((R-L)>>1);
  std::nth_element(L, M, R);
  return M;
}

int LeafBucket::get_piece_by_value(int &v, int &L, int &R) {
  flush_pending_inserts();
  // assert(check(D[0],false,D[0],false));
  int i = 0;
  while (i < nC && (v >= V[i])) i++;  // Find the cracker indices that covers v
  L = (i == 0) ? 0 : C[i - 1];          // Left crack boundary.
  R = (i == nC) ? size() : C[i];        // Right crack boundary.
  while (R - L > CRACK_AT) {            // Narrow down the piece using DDR.
    int M = rough_middle_partition(D + L + (i ? 1 : 0), D + R) - D;
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
  // assert(leaf_check());
  assert(!next);                     // It doesn't make sense crack a chained bucket!
  // assert(check(D[0], false, D[0], false));
  i = get_piece_by_value(v, L, R);    // Find the piece [L,R) containing v.
  // fprintf(stderr, "i = %d, v = %d, %d %d, N = %d, pi = %d\n", i,v,L,R,N, pending_insert);
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

pair<bool,int> LeafBucket::leaf_erase_largest() {
  assert(!next);
  flush_pending_inserts();
  // assert(leaf_check());
  // fprintf(stderr, "leaf1");
  if (N == 0) return make_pair(false, D[0]);
  // fprintf(stderr, "leaf2");
  int pos = 0;
  if (nC) {
    pos = C[nC - 1] + 1;
    if (N - pos < DECRACK_AT) {
      remove_cracker_index(nC - 1);
    }
    // TODO: optimize.
    piece_set_sorted(nC, false);
  } else {
    pending_insert = 0;
  }
  int largest_pos = pos++;
  while (pos < N) {
    if (D[pos] > D[largest_pos])
      largest_pos = pos;
    pos++;
  }
  // fprintf(stderr, "pos %d %d\n", pos, largest_pos);
  swap(D[largest_pos], D[--N]);
  // assert(leaf_check());
  return make_pair(true, D[N]);
}

bool LeafBucket::leaf_erase(int &v) {
  // assert(leaf_check());
  int i, L, R, at = crack(v, i, L, R, false);
  // fprintf(stderr, "at = %d < %d\n", at, R);
  if (at >= R || D[at] != v) return false;  // The element to be erased is not found!

  // Decrack this cracker piece (it becomes too small) or
  // if the deleted element index is a cracker index
  if (nC && R - L <= DECRACK_AT) {
    remove_cracker_index((i > 0) ? --i : i);
  } else if (i > 0 && at == L) {
    at = L + 1;
    for (int j = L + 2; j < R; j++)    // find a replacement element for the cracker index
      if (D[j] < D[at]) at = j;  // that is the smallest element in the piece (L,R)
    std::swap(D[L], D[at]);
    piece_set_sorted(i, false);
    V[i - 1] = D[L];
  }

  assert(at < R && D[at] == v);    // the element v must be found!

  // IMPROVE: use pending delete? antimatter?
  piece_set_unsorted_onwards(i);      // unset the sorted bit i onwards
  assert(i == nC || D[at] < V[i]);
  for (int j = i; j < nC; j++) {        // shuffle out the deleted element
    R = C[j]--;
    D[at] = D[R - 1];
    D[R - 1] = D[R];
    at = R;
  }
  D[at] = D[--N];   // The deleted element has been shuffled out from the bucket.

  // assert(leaf_check());
  return true;
}

bool LeafBucket::leaf_debug(const char *msg, int i, int j) const {
  for (int k=0; k < nC; k++)
    fprintf(stderr,"C[%d/%d] = %d, %d (sorted = %d)\n",
      k, nC, C[k], (int)D[C[k]], piece_is_sorted(k));
  fprintf(stderr,"%s : i=%d/N=%d, j=%d/nC=%d, D[i,i+1] = %d, %d; I=%d, N=%d, next=%p\n",
    msg, i,size(), j,nC, (int)D[i],(int)D[i+1], pending_insert,size(), next);
  return false;
}

bool LeafBucket::leaf_check() const {
  return leaf_check(D[0],false,D[0],false);
}

bool LeafBucket::leaf_check(int lo, bool useLo, int hi, bool useHi) const {
  if (useLo) for (int i = 0; i < size(); i++) if ((D[i] < lo)) {
    fprintf(stderr,"D[%d] = %d, lo = %d\n", i, D[i], lo);
    return leaf_debug("useLo failed", i, 0);
  }
  if (useHi) for (int i = 0; i < size(); i++) if ((D[i] > hi)) {
    fprintf(stderr,"D[%d] = %d, hi = %d\n", i, D[i], hi);
    return leaf_debug("useHi failed", i, 0);
  }
  for (int i=0,j=0; i<N - pending_insert; i++){                    // check cracker indices
    if (j<nC && C[j]==i) assert((V[j] == D[i])), lo = D[i], j++;
    if (piece_is_sorted(j) && (j==nC? (i+1<N - pending_insert) : (i<C[j])) && (D[i+1] < D[i]))
      return leaf_debug("sortedness violation", i,j);
    if (j>0 && (D[i] < D[C[j-1]])) return leaf_debug("lower bound fail", i,j);
    if (j<nC && (D[C[j]] < D[i])) return leaf_debug("upper bound fail", i,j);
    // if (j<nC && !((D[i] < D[C[j]]))) return debug("upper bound fail", i,j);
  }
  return true;
}




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
  clear_indexes();
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
  int oldN = N;
  N = 0;
  LeafBucket *arr[2] { this, b };
  for (int i = 0; i < oldN; i++) {
    int j = !(D[i] < pivot);
    arr[j]->leaf_insert(D[i]);
  }
  return b;
}

void LeafBucket::leaf_split(int &promotedValue, LeafBucket *&new_bucket) {
  // assert(leaf_check());
  assert(next);
  assert(cap == LEAF_BSIZE); // The first bucket must be the smallest capacity.
  new_bucket = NULL;

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

      LeafBucket *nb = transfer_to(new_leaf(parent, LEAF_BSIZE), pivot);
      b->transfer_to(nb, pivot);
      for (int i = 0; i < b->N; i++) {
        leaf_insert(b->D[i]);
      }
      promotedValue = pivot;
      new_bucket = nb;
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

    LeafBucket *chain[2] { this, new_leaf(parent, LEAF_BSIZE) };

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
    promotedValue = pivot;
    new_bucket = chain[1];
  }
  // assert(leaf_check());
}

int LeafBucket::leaf_promote_last() {
  nth_element(D, D + N - 1, D + N);
  return D[--N];
}

void LeafBucket::leaf_optimize() {
  assert(pending_insert >= 0);
  sort(D, D + N);
  pending_insert = 0;
}

int LeafBucket::leaf_lower_pos(int value) {
  // assert(leaf_check());
  int i, L, R, pos = crack(value, i, L, R, true);
  // assert(leaf_check());
  // leaf_debug();
  // fprintf(stderr, "leaf_lower_pos = %d, pos = %d, %d %d\n", value, pos, L, R);
  return pos;
}

InternalBucket::InternalBucket(Bucket *parent, Bucket *left_child) : LeafBucket(parent, INTERNAL_BSIZE) {
  nInternals++;
  pending_insert = -1;
  C[0] = left_child;
  left_child->set_parent(this);
}

InternalBucket::~InternalBucket() {
  nInternals--;
}

InternalBucket* InternalBucket::internal_split() {
  InternalBucket *nb = new InternalBucket(parent, C[N / 2]);
  for (int i = N / 2, j = 0; i < N; i++) {
    nb->D[j++] = D[i];
    nb->C[j] = C[i + 1];
    nb->C[j]->set_parent(nb);
  }
  N /= 2;
  nb->N = N;
  // assert(check());
  return nb;
}


void InternalBucket::internal_insert(int value, Bucket *b) {
  // assert(check());
  assert(!is_full());
  int i = N - 1;
  while (i >= 0 && D[i] > value) {
    D[i + 1] = D[i];
    C[i + 2] = C[i + 1];
    i--;
  }
  D[i + 1] = value;
  C[i + 2] = b;
  N++;
  b->set_parent(this);
  // assert(check());
}

int InternalBucket::internal_lower_pos(int value) {
  int pos = 0;
  while (pos < N && D[pos] < value) pos++;
  return pos;
}

int InternalBucket::internal_promote_last() {
 return D[--N];
}

Bucket*& InternalBucket::child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}

void InternalBucket::internal_erase(int pos) {
  N--;
  while (pos < N) {
    D[pos] = D[pos + 1];
    C[pos] = C[pos + 1];
    pos++;
  }
  C[pos] = C[pos + 1];
}


class CTree {
  Bucket *root;

 public:

  const char *version = "Crack 2048";

  CTree() {
    root = new_leaf(NULL, LEAF_BSIZE);
  }

  void debug() {
    root->debug(0);
    // fprintf(stderr, "\n");
  }

  int max_depth(Bucket *b = NULL) {
    if (!b) b = root;
    if (b->is_leaf()) return 1;
    int ret = -1;
    for (int i = 0; i <= b->size(); i++) {
      int d = max_depth(((InternalBucket*) b)->child(i));
      assert(ret == -1 || ret == d);
      ret = d;
    }
    return ret + 1;
  }

  int slack(Bucket *b = NULL) {
    if (!b) b = root;
    int ret = b->get_cap() - b->size();
    if (b->is_leaf()) return ret;
    for (int i = 0; i <= b->size(); i++) {
      ret += slack(((InternalBucket*) b)->child(i));
    }
    return ret;
  }

  bool split_chain(LeafBucket *b) {
    // fprintf(stderr, "split_chain %d, %d\n", b->size(), b->next_bucket());
    if (!b->next_bucket()) return false;

    int promotedValue;
    LeafBucket *nb;
    b->leaf_split(promotedValue, nb);
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);
    InternalBucket *parent = (InternalBucket*) b->get_parent();

    while (parent && nb) {
      if (parent->is_full()) {
        // fprintf(stderr, "parful\n");
        InternalBucket *inb = parent->internal_split();
        int promotedValueInternal = parent->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          inb->internal_insert(promotedValue, nb);
        } else {
          parent->internal_insert(promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = (InternalBucket*) nb->get_parent();
      } else {
        // fprintf(stderr, "internal\n");
        parent->internal_insert(promotedValue, nb);
        nb = NULL;
        break;
      }
    }
    if (nb) {
      // Replace root
        // fprintf(stderr, "replace root\n");
      assert(parent == NULL);
      root = new InternalBucket(NULL, root);
      ((InternalBucket*) root)->internal_insert(promotedValue, nb);
    }
    return true;
  }

  pair<Bucket*, int> find_bucket(int value, bool include_internal) {
    Bucket *b = root;
    while (true) {
      if (b->is_leaf()) {
        LeafBucket *Lb = (LeafBucket*) b;
        if (!split_chain(Lb)) break;
        assert(Lb->get_parent());
        b = Lb->get_parent();
      } else {
        InternalBucket *ib = (InternalBucket*) b;
        int pos = ib->internal_lower_pos(value);
        if (include_internal && pos < ib->size() && ib->data(pos) == value) {
          return make_pair(ib, pos); // Found in the internal bucket.
        }
        b = ib->child(pos);    // Search the child.
      }
    }
    return make_pair(b, 0);
  }

  pair<bool, int> lower_bound(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    pair<Bucket*, int> p = find_bucket(value, true);

    // Found in internal bucket.
    if (!p.first->is_leaf()) return make_pair(true, value);

    LeafBucket *b = (LeafBucket*) p.first;
    int pos = b->leaf_lower_pos(value);
    if (pos < b->size()) return make_pair(true, b->data(pos));

    InternalBucket *ib = (InternalBucket*) b->get_parent();
    while (ib) {
      pos = ib->internal_lower_pos(value);
      if (pos < ib->size()) return make_pair(true, ib->data(pos));
      ib = (InternalBucket*) ib->get_parent();
    }
    return make_pair(false, 0);
  }

  void insert(int value) {
    // if (t1 > 0) fprintf(stderr, "ins %d\n", value);
    // if (value == 711)  debug();
    Bucket *b = root;
    while (!b->is_leaf()) b = ((InternalBucket*) b)->child_bucket(value);
    ((LeafBucket*) b)->leaf_insert(value);
    // root->debug(0);
  }

  pair<bool, int> erase_largest(Bucket *b) {
    assert(b->is_leaf());
    auto res = ((LeafBucket*) b)->leaf_erase_largest();
    if (res.first) return res;
    while (b) {
      assert(!b->size());
      InternalBucket *parent = (InternalBucket*) b->get_parent();
      if (b->is_leaf()) {
        delete_leaf((LeafBucket*) b);
      } else {
        delete b;
      }
      if (parent->size() > 0) return make_pair(true, parent->internal_promote_last());
      b = parent;
    }
    return res;
  }

  bool erase(int value) {
    // assert(check());
    // fprintf(stderr, "ERASE %d\n", value);
    // debug();

    pair<Bucket*, int> p = find_bucket(value, true);
    if (p.first->is_leaf()) {
      return ((LeafBucket*) p.first)->leaf_erase(value); // May return false.
    }

    Bucket* upper = find_bucket(value, false).first;

    // The buckets may have moved.
    p = find_bucket(value, true);

    InternalBucket *ib = (InternalBucket*) p.first;
    int pos = p.second;
    assert(ib->data(pos) == value);

    auto res = erase_largest(upper);
    // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
    if (!res.first) {
      ib->internal_erase(pos);
    } else {
      ib->set_data(pos, res.second);
    }
    // assert(check());
    return true;

  }

  bool check() {
    return root->check(-2147483648);
  }
};
}

#endif
