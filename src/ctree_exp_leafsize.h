/*
tansfer leaf
use leaf chained as separate allocator
*/

#ifndef _CTREE_H_
#define _CTREE_H_

#include <cstdio>
#include <cstring>
#include <cassert>
#include <queue>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace chrono;

int nLeaves, nInternals, nCap, nDes, locked;

namespace ctree {

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

#define INTERNAL_BSIZE       64  // Must be power of two.
#define LEAF_BSIZE           64  // Must be power of two.
#define LEAF_CHAINED_BSIZE  4098  // Must be power of two.

#define BUCKET(b) bucket_allocator.get(b)
#define CBUCKET(b) chained_bucket_allocator.get(b)
#define CHILDREN(b) (*child_allocator.get(BUCKET(b)->C))

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

template<typename T>
class Allocator {
 public:

  priority_queue<int, vector<int>, greater<int>> free_indices;
  int cap;
  T *D;
  int N;

  void init(int initial_cap) {
    D = new T[cap = initial_cap];
    N = 0;
  }

  void save(string fn) {
    FILE *out = fopen(fn.c_str(), "wb");
    fwrite(&N, sizeof(int), 1, out);
    fwrite(D, sizeof(T), N, out);
    fclose(out);
  }

  void load(string fn) {
    FILE *f = fopen(fn.c_str(), "rb");
    fread(&N, sizeof(int), 1, f);
    D = new T[cap = N];
    fread(D, sizeof(T), N, f);
    fclose(f);
  }

  int alloc() {
    if (free_indices.empty()) {
      if (N == cap) {
        fprintf(stderr, "double %d\n", cap);
        T *newD = new T[cap * 4];
        memcpy(newD, D, sizeof(T) * cap);
        cap *= 4;
        delete[] D;
        D = newD;
      }
      return N++;
    }
    int idx = free_indices.top();
    free_indices.pop();
    return idx;
  }

  void destroy(int idx) {
    free_indices.push(idx);
  }

  T* get(int idx) const {
    assert(idx != -1);
    assert(idx >= 0);
    assert(idx < N);
    return &D[idx];
  }
};

struct ChainedBucket {
  int N;
  int next;
  int D[LEAF_CHAINED_BSIZE];

  int detach_and_get_next() {
    int ret = next;
    next = -1;
    return ret;
  }

  bool is_full() const { return N == LEAF_CHAINED_BSIZE; }

  bool append(int value) {
    if (is_full()) return false;
    D[N++] = value;
    return true;
  }
};

struct Bucket {
  int C;       // Pointer to children.
  int P;       // Pending insert if positive or pending delete if negative.
  int N;       // Number of data elements in this bucket pointed by D.
  int cap;     // The capacity of this bucket pointed by D.
  int parent;  // Pointer to the parent bucket in bucket_allocator.
  int next;    // Pointer to the next chained bucket in bucket_allocator.
  int tail;    // Pointer to the last chained bucket in bucket_allocator.
  int D[LEAF_BSIZE];       // Data values.

  bool is_full() const { return N == cap; }
  bool is_leaf() const { return C == -1; }
  int slack() const { return cap - N; }

  bool append(int value) {
    assert(is_leaf());
    if (is_full()) return false;
    D[N++] = value;
    P++;
    return true;
  }

  void leaf_optimize() {
    assert(P >= 0);
    sort(D, D + N);
    P = 0;
  }

  int leaf_lower_pos(int value) {
    if (P) leaf_optimize();
    // assert(leaf_check());
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    // assert(leaf_check());
    // leaf_debug();
    // fprintf(stderr, "leaf_lower_pos = %d, pos = %d, %d %d\n", value, pos, L, R);
    return pos;
  }

  int internal_lower_pos(int value) {
    int pos = 0;
    while (pos < N && D[pos] < value) pos++;
    return pos;
  }


  int leaf_promote_first() {
    // TODO: optimize
    P = 1;
    int smallest_pos = 0;
    int pos = 1;
    while (pos < N) {
      if (D[pos] < D[smallest_pos]) smallest_pos = pos;
      pos++;
    }
    swap(D[smallest_pos], D[--N]);
    return D[N];
  }

  int internal_promote_first(int *C) {
    int ret = D[0];
    N--;
    for (int i = 0; i < N; i++) {
      D[i] = D[i + 1];
      C[i] = C[i + 1];
    }
    C[N] = C[N + 1];
   return ret;
  }

  int internal_promote_last() {
   return D[--N];
  }

  int detach_and_get_next() {
    int ret = next;
    next = tail = -1;
    return ret;
  }

  void destroy() {
    assert(cap > 0);
    nCap -= cap;
    if (is_leaf()) {
      nLeaves--;
    } else {
      nInternals--;
    }
    nDes++;
    cap = 0;
  }
};

const char *version = "57/2048";

class CTree {
  Allocator<Bucket> bucket_allocator;
  Allocator<ChainedBucket> chained_bucket_allocator;
  Allocator<int[INTERNAL_BSIZE + 1]> child_allocator;
  int root;

 public:

  CTree() {
    bucket_allocator.init(100000000 / INTERNAL_BSIZE * 2);
    chained_bucket_allocator.init(100000000 / LEAF_CHAINED_BSIZE + 1000000);
    child_allocator.init(100000000 / INTERNAL_BSIZE * 2);

    root = bucket_allocator.alloc();
    leaf_init(root, -1, LEAF_BSIZE);
  }

  void chained_leaf_init(int b) {
    CBUCKET(b)->N = 0;
    CBUCKET(b)->next = -1;
    nCap += LEAF_CHAINED_BSIZE;
    nLeaves++;
  }

  void leaf_init(int b, int parent, int cap) {
    BUCKET(b)->C = -1;
    BUCKET(b)->cap = cap;
    BUCKET(b)->parent = parent;
    BUCKET(b)->P = 0;
    BUCKET(b)->N = 0;
    BUCKET(b)->next = -1;
    BUCKET(b)->tail = -1;
    nCap += cap;
    nLeaves++;
    assert(BUCKET(b)->is_leaf());
  }

  void internal_init(int b, int parent, int cap, int left_child) {
    BUCKET(b)->C = child_allocator.alloc();
    BUCKET(b)->cap = cap;
    BUCKET(b)->parent = parent;
    BUCKET(b)->P = 0;
    BUCKET(b)->N = 0;
    BUCKET(b)->next = -1;
    BUCKET(b)->tail = -1;

    nCap += cap;
    assert(!BUCKET(b)->is_leaf());
    nInternals++;

    CHILDREN(b)[0] = left_child;
    BUCKET(left_child)->parent = b;
  }

  void delete_bucket(int b) {
    BUCKET(b)->destroy();
    bucket_allocator.destroy(b);
  }

  void delete_chained_bucket(int b) {
    chained_bucket_allocator.destroy(b);
  }

  static void mark_hi(int *D, int N, int P, int *hi, int &nhi) {
    for (int i = 0; i < N; i++) {
      hi[nhi] = i;
      nhi += D[i] >= P;
    }
  }

  static void mark_lo(int *D, int N, int P, int *lo, int &nlo) {
    for (int i = 0; i < N; i++) {
      lo[nlo] = i;
      nlo += D[i] < P;
    }
  }

  static void fusion(int *Lp, int *Rp, int *hi, int *lo, int &nhi, int &nlo) {
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  void distribute_values(int b, int pivot, int chain[2]) {
    // fprintf(stderr, "has left\n");
    while (CBUCKET(b)->N) {
      int i = !(CBUCKET(b)->D[--CBUCKET(b)->N] < pivot);
      // fprintf(stderr, "proc left %d, i = %d\n", N, i);
      leaf_insert(chain[i], CBUCKET(b)->D[CBUCKET(b)->N]);
    }
  }

/*
  int transfer_to(int src, int des, int pivot) {
    int oldN = BUCKET(src)->N;
    BUCKET(src)->N = 0;
    int arr[2] { src, des };
    for (int i = 0; i < oldN; i++) {
      int j = !(BUCKET(src)->D[i] < pivot);
      leaf_insert(arr[j], BUCKET(src)->D[i]);
    }
    return des;
  }
*/
  void leaf_split(int b, int &promotedValue, int &new_bucket) {
    // assert(leaf_check());
    // fprintf(stderr, "b = %d / %d\n", b, bucket_allocator.N);
    assert(BUCKET(b)->next != -1);
    assert(BUCKET(b)->cap == INTERNAL_BSIZE); // The first bucket must be the smallest capacity.
    new_bucket = -1;

/*
    if (CBUCKET(BUCKET(b)->next)->next == -1 && 0) {
      int nb = BUCKET(b)->detach_and_get_next(); BUCKET(nb)->detach_and_get_next();

      if (BUCKET(b)->N + BUCKET(nb)->N <= INTERNAL_BSIZE) {
        for (int i = 0; i < BUCKET(nb)->N; i++)
          BUCKET(b)->D[BUCKET(b)->N++] = BUCKET(nb)->D[i];
      } else {
        // assert(BUCKET(nb)->cap == cap);

        // Ensure both have at least 5 elements.
        assert(BUCKET(b)->N >= 5 || BUCKET(nb)->N >= 5);
        while (BUCKET(b)->N < 5) BUCKET(b)->D[BUCKET(b)->N++] = BUCKET(nb)->D[--BUCKET(nb)->N];
        assert(BUCKET(b)->N >= 5 || BUCKET(nb)->N >= 5);
        while (BUCKET(nb)->N < 5) BUCKET(nb)->D[BUCKET(nb)->N++] = BUCKET(b)->D[--BUCKET(b)->N];
        assert(BUCKET(b)->N >= 5 && BUCKET(nb)->N >= 5);

        int R[5];
        Random rng(140384); // TODO: use randomized seed.
        for (int i = 0; i < 3; i++) {
          assert(BUCKET(b)->N > 0);
          int j = rng.nextInt(BUCKET(b)->N);
          R[i] = BUCKET(b)->D[j];
          BUCKET(b)->D[j] = BUCKET(b)->D[--BUCKET(b)->N];
        }
        for (int i = 0; i < 2; i++) {
          assert(BUCKET(nb)->N > 0);
          int j = rng.nextInt(BUCKET(nb)->N);
          R[i + 3] = BUCKET(nb)->D[j];
          BUCKET(nb)->D[j] = BUCKET(nb)->D[--BUCKET(nb)->N];
        }

        for (int i = 0; i < 5; i++) {
          // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
        }

        std::nth_element(R, R + 2, R + 5);
        int pivot = R[2];
        BUCKET(b)->D[BUCKET(b)->N++] = R[0];
        BUCKET(b)->D[BUCKET(b)->N++] = R[1];
        BUCKET(b)->D[BUCKET(b)->N++] = R[3];
        BUCKET(nb)->D[BUCKET(nb)->N++] = R[4];

        int nb2 = bucket_allocator.alloc();
        leaf_init(nb2, BUCKET(b)->parent, INTERNAL_BSIZE);
        transfer_to(b, nb2, pivot);
        transfer_to(nb, nb2, pivot);
        for (int i = 0; i < BUCKET(nb)->N; i++) {
          leaf_insert(b, BUCKET(nb)->D[i]);
        }
        promotedValue = pivot;
        new_bucket = nb2;
      }

      delete_bucket(nb);

    } else {
*/
      // fprintf(stderr, "split N = %d\n", BUCKET(b)->N);
      // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
      assert(BUCKET(b)->N + CBUCKET(BUCKET(b)->tail)->N >= 11);
      while (BUCKET(b)->N < 11) {
        BUCKET(b)->D[BUCKET(b)->N++] = CBUCKET(BUCKET(b)->tail)->D[--CBUCKET(BUCKET(b)->tail)->N];
      }
      // fprintf(stderr, "tail N = %d\n", BUCKET(b)->N);
      int R[11];
      Random rng(140384); // TODO: use randomized seed.
      for (int i = 0; i < 11; i++) {
        assert(BUCKET(b)->N > 0);
        int j = rng.nextInt(BUCKET(b)->N);
        R[i] = BUCKET(b)->D[j];
        BUCKET(b)->D[j] = BUCKET(b)->D[--BUCKET(b)->N];
      }
      assert(BUCKET(b)->N >= 0);

      // fprintf(stderr, "random N = %d\n", BUCKET(b)->N);

      int Nb = BUCKET(b)->next;

      // Replace elements with gradually decreasing probability.
      for (int i = 1; Nb != -1; i++) {
        assert(i > 0);
        int j = rng.nextInt(i);
        if (j < 11) {
          assert(CBUCKET(Nb)->N > 0);
          int k = rng.nextInt(CBUCKET(Nb)->N);
          // fprintf(stderr, "swap %d  <>  %d,   %d %d\n", R[j], Nb->D[k], j, k);
          swap(R[j], CBUCKET(Nb)->D[k]);
        }
        if (CBUCKET(Nb)->next == -1) break;
        Nb = CBUCKET(Nb)->next;
      }

      for (int i = 0; i < 11; i++) {
        // fprintf(stderr, "R[%d] = %d\n", i, R[i]);
      }

      std::nth_element(R, R + 5, R + 11);
      int pivot = R[5];
      R[5] = R[10];
      for (int i = 0; i < 10; i++) {
        BUCKET(b)->D[BUCKET(b)->N++] = R[i];
      }
      BUCKET(b)->D[BUCKET(b)->N++] = CBUCKET(Nb)->D[--CBUCKET(Nb)->N];
      // fprintf(stderr, "pivot = %d\n", pivot);
      // fprintf(stderr, "split3 N = %d, pivot = %d\n", next->N, pivot);

      // debug(10);

      new_bucket = bucket_allocator.alloc();
      leaf_init(new_bucket, BUCKET(b)->parent, INTERNAL_BSIZE);
      int chain[2] { b, new_bucket };

      // Split the first bucket (this bucket).
      for (int i = 0; i < BUCKET(b)->N; i++) {
        // fprintf(stderr, "hihi %d < %d\n", i, BUCKET(b)->N);
        if (BUCKET(b)->D[i] >= pivot) {
          leaf_insert(chain[1], BUCKET(b)->D[i]);
          BUCKET(b)->D[i--] = BUCKET(b)->D[--BUCKET(b)->N];
        }
      }

      // fprintf(stderr, "fision\n");
      Nb = BUCKET(b)->detach_and_get_next();

      int Lb = -1, Rb = -1;
      // TODO: optimize locality.
      int hi[LEAF_CHAINED_BSIZE], nhi = 0;
      int lo[LEAF_CHAINED_BSIZE], nlo = 0;
      while (true) {
        // fprintf(stderr, "Lb = %d, Rb = %d\n", Lb, Rb);
        if (nhi && nlo) {
          // fprintf(stderr, "a\n");
          assert(Lb != -1 && Rb != -1);
          fusion(CBUCKET(Lb)->D, CBUCKET(Rb)->D, hi, lo, nhi, nlo);
          if (!nhi) { add_chain(chain[0], Lb); Lb = -1; }
          if (!nlo) { add_chain(chain[1], Rb); Rb = -1; }
        } else if (Lb == -1) {
          // fprintf(stderr, "b\n");
          if (Nb == -1) break;
          Lb = Nb;
          Nb = CBUCKET(Nb)->detach_and_get_next();
          if (!CBUCKET(Lb)->is_full()) break;
        } else if (!nhi) {
          assert(Lb != -1);
          // fprintf(stderr, "c %d\n", CBUCKET(Lb)->N);
          mark_hi(CBUCKET(Lb)->D, CBUCKET(Lb)->N, pivot, hi, nhi);
          if (!nhi){ add_chain(chain[0], Lb); Lb = -1; }
        } else if (Rb == -1) {
          // fprintf(stderr, "d\n");
          if (Nb == -1) break;
          Rb = Nb;
          Nb = CBUCKET(Nb)->detach_and_get_next();
          if (!CBUCKET(Rb)->is_full()) break;
        } else if (!nlo) {
          // fprintf(stderr, "e\n");
          assert(Rb != -1);
          mark_lo(CBUCKET(Rb)->D, CBUCKET(Rb)->N, pivot, lo, nlo);
          if (!nlo){ add_chain(chain[1], Rb); Rb = -1; }
        } else {
          assert(0);
        }
      }
      assert(Nb == -1);

      // fprintf(stderr, "splited\n");
      if (Lb != -1) distribute_values(Lb, pivot, chain), delete_chained_bucket(Lb);
      if (Rb != -1) distribute_values(Rb, pivot, chain), delete_chained_bucket(Rb);
      promotedValue = pivot;
    // }
    // assert(leaf_check());
  }

  void internal_insert(int b, int value, int nb, int left = 0) {
    // fprintf(stderr, "ii %d\n", b);
    int *C = CHILDREN(b);
    assert(!BUCKET(b)->is_leaf());
    assert(!BUCKET(b)->is_full());
    int i = BUCKET(b)->N - 1;
    while (i >= 0 && BUCKET(b)->D[i] > value) {
      BUCKET(b)->D[i + 1] = BUCKET(b)->D[i];
      C[i + 2] = C[i + 1];
      i--;
    }
    BUCKET(b)->D[i + 1] = value;
    if (left == -1) {
      C[i + 2] = C[i + 1];
      C[i + 1] = nb;
    } else {
      C[i + 2] = nb;
    }
    BUCKET(b)->N++;
    // fprintf(stderr, "set parent %d  -> %d\n", nb, b);
    BUCKET(nb)->parent = b;
    // assert(check());
  }

  void internal_erase(int b, int *C, int pos, int stride) {
    BUCKET(b)->N--;
    while (pos < BUCKET(b)->N) {
      BUCKET(b)->D[pos] = BUCKET(b)->D[pos + 1];
      C[pos + stride] = C[pos + stride + 1];
      pos++;
    }
    if (!stride) C[pos] = C[pos + 1];
  }

/*
  int leaf_shift_left(int b, int pos) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(BUCKET(L)->is_leaf());
    assert(BUCKET(R)->is_leaf());
    if (BUCKET(L)->next != -1) return false;
    if (BUCKET(R)->next != -1) return false;

    // Move from M to L as many as possible.
    int changed = 0;
    while (!BUCKET(L)->is_full() && BUCKET(R)->N) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      BUCKET(b)->D[pos] = BUCKET(R)->leaf_promote_first();
      changed = 1;
    }
    if (!BUCKET(L)->is_full() && !BUCKET(R)->N) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      internal_erase(b, CHILDREN(b), pos, 1);
      delete_bucket(R);
      return 2; // R is empty, compaction is done.
    }
    return changed;
  }

  bool internal_shift_left(int b, int pos, int numMove = INTERNAL_BSIZE + 1) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(!BUCKET(L)->is_leaf());
    assert(!BUCKET(R)->is_leaf());
    assert(BUCKET(L)->next == -1);
    assert(BUCKET(R)->next == -1);

    // Move from R to L as many as possible.
    bool changed = false;
    while (!BUCKET(L)->is_full() && BUCKET(R)->N && numMove-- > 0) {
      internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
      BUCKET(b)->D[pos] = BUCKET(R)->internal_promote_first(CHILDREN(R));
      changed = true;
    }
    // if (!BUCKET(L)->is_full() && !BUCKET(R)->N) {
    //   internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
    //   internal_erase(b, CHILDREN(b), pos, 1);
    //   delete_bucket(R);
    //   return true; // R is empty, compaction is done.
    // }
    return changed;
  }

  bool internal_shift_right(int b, int pos, int numMove = INTERNAL_BSIZE + 1) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(!BUCKET(L)->is_leaf());
    assert(!BUCKET(R)->is_leaf());
    assert(BUCKET(L)->next == -1);
    assert(BUCKET(R)->next == -1);

    // Move from L to R as many as possible.
    bool changed = false;
    while (!BUCKET(R)->is_full() && BUCKET(L)->N && numMove-- > 0) {
      internal_insert(R, BUCKET(b)->D[pos], child(L, BUCKET(L)->N), -1);
      BUCKET(b)->D[pos] = BUCKET(L)->internal_promote_last();
      changed = true;
    }
    // if (!BUCKET(R)->is_full());
    // R->internal_insert(ib->data(rpos), L->child(L->size()), -1);
    // ib->internal_erase(rpos, 0);
    // delete L;
    return changed;
  }
*/
  int internal_split(int b) {
    int nb = bucket_allocator.alloc();
    int *C = CHILDREN(b);
    internal_init(nb, BUCKET(b)->parent, LEAF_BSIZE, C[BUCKET(b)->N / 2]);
    int *nbC = CHILDREN(nb);
    for (int i = BUCKET(b)->N / 2, j = 0; i < BUCKET(b)->N; i++) {
      BUCKET(nb)->D[j++] = BUCKET(b)->D[i];
      nbC[j] = C[i + 1];
      BUCKET(nbC[j])->parent = nb;
    }
    BUCKET(b)->N /= 2;
    BUCKET(nb)->N = BUCKET(b)->N;
    // assert(check());
    return nb;
  }

  int internal_find_child_pos(int b, int c) {
    int *C = CHILDREN(b);
    for (int i = 0; i <= BUCKET(b)->N; i++) {
      if (C[i] == c) return i;
    }
    assert(0);
    return 0;
  }

  bool split_chain(int b) {
    // assert(check());
    if (BUCKET(b)->next == -1) return false;
    assert(!locked);

    int promotedValue;
    int nb;
    // fprintf(stderr, "split_chain %d, %d\n", BUCKET(b)->N, BUCKET(b)->next);
    // assert(check());
    leaf_split(b, promotedValue, nb);
    // assert(check());
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);

    int parent = BUCKET(b)->parent;
    while (parent != -1 && nb != -1) {
      if (BUCKET(parent)->is_full()) {
        // fprintf(stderr, "parful\n");

        // Optional optimization:
        
        // assert(!BUCKET(parent)->is_leaf());
        // int pp = BUCKET(parent)->parent;
        // if (pp != -1) {
        //   assert(!BUCKET(pp)->is_leaf());
        //   int pos = internal_find_child_pos(pp, parent);
        //   // assert(check());

        //   if (pos > 0 && BUCKET(parent)->D[0] < promotedValue && internal_shift_left(pp, pos - 1, 1)) {
        //     // fprintf(stderr, "shift left %d, p = %d, b = %d, root = %d, pos = %d / %d\n", pp, parent, b, root, pos, BUCKET(pp)->N);
        //     assert(!BUCKET(parent)->is_full());
        //     // assert(check());
        //     internal_insert(parent, promotedValue, nb);
        //     nb = -1;
        //     // assert(check());
        //     break;
        //   } else if (pos < BUCKET(pp)->N && promotedValue < BUCKET(parent)->D[BUCKET(parent)->N - 1] && internal_shift_right(pp, pos, 1)) {
        //     // assert(check());
        //     assert(!BUCKET(parent)->is_full());
        //     internal_insert(parent, promotedValue, nb);
        //     nb = -1;
        //     // assert(check());
        //     break;
        //   }
        // }
        
        int inb = internal_split(parent);
        int promotedValueInternal = BUCKET(parent)->internal_promote_last();
        if (promotedValue >= promotedValueInternal) {
          internal_insert(inb, promotedValue, nb);
        } else {
          internal_insert(parent, promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = BUCKET(nb)->parent;
      } else {
        // fprintf(stderr, "internal\n");
        internal_insert(parent, promotedValue, nb);
        nb = -1;
      }
    }
    // assert(check());
    // fprintf(stderr, "nb = %d\n", nb);
    if (nb != -1) {
      // Replace root
      // fprintf(stderr, "OLD ROOT %d\n", root);
      assert(parent == -1);
      assert(BUCKET(root)->parent == -1);
      int new_root = bucket_allocator.alloc();
      internal_init(new_root, -1, INTERNAL_BSIZE, root); // New internal bucket.
      internal_insert(new_root, promotedValue, nb);
      root = new_root;
      // fprintf(stderr, "NEW ROOT %d\n", root);
    }
    // fprintf(stderr, "\n\n\ndone split %d\n", b);
    // debug();
    // assert(check());
    return true;
  }

  bool optimize(int b = -1) {
    if (b == -1) {
      if (root == -1) return false;
      while (optimize(root)) {
        if (BUCKET(root)->N == 0) {
          fprintf(stderr, "PROMOTE ROOT\n");
          assert(!BUCKET(root)->is_leaf());
          int c = CHILDREN(root)[0];
          delete_bucket(root);
          root = c;
          BUCKET(root)->parent = -1;
        }
        fprintf(stderr, "optimize\n");
      }
      return false;
    }

    if (BUCKET(b)->is_leaf()) {
      // fprintf(stderr, "split_chain %d\n", b);
      bool splitted = split_chain(b);
      // fprintf(stderr, "split_chain %d\n", b);
      while (split_chain(b));
      // fprintf(stderr, "split_chain %d %d\n", b, splitted);
      return splitted;
    }

    bool changed = false;
    // fprintf(stderr, "internal %d\n", b);
    for (int i = 0; i <= BUCKET(b)->N; i++) {
      if (optimize(CHILDREN(b)[i])) {
        i = -1;
        changed = true;
      }
    }

    // OPTIONAL:
    // for (int i = 0; i < BUCKET(b)->N; i++) {
    //   int L = child(b, i);
    //   int R = child(b, i + 1);
    //   assert(BUCKET(L)->is_leaf() == BUCKET(R)->is_leaf());
    //   if (BUCKET(L)->is_leaf()) {
    //     if (leaf_shift_left(b, i)) changed = 1;
    //   } else {
    //     if (internal_shift_left(b, i)) changed = 1;
    //   }
    // }

    // fprintf(stderr, "internal %d %d\n", b, changed);
    return changed;
  }

  void debug_data(int b) {
    if (BUCKET(b)->is_leaf()) {
      fprintf(stderr, "N = %d (p=%d, LEAF), ", BUCKET(b)->N, BUCKET(b)->P);
    } else {
      fprintf(stderr, "N = %d, ", BUCKET(b)->N);
    }
    fprintf(stderr, "[");
    for (int i = 0; i < BUCKET(b)->N; i++) {
      fprintf(stderr, "%d ", BUCKET(b)->D[i]);
    }
    fprintf(stderr, "]");
  }

  int debug(int b = -1, int depth = 0) {
    for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
    if (b == -1) b = root;
    int sz = 0;
    if (BUCKET(b)->is_leaf()) {
      while (true) {
        sz += BUCKET(b)->N;
        debug_data(b);
        if (BUCKET(b)->next == -1) break;
        fprintf(stderr, "------ ");
        b = BUCKET(b)->next;
      }
      fprintf(stderr, "\n");
    } else {
      debug_data(b);
      fprintf(stderr, "\n");
      for (int i = 0; i <= BUCKET(b)->N; i++) {
        sz += debug(CHILDREN(b)[i], depth + 1);
      }
      for (int i = 0; i < depth; i++) fprintf(stderr, "  ");
      fprintf(stderr, "size = %d\n", sz);
    }
    return sz;
  }

  int max_depth(int b = -1) {
    if (b == -1) b = root;
    if (BUCKET(b)->is_leaf()) return 1;
    int ret = -1;
    for (int i = 0; i <= BUCKET(b)->N; i++) {
      int d = max_depth(CHILDREN(b)[i]);
      assert(ret == -1 || ret == d);
      ret = d;
    }
    return ret + 1;
  }

  int slack(int b = -1, int last = 0) {
    if (b == -1) b = root;
    int ret = BUCKET(b)->slack();
    // if (ret > 10) fprintf(stderr, "slack = %d, for leaf = %d, last = %d\n", ret, b->is_leaf(), last);
    if (BUCKET(b)->is_leaf()) return ret; else ret = 0;
    for (int i = 0; i <= BUCKET(b)->N; i++) {
      ret += slack(CHILDREN(b)[i], i >= BUCKET(b)->N);
    }
    return ret;
  }

  int size() {
    int ret = 0;
    int b = root;
    while (b != -1) {
      Bucket *B = BUCKET(b);
      ret += B->N;
      b = B->next;
    }
    return ret;
  }

/*
  bool leaf_compact(int b, int start, int end) {
    int last = 0;
    for (int i = start; i < end; i++) {
      last = leaf_shift_left(b, i);
    }
    assert(last == 2);
    fprintf(stderr, "saved 1 leaf %d\n", b);
    return true;
  }

  bool leaf_compact(int b) {
    int slack = 0, start = 1;
    Bucket *B = BUCKET(CHILDREN(b)[0]);
    if (B->next == -1) slack += B->slack(), start = 0;
    for (int i = 1; i <= BUCKET(b)->N; i++) {
      B = BUCKET(CHILDREN(b)[i]);
      if (B->next != -1) {
        slack += B->slack();
        if (slack >= INTERNAL_BSIZE)
          return leaf_compact(b, start, i);
      } else {
        return false;
        // if (slack) fprintf(stderr, "gathered slack = %d, %d\n", slack, i - start);
        slack = 0;
        start = i + 1;
      }
    }
    fprintf(stderr, "gathered slack = %d\n", slack);
    return false;
  }
*/
  pair<int, int> find_bucket(int value, bool include_internal) {
    int b = root;
    while (true) {
      // fprintf(stderr, "find_bucket %d\n", b);
      if (BUCKET(b)->is_leaf()) {
        if (!split_chain(b)) break;
        b = BUCKET(b)->parent;
        assert(b != -1);
      } else {
        int pos = BUCKET(b)->internal_lower_pos(value);
        if (include_internal && pos < BUCKET(b)->N && BUCKET(b)->D[pos] == value) {
          return make_pair(b, pos); // Found in the internal bucket.
        }
        b = CHILDREN(b)[pos];    // Search the child.
      }
    }
    return make_pair(b, 0);
  }

  pair<bool, int> lower_bound(int value) {
    // assert(check());
    // fprintf(stderr, "lower_bound %d\n", value);
    pair<int, int> p = find_bucket(value, true);
    // fprintf(stderr, "lower_bound1 %d\n", value);

    // Found in internal bucket.
    pair<bool, int> ret = make_pair(false, 0);
    if (!BUCKET(p.first)->is_leaf()) {
      ret = make_pair(true, value);
    } else {
      int pos = BUCKET(p.first)->leaf_lower_pos(value);
      if (pos < BUCKET(p.first)->N) {
        ret = make_pair(true, BUCKET(p.first)->D[pos]);

        // OPTIONAL optimization:
        // int parent = BUCKET(p.first)->parent;
        // if (parent != -1) {
        //   leaf_compact(parent);
        // }
      } else {
        int b = BUCKET(p.first)->parent;
        while (b != -1) {
          pos = BUCKET(b)->internal_lower_pos(value);
          if (pos < BUCKET(b)->N) {
            ret = make_pair(true, BUCKET(b)->D[pos]);
            // fprintf(stderr, "x");
            break;
          }
          b = BUCKET(b)->parent;
        }
      }
    }
    // fprintf(stderr, "lower_bound4 %d\n", value);

    return ret;
  }


  pair<bool, int> lower_bound2(int value) {
    // fprintf(stderr, "lower_bound %d\n", value);
    return lower_bound_rec(root, value);
  }

  pair<bool, int> lower_bound_rec(int b, int value) {
    if (BUCKET(b)->is_leaf()) {
      if (split_chain(b)) {
        b = BUCKET(b)->parent;
      } else {
        int pos = BUCKET(b)->leaf_lower_pos(value);
        return  (pos < BUCKET(b)->N) ? make_pair(true, BUCKET(b)->D[pos]) : make_pair(false, 0);
      }
    }
    int pos = BUCKET(b)->internal_lower_pos(value);
    if (pos < BUCKET(b)->N && BUCKET(b)->D[pos] == value) {
      return make_pair(true, value); // Found in the internal bucket.
    }
    auto res = lower_bound_rec(CHILDREN(b)[pos], value);
    if (!res.first && pos < BUCKET(b)->N) {
      res = make_pair(true, BUCKET(b)->D[pos]);
    }
    return res;
  }

  void add_chain(int head, int next) {
    // assert(BUCKET(next)->is_leaf());
    assert(head != -1);
    Bucket *b = BUCKET(head);
    assert(b->is_leaf());
    if (b->next == -1) {
      b->next = b->tail = next;
    } else {
      CBUCKET(b->tail)->next = next;
      b->tail = next;
    }
  }

  void batch_insert(int *arr, int N) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0;
    while (i + LEAF_CHAINED_BSIZE <= N) {
      int idx = chained_bucket_allocator.alloc();
      chained_leaf_init(idx);
      ChainedBucket *b = CBUCKET(idx);
      for (int j = 0; j < LEAF_CHAINED_BSIZE; j++) {
        b->D[b->N++] = arr[i++];
      }
      // fprintf(stderr, "chain %d\n", i);
      add_chain(root, idx);
    }
    // fprintf(stderr, "done %d\n", i);
    while (i < N) {
      insert(arr[i++]);
    }
    // fprintf(stderr, "done %d\n", i);
  }

  void insert(int value) {
    leaf_insert(root, value);
  }

  void insert2(int &b, int value) {
    // fprintf(stderr, "ins %d\n", value);
    // if (value == 711)  debug();
    assert(b != -1);
    while (true) {
      Bucket *B = BUCKET(b);
      if (B->is_leaf()) break;
      // int *C = *child_allocator.get(B->C);
      // b = B->child(value, child_allocator);
      assert(0);
    }

    leaf_insert(b, value);
    // root->debug(0);
  }

  void leaf_insert(int b, int value) {
    if (BUCKET(b)->append(value)) return;
    // assert(BUCKET(b)->cap == INTERNAL_BSIZE);
    // assert(BUCKET(BUCKET(b)->tail)->next == -1);
    int tail = BUCKET(b)->tail;
    if (tail == -1 || CBUCKET(tail)->is_full()) {
      tail = chained_bucket_allocator.alloc();
      chained_leaf_init(tail);
      CBUCKET(tail)->append(value);
      add_chain(b, tail);
    } else {
      CBUCKET(tail)->append(value);
    }
  }

  bool erase(int value) {
    return true;
  }

  bool leaf_debug(const char *msg, int i, int j) {
    fprintf(stderr, "%s, i = %d, j = %d, root = %d\n", msg, i, j, root);
    return false;
  }

  bool leaf_check(int b, int lo, bool useLo, int hi, bool useHi) {
    if (useLo) for (int i = 0; i < BUCKET(b)->N; i++) if ((BUCKET(b)->D[i] < lo)) {
      fprintf(stderr,"D[%d] = %d, lo = %d\n", i, BUCKET(b)->D[i], lo);
      return leaf_debug("useLo failed", i, 0);
    }
    if (useHi) for (int i = 0; i < BUCKET(b)->N; i++) if ((BUCKET(b)->D[i] > hi)) {
      fprintf(stderr,"D[%d] = %d, hi = %d\n", i, BUCKET(b)->D[i], hi);
      return leaf_debug("useHi failed", i, 0);
    }
    return true;
  }

  bool check(int b, int lo, int hi) {
    return BUCKET(b)->is_leaf() ? leaf_check(b, lo, true, hi, true) : check(b, lo);
  }

  bool check(int b = -1, int lo = -2147483648) {
    if (b == -1) b = root;
    // fprintf(stderr, "check %d, leaf = %d\n", b, BUCKET(b)->is_leaf());
    assert(BUCKET(b)->N >=0 && BUCKET(b)->N <= BUCKET(b)->cap);
    if (BUCKET(b)->is_leaf()) return leaf_check(b, lo, true, 0, false);
    if (BUCKET(CHILDREN(b)[0])->parent != b) return leaf_debug("parent mismatch", BUCKET(CHILDREN(b)[0])->parent, b);
    if (BUCKET(b)->N && !check(CHILDREN(b)[0], lo, BUCKET(b)->D[0])) return false;
    for (int i = 0; i < BUCKET(b)->N; i++) {
      // fprintf(stderr, "%d cek anak ke-%d\n", b, i);
      assert(i == 0 || BUCKET(b)->D[i - 1] <= BUCKET(b)->D[i]);
      if (BUCKET(CHILDREN(b)[i + 1])->parent != b) return leaf_debug("parent mismatch internal", BUCKET(CHILDREN(b)[i + 1])->parent, b);
      if (!check(CHILDREN(b)[i + 1], BUCKET(b)->D[i], (i + 1 < BUCKET(b)->N) ? BUCKET(b)->D[i + 1] : 2147483647)) return false;
    }
    return true;
  }

  void save(string fn) {
    child_allocator.save(fn + ".child");
    bucket_allocator.save(fn + ".bucket");
    fprintf(stderr, "chained size = %d, free = %lu\n", chained_bucket_allocator.N, chained_bucket_allocator.free_indices.size());
  }

  void load(string fn) {
    child_allocator.load(fn + ".child");
    bucket_allocator.load(fn + ".bucket");
  }
};


/*

pair<bool,int> leaf_erase_largest() {
  assert(!next);
  if (N == 0) return make_pair(false, D[0]);
  int pos = 0;
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

bool leaf_erase(int &v) {
  // assert(leaf_check());
  int pos = 0;
  while (pos < N) {
    if (D[pos] == v) {
      D[pos] = D[--N];
      return true;
    }
  }
  return false;
}



int leaf_promote_last() {
  if (P > 0) leaf_optimize();
  return D[--N];
}

Bucket*& child_bucket(int value) {
  int pos = 0;
  while (pos < N && !(value < D[pos])) pos++;
  return child(pos);
}


class CTree {
  Bucket *root;

 public:

  const char *version = "Exp LEAF_BSIZE 2048";

  double t1 = 0, t2 = 0, t3 = 0;

  bool compact_leaves(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    while (!L->is_full() && M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->set_data(lpos, M->leaf_promote_first());
    }
    if (!L->is_full() && !M->size()) {
      L->leaf_insert(ib->data(lpos));
      ib->internal_erase(lpos, 1);
      delete_bucket(M);
      return true; // M is empty, compaction is done.
    }

    Bucket *R = (Bucket*) ib->child(rpos + 1);
    assert(!R->next_bucket());

    // Move from M to R as many as possible.
    while (M->size()) {
      assert(!R->is_full());
      R->leaf_insert(ib->data(rpos));
      ib->set_data(rpos, M->leaf_promote_last());
    }
    assert(!R->is_full());
    R->leaf_insert(ib->data(rpos));
    ib->internal_erase(rpos, 0);
    delete_bucket(M);
    return true;
  }

  bool compact_internals(Bucket *ib, int rpos) {
    int lpos = rpos - 1;
    Bucket *L = (Bucket*) ib->child(lpos);
    Bucket *M = (Bucket*) ib->child(rpos);
    assert(!L->next_bucket() && !M->next_bucket());

    // Move from M to L as many as possible.
    while (!L->is_full() && M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->set_data(lpos, M->internal_promote_first());
    }
    if (!L->is_full() && !M->size()) {
      L->internal_insert(ib->data(lpos), M->child(0));
      ib->internal_erase(lpos, 1);
      delete M;
      return true; // M is empty, compaction is done.
    }

    Bucket *R = (Bucket*) ib->child(rpos + 1);
    assert(!R->next_bucket());

    // Move from M to R as many as possible.
    while (M->size()) {
      assert(!R->is_full());
      R->internal_insert(ib->data(rpos), M->child(M->size()), -1);
      ib->set_data(rpos, M->internal_promote_last());
    }
    assert(!R->is_full());
    R->internal_insert(ib->data(rpos), M->child(M->size()), -1);
    ib->internal_erase(rpos, 0);
    delete M;
    return true;
  }



  pair<bool, int> erase_largest(Bucket *b) {
    assert(b->is_leaf());
    auto res = ((Bucket*) b)->leaf_erase_largest();
    if (res.first) return res;
    while (b) {
      assert(!b->size());
      Bucket *parent = (Bucket*) b->get_parent();
      if (b->is_leaf()) {
        delete_bucket((Bucket*) b);
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
      return ((Bucket*) p.first)->leaf_erase(value); // May return false.
    }

    Bucket* upper = find_bucket(value, false).first;

    // The buckets may have moved.
    p = find_bucket(value, true);

    Bucket *ib = (Bucket*) p.first;
    int pos = p.second;
    assert(ib->data(pos) == value);

    auto res = erase_largest(upper);
    // fprintf(stderr, "ii res = %d, largest = %d\n", res.first, res.second);
    if (!res.first) {
      ib->internal_erase(pos, 0);
    } else {
      ib->set_data(pos, res.second);
    }
    // assert(check());
    return true;
  }
};
*/
}


#endif
