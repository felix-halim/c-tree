#include <map>
#include <set>
#include <vector>

#ifndef REP
#define REP(i,n) for (int i=0,_n=n; i<_n; i++)
#endif

#ifndef FORE
#define FORE(i,a) for (__typeof(a.begin()) i=a.begin(); i!=a.end(); i++)
#endif

using namespace std;

#ifndef _SCRACK_RANDOM_H_
#define _SCRACK_RANDOM_H_

#include <assert.h>

class Random {
    static const long long multiplier = 0x5DEECE66DLL;
    static const long long addend = 0xBLL;
    static const long long mask = (1LL << 48) - 1;
  unsigned long long seed;

public:
  Random(){ setSeed(140384); }
  Random(long long s){ setSeed(s); }
  void setSeed(long long s){ seed = (s ^ multiplier) & mask; }
  int next(int bits) {
    seed = (seed * multiplier + addend) & mask;
    return (int)(seed >> (48 - bits));
  }
  int nextInt(int n) {
    assert(n > 0);
    // i.e., n is a power of 2
    if ((n & -n) == n) return (int)((n * (long long)next(31)) >> 31);
    int bits, val;
    do { bits = next(31); val = bits % n; }
    while (bits - val + (n-1) < 0);
    assert(val >= 0 && val < n);
    return val;
  }
  int nextInt(){ return next(32); }
  long long nextLong(){ return ((long long)(next(32)) << 32) + next(32); }
  bool nextBoolean(){ return next(1) != 0; }
  float nextFloat(){ return next(24) / ((float)(1 << 24)); }
  double nextDouble(){ return (((long long)(next(26)) << 27) + next(27)) / (double)(1LL << 53); }
};

#endif
struct CIndex {
  int pos;      // the cracker index position
  int holes;    // the number of holes in front
  bool sorted;  // is this piece sorted
  int prev_pos() const { return pos - holes; }
};

typedef int value_type;                   // the data element value type (can be int or long long)
typedef map<value_type, CIndex> ci_type;  // cracker[cracker_value] = (cracker_index, cracker_holes)
typedef ci_type::iterator ci_iter;        // the iterator type for the cracker

multiset<int> pins, pdel; // pending updates (insert / delete)
int *arr, N;              // the dataset array
int *marr, msize;         // for materialization
ci_type ci;               // the cracker index

int partition(value_type *arr, value_type v, int L, int R){
  return std::partition(arr+L, arr+R, bind2nd(less<value_type>(), v)) - arr;
}

// return an iterator at R
ci_iter find_piece(ci_type &ci, int N, value_type v, int &L, int &R){
  L = 0; R = N;
  ci_iter it = ci.lower_bound(v);

  // handling off by one search
  if (it == ci.end()){
    if (it != ci.begin()){
      it--;
      L = it->second.pos;
      it++;
    }
  } else if (it == ci.begin()){
    if (v < it->first){
      R = it->second.prev_pos();
    } else {
      L = it->second.pos;
      it++;
      if (it != ci.end()) R = it->second.prev_pos();
    }
  } else if (v < it->first){
    R = it->second.prev_pos();
    it--;
    L = it->second.pos;
    it++;
  } else {
    L = it->second.pos;
    it++;
    if (it != ci.end()) R = it->second.prev_pos();
  }
  return it;
}

int add_crack(ci_type &ci, int &N, value_type v, int p){
  if (p==0 || p >= N) return p;
  ci_iter i = ci.lower_bound(v);
  if (i != ci.end()){
    if (i->second.pos == p) return p;
    ci_iter j = i;
    if (j->first == v) j++;
    if (j != ci.end())
      if (j->second.prev_pos() == p)
        return p;
  }
  if (i != ci.begin()){
    i--;
    if (i->second.pos == p) return p;
  }
  if (ci.count(v)){
    assert(ci[v].pos == p);
    return p;
  }
//  fprintf(stderr, "add crack v = %d, p = %d\n", v,p);
  ci[v] = (CIndex){ p,0, false };
  return p;
}

void print_crackers(ci_type &ci){
  FORE(i,ci) fprintf(stderr,"ch[%d] = %d, %d (%d)\n",
    i->first,i->second.pos,i->second.holes,i->second.sorted);
}

void print_all(ci_type &ci, value_type *arr, int &N, multiset<int> &pending_insert, multiset<int> &pending_delete, int line){
  fprintf(stderr,"PRINT ALL %d\n",line);
  REP(i,N) fprintf(stderr,"arr[%d] = %d\n",i,arr[i]);
  print_crackers(ci);
//  FORE(i,pending_delete) fprintf(stderr,"pending_delete = %d\n",*i);
//  FORE(i,pending_insert) fprintf(stderr,"pending_insert = %d\n",*i);
  fprintf(stderr,"----\n");
}

void check(ci_type &ci, value_type *arr, int &N, multiset<int> &pending_insert, multiset<int> &pending_delete){
//  fprintf(stderr,"CHECKING\n");
  int idx = 0, lo = 0;
  assert(pending_insert.empty() || *pending_insert.begin() != -1);
  assert(pending_delete.empty() || *pending_delete.begin() != -1);
  for (ci_iter it=ci.begin(); it!=ci.end(); it++){
    ci_iter it2 = it; it2++;
    if (it2 != ci.end()){
      if (!(it->second.pos != it2->second.prev_pos())){
        //print_all(881);
        fprintf(stderr,"at %d, %d %d\n",it->second.pos,it2->second.pos,it2->second.holes);
      }
      assert(it->second.pos != it2->second.prev_pos());
    }
    if (!(it->second.pos < N))
      fprintf(stderr,"%d %d %d, N = %d\n",it->first,it->second.pos,it->second.holes,N);
    assert(it->second.pos < N);
    REP(j,it->second.prev_pos() - idx){
      if (it->second.sorted && j && arr[idx-1]!=-1 && arr[idx]!=-1){
        if (!(arr[idx-1] < arr[idx])){
          print_all(ci,arr,N,pending_insert,pending_delete,99);
          fprintf(stderr,"idx = %d",idx);
        }
        assert(arr[idx-1] < arr[idx]);
      }
      if (arr[idx] < lo){
        fprintf(stderr,"arr[%d] < lo:: %d < %d\n",idx,arr[idx],lo);
//        print_all();
      }
      assert(arr[idx] >= lo);
      if (arr[idx] >= it->first){
        fprintf(stderr,"idx=%d, %d != %d\n",idx,arr[idx],it->first);
//        print_all();
      }
      assert(arr[idx] < it->first);
      idx++;
    }
    REP(j,it->second.holes){
      assert(arr[it->second.pos - j - 1] == -1);
    }
    lo = it->first;
    idx = it->second.pos;
//    fprintf(stderr,"lo = %d, idx = %d\n",lo,idx);
  }
  while (idx < N){
    if (arr[idx] < lo) fprintf(stderr,"arr[%d] = %d  <  lo = %d\n",idx,arr[idx],lo);
    assert(arr[idx] >= lo);
    idx++;
  }
}

// two way split in one scan
void split_ab(value_type *arr, int L, int R, value_type a, value_type b, int &i1, int &i2){
  i1 = L; R--;
  while (L<=R){
    while (L<=R && arr[L] < b){
      if (arr[L] < a) swap(arr[L], arr[i1++]);
      L++;
    }
    while (L<=R && arr[R] >= b) R--;
    if (L<R) swap(arr[L], arr[R]);
  }
  i2 = L;
}

// symmetric variant 11112222333333333333333333333        1111 222 3123323133233331 222 33
void symmetric_crack3(value_type *b, int L, int R, value_type lo, value_type hi, int &ll, int &hh){
  int lh = ll = L, hl = hh = R-1;

  while (lh < hl) {
    if (ll == lh) {
      while (lh < hl && b[lh] <  lo) lh++;  // advance ll below lo
      ll = lh;
    } else {
      while (lh < hl && b[lh] <  lo) swap(b[ll++], b[lh++]);    // drag lo
    }
    if (hh == hl) {
      while (lh < hl && b[hl] >= hi) hl--;  // advance hh above hi
      hh = hl;
    } else {
      while (lh < hl && b[hl] >= hi) swap(b[hh--], b[hl--]);    // drag hi
    }
    while (lh < hl && b[lh] >= lo && b[lh] <  hi) lh++; // advance lh
    while (lh < hl && b[hl] >= lo && b[hl] <  hi) hl--; // advance hl
    while (lh < hl && b[lh] >= hi && b[hl] <  lo) {
      swap(b[ll++], b[hl--]);
      swap(b[hh--], b[lh++]);
    }
    while (lh < hl && hl < hh && b[lh] >= hi) swap(b[hh--], b[lh++]);
    while (lh < hl && ll < lh && b[hl] <  lo) swap(b[ll++], b[hl--]);
  }
}


// three way cracks
void split_abc(value_type *arr, int L, int R,
    value_type a, value_type b, value_type c, int &i1, int &i2, int &i3){
  assert(a<=b && b<=c);
  i1 = i2 = L; R--;
  while (L<=R){
    while (L<=R && arr[L] < c){
      if (arr[L] < b){
        swap(arr[L], arr[i2++]);
        if (arr[i2-1] < a) swap(arr[i2-1], arr[i1++]);
      }
      L++;
    }
    while (L<=R && arr[R] >= c) R--;
    if (L<R) swap(arr[L], arr[R]);
  }
  i3 = L;
}

void materialize_piece(value_type *arr, value_type a, value_type b, int L, int R, value_type *marr, int &msize){
  for (int i=L; i<R; i++)
    if (arr[i] >= a && arr[i] < b)
      marr[msize++] = arr[i];
}

// materialize the query result leveraging existing cracks and don't add any cracks
void materialize_it(ci_type &ci, value_type *arr, int &N, value_type a, value_type b, int &i1, int &i2, value_type *marr, int &msize){
  int L1,R1; find_piece(ci, N, a, L1,R1);
  int L2,R2; find_piece(ci, N, b, L2,R2);

  // materialize for v1's piece
  materialize_piece(arr, a, b, L1, R1, marr, msize);

  if (L1!=L2){
    // v1 and v2 are on different pieces materialize for v2's piece
    materialize_piece(arr, a, b, L2,R2, marr, msize);
  }

  i1 = R1;  // make the intermediate pieces as views
  i2 = L2;
}

/*
int fifty4(Query &q){ return (nth % 4 == 0)? rcmat(q) : crack(q); }
int fifty8(Query &q){ return (nth % 8 == 0)? rcmat(q) : crack(q); }
int fifty16(Query &q){ return (nth % 16 == 0)? rcmat(q) : crack(q); }
int fifty32(Query &q){ return (nth % 32 == 0)? rcmat(q) : crack(q); }
*/


#include "hash.h"

void get_pendings(ci_type &ci, value_type a, value_type b,
          ci_iter it1, ci_iter it2,
          multiset<value_type> &pend,
          vector<value_type> &res){

  multiset<value_type>::iterator pit;

  if (it1 == ci.end() || a < it1->first){
    if (it1 != ci.begin()){
      it1--;
      pit = pend.lower_bound(it1->first);
      it1++;
    } else {
      pit = pend.begin();
    }
  } else {
    pit = pend.lower_bound(it1->first);
  }

  assert(it2 == ci.end() || b < it2->first);

  while (pit != pend.end()){
    if (it2 != ci.end() && *pit >= it2->first) break;
    if (*pit >= b) break; // HERE
    if (*pit < a){ pit++; continue; } // HERE
    res.push_back(*pit);
    pend.erase(pit++);
  }
}

bool piece_is_empty(ci_type &ci, ci_iter &it2){
  if (it2 == ci.end()) return false;
  ci_iter it3 = it2; it3++;
  if (it3 == ci.end()) return false;
  int R2 = it3->second.prev_pos();

  if (it2->second.pos < R2) return false;
  assert(it2->second.pos == R2);
  int holes = it2->second.holes;
  //fprintf(stderr,"ERASE %d\n",R2);
  ci.erase(it2++);
  if (it2 != ci.end()){
    it2->second.holes += holes;
    ci_iter it3 = it2; it3++;
  }
  return true;
}

void do_merge_ripple_insert(ci_type &ci, value_type *arr, int &N, multiset<int> &pending_insert, multiset<int> &pending_delete, value_type a, value_type b, value_type &new_hi){
//  fprintf(stderr,"domri\n");
//  check(ci,arr,N,pending_insert,pending_delete);
  int L1,R1; ci_iter it1 = find_piece(ci, N, a, L1,R1);
  int L2,R2; ci_iter it2 = find_piece(ci, N, b, L2,R2);

  vector<value_type> pIns;
  int rippled = 0, sepak = 0;
  get_pendings(ci,a,b,it1,it2,pending_insert,pIns);
  rippled += pIns.size();
  if (pIns.size()) new_hi = pIns.back();
//  REP(i,pIns.size()) fprintf(stderr,"(+%d) ",pIns[i]); fprintf(stderr,"\n");

  while (pIns.size()){
    int idx = -1;
    if (it2 != ci.end()){
      if (it2->second.holes){
        idx = it2->second.pos - (it2->second.holes--);
        it2->second.sorted = false;
        arr[idx] = pIns.back();
        pIns.pop_back();
      } else {
        assert(it2 != ci.end());
        if (piece_is_empty(ci,it2)) continue;
        value_type temp = arr[idx = it2->second.pos];
        arr[it2->second.pos++] = pIns.back();
        pIns.pop_back();
        it2->second.sorted = false;
        ci_iter it3 = it2; it3++;  // mark this piece as not sorted
        if (it3 != ci.end()) it3->second.sorted = false;

//        do_insert(q,temp,NULL); // insert to pending insert
        if (pending_delete.count(temp)){
          pending_delete.erase(pending_delete.lower_bound(temp));
        } else {
          pending_insert.insert(temp);
        }

        if (it2->second.pos >= N) ci.erase(it2++);
        piece_is_empty(ci,it2);
      }
    } else {
      arr[idx = N++] = pIns.back();
      pIns.pop_back();
    }
    // MCI the arr[idx]
    ci_iter i = it2, pi = it2;
    if (i != ci.begin()){
      int nth = 0, R = N;
      if (i != ci.end()){
        R = i->second.prev_pos();
        i->second.sorted = false;
      }
      for (i--; ; pi=i--,nth++){
//        fprintf(stderr,"%d. ripple arr[%d] = %d, i->first = %d (%d %d), R = %d, %d\n",
//          nth,idx,arr[idx],i->first,i->second.pos,i->second.holes,R,i->first <= arr[idx]);
        if (i->first <= arr[idx]) break; // done

        #ifndef MAX_CRACKERS
          #define MAX_CRACKERS 2000000000
        #endif
        if (ci.size() > MAX_CRACKERS){
          // delete cracker indexes that gets in the way of shuffling
          if (pi == ci.end()){
            assert(N==R);
            assert(idx==N-1);
            if (i->second.holes){
              assert(i->second.pos < N);
              swap(arr[idx = i->second.prev_pos()], arr[--N]);
              i->second.holes--;
            }
            while (i->second.holes && i->second.pos < N){
              swap(arr[i->second.prev_pos()], arr[--N]);
              i->second.holes--;
            }
            N -= i->second.holes;
            if (idx != N-1)
              swap(arr[idx], arr[N-1]), idx = N-1;
            ci.erase(i++);
          } else {
            assert(idx == pi->second.prev_pos()-1);
            if (i->second.holes){
              assert(i->second.pos < pi->second.prev_pos());
              swap(arr[idx = i->second.prev_pos()], arr[pi->second.prev_pos()-1]);
              i->second.holes--;
              pi->second.holes++;
            }
            while (i->second.holes && i->second.pos < pi->second.prev_pos()){
              swap(arr[i->second.prev_pos()], arr[pi->second.prev_pos()-1]);
              i->second.holes--;
              pi->second.holes++;
            }
            pi->second.holes += i->second.holes;
            assert(idx < pi->second.prev_pos());
            if (idx != pi->second.prev_pos()-1)
              swap(arr[idx], arr[pi->second.prev_pos()-1]), idx = pi->second.prev_pos()-1;
            ci.erase(i++);
          }
        } else {
          i->second.sorted = false;
          rippled++;
          swap(arr[i->second.pos], arr[idx]);
          if (i->second.holes){
            idx = i->second.prev_pos();
            swap(arr[i->second.pos], arr[idx]);
          } else {
            idx = i->second.pos;
          }
          i->second.pos++;
          assert(i->second.pos < N);
          if (i->second.pos == R){
            ci_iter j = i; j++;
            if (j != ci.end()) j->second.holes += i->second.holes;
            ci.erase(i++);
          }
        }
        if (i == ci.begin()) break;
        if (i == ci.end()) R = N;
        else R = i->second.prev_pos();
      }
      // check();
    }
  }
  ci_iter i = it2;
  if (i!=ci.end()) i++;
  while (i!=ci.end()){
    if (it2->second.pos < i->second.pos) break;
    ci.erase(i++);
    sepak++;
  }
//  if (rippled > 0 || sepak > 0) fprintf(stderr,"rippled = %d, sepak = %d\n",rippled,sepak);
  assert(pIns.size() == 0);
//  fprintf(stderr,"/domri\n");
//  check(ci,arr,N,pending_insert,pending_delete);
}

void do_merge_ripple_delete(ci_type &ci, value_type *arr, int &N, multiset<int> &pending_insert, multiset<int> &pending_delete, value_type a, value_type b){
//  check(ci,arr,N,pending_insert,pending_delete);
  int L1,R1; ci_iter it1 = find_piece(ci, N, a, L1,R1);
  int L2,R2; ci_iter it2 = find_piece(ci, N, b, L2,R2);

  // fetch all the pending insert and delete within range [ arr[L1], arr[R2] )
  vector<value_type> pDel;
  get_pendings(ci,a,b,it1,it2,pending_delete,pDel);

  Hash h(pDel.size()); // prepare for Hash JOIN

  // MRD all tuples in [ arr[L1], arr[R2] )
  int deli = 0;
  while (it1 != ci.end() && deli < int(pDel.size())){
    if (pDel[0] < it1->first){ // proceed if there is something to merge
      int pL1 = L1;

      // hash the pending deletes for this piece
      int toDelete = 0;
      while (deli < int(pDel.size()) && pDel[deli] < it1->first){
        int v = h.get(pDel[deli]);
        h.set(pDel[deli], (v>0)? (v+1) : 1);
        toDelete++;
        deli++;
      }

      // scan the tuples in the piece, deleting matched elements
      int R = it1->second.prev_pos();
      while (L1 < R && toDelete){
        assert(arr[L1]>=0);
        int v = h.get(arr[L1]);
        if (v > 0){
          h.set(arr[L1], v-1);
          arr[L1] = arr[--R];
          arr[R] = -1;
          it1->second.holes++;  // increase hole size
          it1->second.sorted = false;
          toDelete--;
        } else {
          L1++;
        }
      }

      if (it1 != ci.begin()){  // check if the deleted piece becomes empty
        if (it1->second.prev_pos() == pL1){  // cracker is empty
          ci_iter prev = it1;
          prev--;
          it1->second.holes += prev->second.holes; // transfer holes
          it1->second.sorted = false;
          ci.erase(prev++);
          it1 = prev;
        }
      }
      assert(toDelete == 0);
    }
    L1 = (it1++)->second.pos;
  }

  if (deli < int(pDel.size())){
    // hash the pending deletes for this piece
    int toDelete = 0;
    while (deli < int(pDel.size())){
      int v = h.get(pDel[deli]);
      h.set(pDel[deli], (v>0)? (v+1) : 1);
      toDelete++;
      deli++;
    }

    // scan the tuples in the piece, deleting matched elements
    while (L1 < N && toDelete){
      if (arr[L1] <=0){
        REP(i,N) fprintf(stderr,"%d %d\n",i,arr[i]);
        fprintf(stderr,"L1 = %d\n",L1);
      }
      assert(arr[L1]>0);
      int v = h.get(arr[L1]);
      if (v > 0){
        h.set(arr[L1], v-1);
        arr[L1] = arr[--N];
        arr[N] = -1;
        toDelete--;

        while (!ci.empty()){
          ci_iter li = ci.end(); li--;
          if (li->second.pos < N) break;
          assert(li->second.pos == N);
          N -= li->second.holes;
          ci.erase(li);
        }
      } else {
        L1++;
      }
    }
    assert(toDelete == 0);
  }
  assert(int(pDel.size()) == deli);
//  fprintf(stderr,"DEL CHECKING... "); multiset<int> asdf; check(ci,arr,N,asdf,pending_delete); fprintf(stderr,"CHECKED\n");
}

// this is to make all qualifying tuples in range [a,b) are inside the main array
int merge_ripple(ci_type &ci, value_type *arr, int &N, multiset<int> &pins, multiset<int> &pdel, value_type a, value_type b){
  if (!pins.size() && !pdel.size()) return -1;
  value_type new_hi = -1;

  do_merge_ripple_insert(ci,arr,N,pins,pdel,a,b,new_hi);

  do_merge_ripple_delete(ci,arr,N,pins,pdel,a,b);

//  check(q);
  return new_hi;
}

int targeted_random_crack(ci_type &ci, value_type v, value_type *arr, int &N, int L, int R, int ncracks, int crack_at){
  while (ncracks-- > 0 && R - L > crack_at){    // split if the piece size is > CRACK_AT
    value_type X = arr[L + rand()%(R-L)];    // split in random position
    value_type X2 = arr[L + rand()%(R-L)];    // split in random position
    value_type X3 = arr[L + rand()%(R-L)];    // split in random position
    if (X==X2 && X==X3) break;          // guard agains too many duplicates

    // make X the median of the samples
    if (X2 > X) swap(X,X2);
    if (X > X3) swap(X,X3);
    if (X2 > X) swap(X,X2);

    int M = partition(arr, X, L,R);        // add crack on X
    add_crack(ci, N, X, M);
    if (v < X) R = M; else L = M;        // go to the correct sub-piece
  }
  return add_crack(ci, N, v, partition(arr, v, L,R));
}

int targeted_random_crack3(ci_type &ci, value_type v, value_type *arr, int &N, int L, int R, int crack_at){
  if (R - L > crack_at){    // split if the piece size is > CRACK_AT
    value_type X = arr[L + rand()%(R-L)];    // split in random position
    value_type X2 = arr[L + rand()%(R-L)];    // split in random position
    value_type X3 = arr[L + rand()%(R-L)];    // split in random position
    if (!(X==X2 && X==X3)){            // guard agains too many duplicates
      if (X2 > X) swap(X,X2);          // make X the median of the samples
      if (X > X3) swap(X,X3);
      if (X2 > X) swap(X,X2);

      if (v < X){
        int i1, i2;
        symmetric_crack3(arr, L, R, v, X, i1, i2);
        add_crack(ci, N, v, i1);
        add_crack(ci, N, X, i2);
        return i1;
      } else if (v > X){
        int i1, i2;
        symmetric_crack3(arr, L, R, X, v, i1, i2);
        add_crack(ci, N, X, i1);
        add_crack(ci, N, v, i2);
        return i2;
      }
    }
  }
  return add_crack(ci, N, v, partition(arr, v, L,R));
}



// MDD1R : Materialized DD1R
template <int CHECK>
int mdd1r_split_and_materialize(int L, int R, value_type a, value_type b){
  value_type X = arr[L + rand()%(R-L)];      // pick a pivot value X randomly in index [L,R)

  // optimization based on where [a,b] relative to X
  if (b < X){
    if (CHECK==0){                  // left check
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X){
          if (a<=arr[L]) marr[msize++] = (arr[L]);  // check left only
          L++;
        }
        while (L<=R && arr[R] >= X) R--;  // the arr[R] don't need to be checked for materialization!
        if (L<R) swap(arr[L],arr[R]);
      }
    } else if (CHECK==1){           // right check
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X){
          if (arr[L]<b) marr[msize++] = (arr[L]);  // check right only
          L++;
        }
        while (L<=R && arr[R] >= X) R--;  // the arr[R] don't need to be checked for materialization!
        if (L<R) swap(arr[L],arr[R]);
      }
    } else {                        // both checks
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X){
          if (a<=arr[L] && arr[L]<b) marr[msize++] = (arr[L]);  // check both and materialize
          L++;
        }
        while (L<=R && arr[R] >= X) R--;  // the arr[R] don't need to be checked for materialization!
        if (L<R) swap(arr[L],arr[R]);
      }
    }
  } else if (X <= a){
    if (CHECK==0){
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X) L++;    // the arr[L] don't need to be checked for materialization!
        while (L<=R && arr[R] >= X){
          if (a<=arr[R]) marr[msize++] = (arr[R]);  // check left
          R--;
        }
        if (L<R) swap(arr[L],arr[R]);
      }
    } else if (CHECK==1){
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X) L++;    // the arr[L] don't need to be checked for materialization!
        while (L<=R && arr[R] >= X){
          if (arr[R]<b) marr[msize++] = (arr[R]);  // check right
          R--;
        }
        if (L<R) swap(arr[L],arr[R]);
      }
    } else {
      for (R--; L<=R; ){            // split [L,R) based on X and materialize
        while (L<=R && arr[L] < X) L++;    // the arr[L] don't need to be checked for materialization!
        while (L<=R && arr[R] >= X){
          if (a<=arr[R] && arr[R]<b) marr[msize++] = (arr[R]);  // check both
          R--;
        }
        if (L<R) swap(arr[L],arr[R]);
      }
    }
  } else {
    for (R--; L<=R; ){              // split [L,R) based on X and materialize
      while (L<=R && arr[L] < X){
        if (a<=arr[L]) marr[msize++] = (arr[L]);  // materialize
        L++;
      }
      while (L<=R && arr[R] >= X){
        if (arr[R]<b) marr[msize++] = (arr[R]);  // materialize
        R--;
      }
      if (L<R) swap(arr[L],arr[R]);
    }
  }
  return add_crack(ci, N, X, L);          // add crack on random value X
}

int crack(value_type a, value_type b){
  merge_ripple(ci, arr, N, pins, pdel, a, b);  // merge qualified updates

  int L1, R1, i1; find_piece(ci, N, a, L1, R1);
  int L2, R2, i2; find_piece(ci, N, b, L2, R2);

  if (L1==L2){            // a and b is on the same piece
    assert(R1 == R2);
    split_ab(arr, L1,R1, a, b, i1, i2);  // 3-split in one scan 
  } else {                // a and b is on different piece
    i1 = partition(arr, a, L1, R1);
    i2 = partition(arr, b, L2, R2);
  }

  add_crack(ci, N, a, i1);
  add_crack(ci, N, b, i2);

  return i2 - i1;    // return number of qualified tuples
}

// do random crack on the target piece, then materialize
int random_crack_and_materialize(value_type a, value_type b){  
  merge_ripple(ci, arr, N, pins, pdel, a, b);  // merge qualified updates
  int L1,R1,i1; find_piece(ci, N, a, L1,R1);
  int L2,R2,i2; find_piece(ci, N, b, L2,R2);
  assert(L1 != L2 || R1 == R2);
  if (L1 < R1){
    value_type X1 = arr[L1 + rand()%(R1-L1)];
    add_crack(ci, N, X1, partition(arr,X1,L1,R1));
  }
  if (L1 != L2){
    if (L2 < R2){
      value_type X2 = arr[L2 + rand()%(R2-L2)];
      add_crack(ci, N, X2, partition(arr,X2,L2,R2));
    }
  }
  msize = 0;
  materialize_it(ci, arr, N, a, b, i1, i2, marr, msize);
  return msize + max(0, i2 - i1);
}
