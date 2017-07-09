#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <functional>
#include <set>
#include <stack>
#include <vector>

#include "random.h"
#include "time_it.h"

using namespace std;

#ifdef DBG
#define assert_dbg(x) assert(x)
#else
#define assert_dbg(x) ((void)0)
#endif

/*
// Count only using std::partition, without actually moving the data.
Bucket* hypothetical_partition(Bucket *B, int P) {
  for (; B; B = B->next) {
    int pos = partition(B->arr, B->arr + B->n, bind2nd(less<int>(), P)) -
B->arr;
    assert(pos >= 0 && pos <= B->n);
  }
  return NULL;
}

// Use if-else branch and move one-by-one.
Bucket* use_if(Bucket *B, int P) {
  Bucket *F = NULL; // Free-ed bucket.
  Bucket *L = new Bucket(), *Lh = L;
  Bucket *R = new Bucket();
  for (Bucket *next; B; B = next) {
    next = B->next;
    while (B->n) {
      int n = min(B->n, min(L->slack(), R->slack()));
      assert(n > 0);
      int *src = B->arr + B->n - 1;
      B->n -= n;
      while (n--){
        if (*src < P){
          L->arr[L->n++] = *(src--);
        } else {
          R->arr[R->n++] = *(src--);
        }
      }
      if (!B->n) { assert(!F || !next); F = B; }
      if (!L->slack()) { assert(F && !F->n); chain(L, F); F = NULL; }
      if (!R->slack()) { assert(F && !F->n); chain(R, F); F = NULL; }
    }
  }
  return Lh;
}

// Use memmove to move in bulk.
Bucket* use_partition(Bucket *B, int P) {
  Bucket *L = new Bucket(), *Lh = L;
  Bucket *R = new Bucket();
  set<Bucket*> F;
  // int nL = 0;
  for (Bucket *next; B; B = next) {
    next = B->next;
    int pos = partition(B->arr, B->arr + B->n, bind2nd(less<int>(), P)) -
B->arr;
    if (L->slack() >= pos) {
      memmove(L->arr + L->n, B->arr, sizeof(int) * pos);
      L->n += pos;
      B->n = 0;
      F.insert(B);
    } else {
      int m = L->slack();
      memmove(L->arr + L->n, B->arr + pos - m, sizeof(int) * m);
      L->n += m;
      B->n = pos - m;
      chain(L, B);
      // nL += BSIZE;
    }

    if (R->slack() >= BSIZE - pos){
      memmove(R->arr + R->n, B->arr + pos, sizeof(int)*(BSIZE-pos));
      R->n += BSIZE-pos;
      if (L != B){
        F.insert(B);
        B->n = 0;
      }
    } else {
      int m = R->slack();
      memmove(R->arr + R->n, B->arr + pos, sizeof(int)* m);
      R->n += m;
      pos += m;
      assert(!F.empty());
      R = *F.begin();
      F.erase(F.begin());
      assert(R->n == 0);
      memmove(R->arr + R->n, B->arr + pos, sizeof(int)*(BSIZE-pos));
      R->n = BSIZE-pos;
    }
  }
  // nL += L->n;
  // fprintf(stderr, "Ln = %d\n", nL);
  return Lh;
}

// std::partition minimizing data movement.
Bucket* use_partition_minimize_move(Bucket *B, int P) {
  vector<pair<int, Bucket*> > L, R;
  for (Bucket *next; B; B = next) {
    next = B->next;
    int pos = partition(B->arr, B->arr + B->n, bind2nd(less<int>(), P)) -
B->arr;
    if (pos * 2 > B->n) L.push_back(make_pair(pos, B));
    else R.push_back(make_pair(pos, B));
  }

  int nL = 0;
  Bucket *ret = NULL;
  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());
  while (!L.empty() && !R.empty()) {
    pair<int, Bucket*> &pL = L.back(); L.pop_back();
    pair<int, Bucket*> &pR = R.back(); R.pop_back();
    int *a = pL.second->arr, &i = pL.first;
    int *b = pR.second->arr, &j = pR.first;
    int m = min(pL.second->n - i, j);
    while (m--) swap(a[i++], b[j--]);
    if (i < pL.second->n) L.push_back(pL);
    else { nL += BSIZE; assert(pL.second->n == BSIZE); chain(ret, pL.second); }
    if (j > 0) R.push_back(pR);
  }

  for (int k=0; k < (int)L.size(); ) {
    if (k+1 == (int)L.size()){
      nL += L.back().first;
      L.back().second->n = L.back().first;
      chain(ret, L.back().second);
      break;
    }
    pair<int, Bucket*> &pL = L[k], &pR = L.back();
    int *a = pL.second->arr, &i = pL.first;
    int *b = pR.second->arr, &j = pR.first;
    int m = min(pL.second->n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i >= pL.second->n){ k++; nL += BSIZE; assert(pL.second->n == BSIZE);
chain(ret, pL.second); }
    if (j == 0) L.pop_back();
  }

  for (int k=0; k < (int)R.size(); ){
    if (k+1 == (int)R.size()) {
      nL += R.back().first;
      R.back().second->n = R.back().first;
      chain(ret, R.back().second);
      break;
    }
    pair<int, Bucket*> &pL = R[k], &pR = R.back();
    int *a = pL.second->arr, &i = pL.first;
    int *b = pR.second->arr, &j = pR.first;
    int m = min(i, pR.second->n - j);
    while (m--) swap(a[--i], b[j++]);
    if (j >= pR.second->n){ R.pop_back(); nL += BSIZE; assert(pR.second->n ==
BSIZE); chain(ret, pR.second); }
    if (i == 0){ k++;  }
  }

  fprintf(stderr, "nL = %d\n", nL);
  return ret;
}

int use_nb(Bucket *B, int P) {
  int F = -1, L = NBUCKETS, R = L+1, nL = L->n = B[R].n = 0;
  REP(i,NBUCKETS){
    while (B[i].n){
      if (L==-1){ assert(F!=-1); L = F; F = -1; }
      else if (R==-1){ assert(F!=-1); R = F; F = -1; }
      int n = min(B[i].n, min(B[L].slack(), B[R].slack()));
      assert(n > 0);
      int *src = B[i].arr + B[i].n - 1;
      B[i].n -= n;
      int cnt[2] = { B[L].n, B[R].n };
      int *arr[2] = { B[L].arr, B[R].arr };
      while (n--){
        int t = *src >= P;
        arr[t][cnt[t]++] = *(src--);
      }
      B[L].n = cnt[0];
      B[R].n = cnt[1];
      if (!B[i].n){ assert(F==-1); F = i; }
      if (!B[L].slack()) L = -1, nL += BSIZE;
      if (!B[R].slack()) R = -1;
    }
  }
  if (L!=-1) nL += B[L].n;
  return nL;
  // fprintf(stderr, "%d / %d = %.3lf\n", nL,N,100.0*nL/N);
}

int use_nbo(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS){
    int cnt[2] = { 0, B[R].n }, j = 0;
    int *arr[2] = { B[i].arr, B[R].arr };
    while (j < B[i].n){
      // fprintf(stderr,"i = %d, j = %d/%d, R = %d, %d %d,
%d\n",i,j,B[i].n,R,cnt[0],B[i].n,B[R].slack());
      int n = min(B[i].n - j, B[R].slack()); assert(n > 0);
      int *src = B[i].arr + j; j += n;
      while (n--){
        int t = *src >= P;
        arr[t][cnt[t]++] = *(src++);
      }
      B[R].n = cnt[1];
      if (!B[R].slack()){
        R = F++, B[R].n = 0;
        cnt[1] = 0;
        arr[1] = B[R].arr;
      }
      assert(F < 2*NBUCKETS+10);
    }
    assert(j == B[i].n);
    B[i].n = cnt[0];
    nL += cnt[0];
  }
  return nL;
}

int use_cpy(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS) REP(j, B[i].n){
    nL += B[i].arr[j] < P;
    B[F+i].arr[j] = B[i].arr[j];
  }
  return nL;
}

int use_cpyp2(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS){
    int h = B[i].n / 2;
    REP(j, h){
      nL += (B[i].arr[j] < P) + (B[i].arr[h+j] < P);
      B[F+i].arr[j] = B[i].arr[j];
      B[F+i].arr[j+h] = B[i].arr[j+h];
    }
  }
  return nL;
}

int use_cpyp4(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS){
    for (int j=0; j<B[i].n; j+=4){
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) +
(B[i].arr[j+3] < P);
      B[F+i].arr[j] = B[i].arr[j];
      B[F+i].arr[j+1] = B[i].arr[j+1];
      B[F+i].arr[j+2] = B[i].arr[j+2];
      B[F+i].arr[j+3] = B[i].arr[j+3];
    }
  }
  return nL;
}

int use_cpyp8(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS){
    for (int j=0; j<B[i].n; j+=8){
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) +
(B[i].arr[j+3] < P) +
          (B[i].arr[j+4] < P) + (B[i].arr[j+5] < P) + (B[i].arr[j+6] < P) +
(B[i].arr[j+7] < P);
      B[F+i].arr[j] = B[i].arr[j];
      B[F+i].arr[j+1] = B[i].arr[j+1];
      B[F+i].arr[j+2] = B[i].arr[j+2];
      B[F+i].arr[j+3] = B[i].arr[j+3];
      B[F+i].arr[j+4] = B[i].arr[j+4];
      B[F+i].arr[j+5] = B[i].arr[j+5];
      B[F+i].arr[j+6] = B[i].arr[j+6];
      B[F+i].arr[j+7] = B[i].arr[j+7];
    }
  }
  return nL;
}

int use_cpyp16(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  REP(i,NBUCKETS){
    for (int j=0; j<B[i].n; j+=16){
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) +
(B[i].arr[j+3] < P) +
          (B[i].arr[j+4] < P) + (B[i].arr[j+5] < P) + (B[i].arr[j+6] < P) +
(B[i].arr[j+7] < P) +
          (B[i].arr[j+8] < P) + (B[i].arr[j+9] < P) + (B[i].arr[j+10] < P) +
(B[i].arr[j+11] < P) +
          (B[i].arr[j+12] < P) + (B[i].arr[j+13] < P) + (B[i].arr[j+14] < P) +
(B[i].arr[j+15] < P);
      B[F+i].arr[j] = B[i].arr[j];
      B[F+i].arr[j+1] = B[i].arr[j+1];
      B[F+i].arr[j+2] = B[i].arr[j+2];
      B[F+i].arr[j+3] = B[i].arr[j+3];
      B[F+i].arr[j+4] = B[i].arr[j+4];
      B[F+i].arr[j+5] = B[i].arr[j+5];
      B[F+i].arr[j+6] = B[i].arr[j+6];
      B[F+i].arr[j+7] = B[i].arr[j+7];
      B[F+i].arr[j+8] = B[i].arr[j+8];
      B[F+i].arr[j+9] = B[i].arr[j+9];
      B[F+i].arr[j+10] = B[i].arr[j+10];
      B[F+i].arr[j+11] = B[i].arr[j+11];
      B[F+i].arr[j+12] = B[i].arr[j+12];
      B[F+i].arr[j+13] = B[i].arr[j+13];
      B[F+i].arr[j+14] = B[i].arr[j+14];
      B[F+i].arr[j+15] = B[i].arr[j+15];
    }
  }
  return nL;
}

int use_memo(Bucket *B, int P) {
  int F = NBUCKETS, R = F++, nL = B[R].n = 0;
  char amt[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
  REP(i,NBUCKETS){
    for (int j=0; j<B[i].n; j+=16){
      char c1 = ((B[i].arr[j] < P) << 3) | ((B[i].arr[j+1] < P) << 2) |
((B[i].arr[j+2] < P) << 1) | ((B[i].arr[j+3] < P) << 0);
      char c2 = ((B[i].arr[j+4] < P) << 3) | ((B[i].arr[j+5] < P) << 2) |
((B[i].arr[j+6] < P) << 1) | ((B[i].arr[j+7] < P) << 0);
      char c3 = ((B[i].arr[j+8] < P) << 3) | ((B[i].arr[j+9] < P) << 2) |
((B[i].arr[j+10] < P) << 1) | ((B[i].arr[j+11] < P) << 0);
      char c4 = ((B[i].arr[j+12] < P) << 3) | ((B[i].arr[j+13] < P) << 2) |
((B[i].arr[j+14] < P) << 1) | ((B[i].arr[j+15] < P) << 0);
      nL += amt[c1] + amt[c2] + amt[c3] + amt[c4];
    }
  }
  return nL;
}

int use_cnt4(Bucket *B, int P) {
  int n[4] = { 0 };
  REP(i,NBUCKETS){
    for (int j=0; j<B[i].n; j+=4){
      n[0] += (B[i].arr[j] < P);
      n[1] += (B[i].arr[j+1] < P);
      n[2] += (B[i].arr[j+2] < P);
      n[3] += (B[i].arr[j+3] < P);
    }
  }
  return n[0]+n[1]+n[2]+n[3];
}

int use_cpy4(Bucket *B, int P) {
  int nL = 0, *a[4] = { B[NBUCKETS].arr, B[NBUCKETS+1].arr, B[NBUCKETS+2].arr,
B[NBUCKETS+3].arr};
  REP(i,NBUCKETS){
    int n[4] = { 0 };
    for (int j=0; j<B[i].n; j+=4){
      a[0][n[0]++] = B[i].arr[j];
      a[1][n[1]++] = B[i].arr[j+1];
      a[2][n[2]++] = B[i].arr[j+2];
      a[3][n[3]++] = B[i].arr[j+3];
    }
    nL += n[0] + n[1] + n[2] + n[3];
  }
  return nL;
}

int use_noswap(Bucket *B, int P) {
  int nL = 0;
  int hi[BSIZE], lo[BSIZE], nhi=0, nlo=0;
  REP(i,NBUCKETS){
    nlo = nhi = 0;
    for (int j=0; j<B[i].n; j+=2){
      lo[nlo] = B[i].arr[j];
      hi[nhi] = B[i].arr[j+1];
      nlo += B[i].arr[j] < P;
      nhi += B[i].arr[j+1] < P;
    }
    nL += nlo + nhi;
  }
  return nL;
}

int use_nbswap(Bucket *B, int P) {
  int nL = 0, i = 0, L=-1, R=-1;
  int hi[BSIZE], lo[BSIZE], nhi=0, nlo=0;
  int t1=0,t2=0;
  while (true){
    if (nhi && nlo){
      int m = min(nhi, nlo);
      assert(m > 0);
      assert(L != -1);
      assert(R != -1);
      while (m--){
        swap(B[L].arr[hi[--nhi]], B[R].arr[lo[--nlo]]);
      }
      // nhi -= m;
      // nlo -= m;

      if (nhi==0) L = -1;
      if (nlo==0) R = -1;
    } else if (L == -1){
      if (i == NBUCKETS) break;
      nL += BSIZE;
      L = i++;
    } else if (R == -1){
      if (i == NBUCKETS) break;
      R = i++;
    } else if (nhi == 0){
      assert(L!=-1);
      // timit();
      for (int j=0; j<B[L].n; j++){
        hi[nhi] = j;
        nhi += B[L].arr[j] >= P;
      }
      // t1 += timit();
      if (nhi == 0) L = -1;
    } else if (nlo == 0){
      assert(R!=-1);
      // timit();
      for (int j=0; j<B[R].n; j++){
        lo[nlo] = j;
        nlo += B[R].arr[j] < P;
      }
      // t2 += timit();
      if (nlo == 0) R = -1;
    } else {
      assert(0);
    }
  }
  if (L!=-1){ nL -= BSIZE; REP(j,B[L].n) nL += B[L].arr[j] < P; }
  if (R!=-1) REP(j,B[R].n) nL += B[R].arr[j] < P;
  return nL;
}

template <int T>
int parsplit(int *arr, int n, int *a1, int *a2, int &i1, int &i2, int P){
  for (int k=0; k<n; k+=2){
    a1[i1] = k;
    a2[i2] = k+1;
    if (T==0){
      i1 += arr[k] >= P;
      i2 += arr[k+1] >= P;
    } else {
      i1 += arr[k] < P;
      i2 += arr[k+1] < P;
    }
  }
}

template <int T>
int parsplit(int *arr, int n, int *a, int &i, int P){
  for (int k=0; k<n; k++){
    a[i] = k;
    if (T==0){
      i += arr[k] >= P;
    } else {
      i += arr[k] < P;
    }
  }
}

int fusion(int Rb, int *lo, int &nlo, int Lb, int *hi, int &nhi){
  int n = min(nhi, nlo); assert(n > 0);
  while (n--) swap(B[Rb].arr[lo[--nlo]], B[Lb].arr[hi[--nhi]]);
}

int use_nbsortparnb(Bucket *B, int P) {
  int nL = 0;
  t1 = t2 = t3 = t4 = t5 = 0.0;
  vector<pair<int,int> > L, R, aL, aR;

  REP(i,NBUCKETS){
    int nlow = 0;
    REP(j,9) nlow += B[i].arr[rand()%B[i].n] < P;
    if (nlow < 1 || 9 <= nlow){
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P))
- B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else
R.push_back(make_pair(pos, i));
    } else if (nlow >= 5){
      aL.push_back(make_pair(nlow,i));
    } else {
      aR.push_back(make_pair(nlow,i));
    }
  }

  // t1 += timit();

  // fprintf(stderr,"a");
  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());

  while (!L.empty() && !R.empty()){
    pair<int,int> &pL = L.back(); L.pop_back();
    pair<int,int> &pR = R.back(); R.pop_back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[j--]);
    if (i<B[pL.second].n) L.push_back(pL); else nL += BSIZE;
    if (j>0) R.push_back(pR);
  }

  // fprintf(stderr,"b");
  for (int k=0; k < L.size(); ){
    if (k+1 == L.size()){ aL.push_back(L.back()); break; }
    pair<int,int> &pL = L[k], &pR = L.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i>=B[pL.second].n){ k++; nL += BSIZE; }
    if (j==0) L.pop_back();
  }

  // fprintf(stderr,"c");
  for (int k=0; k < R.size(); ){
    if (k+1 == R.size()){ aR.push_back(R.back()); break; }
    pair<int,int> &pL = R[k], &pR = R.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(i, B[pR.second].n - j);
    while (m--) swap(a[--i], b[j++]);
    if (j>=B[pR.second].n){ R.pop_back(); nL += BSIZE; }
    if (i==0) k++;
  }

  // t2 += timit();

  // fprintf(stderr,"d");

  L.swap(aL); sort(L.rbegin(), L.rend());
  R.swap(aR); sort(R.begin(), R.end());
  int hi[2][BSIZE/2], lo[2][BSIZE/2], nhi[2]={0}, nlo[2]={0}, Lb = -1, Rb = -1;
  while ((!L.empty() || nhi[0] || nhi[1]) && (!R.empty() || nlo[0] || nlo[1])){
    if (nhi[0] == 0 && nhi[1] == 0){
      Lb = L.back().second; L.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    }
    if (nlo[0] == 0 && nlo[1] == 0){
      Rb = R.back().second; R.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    }
    REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
      fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
  }
  assert(L.empty() || R.empty());
  // fprintf(stderr,"e");

  // t3 += timit();

  for (int k=0; ; ){
    if (nlo[0] == 0 && nlo[1] == 0){
      if (k == (int)R.size()) break;
      Rb = R[k++].second;
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    } else if (nhi[0] == 0 && nhi[1] == 0){
      if (k == (int)R.size()) break;
      Lb = R.back().second; R.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    } else {
      REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
        fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
    }
  }
  // fprintf(stderr,"f");

  // t4 += timit();

  for (int k=0; ; ){
    if (nhi[0] == 0 && nhi[1] == 0){
      if (k == (int)L.size()) break;
      Lb = L[k++].second;
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    } else if (nlo[0] == 0 && nlo[1] == 0){
      if (k == (int)L.size()) break;
      Rb = L.back().second; L.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    } else {
      REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
        fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
    }
  }
  // fprintf(stderr,"g");
  if (nhi[0]+nhi[1]){ nL -= BSIZE; REP(k,B[Lb].n) nL += B[Lb].arr[k] < P; }
  if (nlo[0]+nlo[1]){ REP(k,B[Rb].n) nL += B[Rb].arr[k] < P; }

  // t5 += timit();

  return nL;
}

int use_nbsortswap2(Bucket *B, int P) {
  int nL = 0;
  t1 = t2 = t3 = t4 = t5 = 0.0;
  vector<pair<int,int> > L, R, aL, aR;

  REP(i,NBUCKETS){
    int nlow = 0;
    REP(j,9) nlow += B[i].arr[rand()%B[i].n] < P;
    if (nlow < 1 || 9 <= nlow){
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P))
- B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else
R.push_back(make_pair(pos, i));
    } else if (nlow >= 5){
      aL.push_back(make_pair(nlow,i));
    } else {
      aR.push_back(make_pair(nlow,i));
    }
  }

  // t1 += timit();

  // fprintf(stderr,"a");
  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());

  while (!L.empty() && !R.empty()){
    pair<int,int> &pL = L.back(); L.pop_back();
    pair<int,int> &pR = R.back(); R.pop_back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[j--]);
    if (i<B[pL.second].n) L.push_back(pL); else nL += BSIZE;
    if (j>0) R.push_back(pR);
  }

  // fprintf(stderr,"b");
  for (int k=0; k < L.size(); ){
    if (k+1 == L.size()){ aL.push_back(L.back()); break; }
    pair<int,int> &pL = L[k], &pR = L.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i>=B[pL.second].n){ k++; nL += BSIZE; }
    if (j==0) L.pop_back();
  }

  // fprintf(stderr,"c");
  for (int k=0; k < R.size(); ){
    if (k+1 == R.size()){ aR.push_back(R.back()); break; }
    pair<int,int> &pL = R[k], &pR = R.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(i, B[pR.second].n - j);
    while (m--) swap(a[--i], b[j++]);
    if (j>=B[pR.second].n){ R.pop_back(); nL += BSIZE; }
    if (i==0) k++;
  }

  // t2 += timit();

  // fprintf(stderr,"d");

  L.swap(aL); sort(L.rbegin(), L.rend());
  R.swap(aR); sort(R.begin(), R.end());
  int hi[2][BSIZE/2], lo[2][BSIZE/2], nhi[2]={0}, nlo[2]={0}, Lb = -1, Rb = -1;
  while ((!L.empty() || nhi[0] || nhi[1]) && (!R.empty() || nlo[0] || nlo[1])){
    if (nhi[0] == 0 && nhi[1] == 0){
      Lb = L.back().second; L.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    }
    if (nlo[0] == 0 && nlo[1] == 0){
      Rb = R.back().second; R.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    }
    REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
      fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
  }
  assert(L.empty() || R.empty());
  // fprintf(stderr,"e");

  // t3 += timit();

  for (int k=0; ; ){
    if (nlo[0] == 0 && nlo[1] == 0){
      if (k == (int)R.size()) break;
      Rb = R[k++].second;
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    } else if (nhi[0] == 0 && nhi[1] == 0){
      if (k == (int)R.size()) break;
      Lb = R.back().second; R.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    } else {
      REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
        fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
    }
  }
  // fprintf(stderr,"f");

  // t4 += timit();

  for (int k=0; ; ){
    if (nhi[0] == 0 && nhi[1] == 0){
      if (k == (int)L.size()) break;
      Lb = L[k++].second;
      parsplit<0>(B[Lb].arr, B[Lb].n, hi[0], hi[1], nhi[0], nhi[1], P);
      nL += BSIZE;
    } else if (nlo[0] == 0 && nlo[1] == 0){
      if (k == (int)L.size()) break;
      Rb = L.back().second; L.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo[0], lo[1], nlo[0], nlo[1], P);
    } else {
      REP(i,2) REP(j,2) if (nhi[j] && nlo[i])
        fusion(Rb,lo[i],nlo[i], Lb,hi[j],nhi[j]);
    }
  }
  // fprintf(stderr,"g");
  if (nhi[0]+nhi[1]){ nL -= BSIZE; REP(k,B[Lb].n) nL += B[Lb].arr[k] < P; }
  if (nlo[0]+nlo[1]){ REP(k,B[Rb].n) nL += B[Rb].arr[k] < P; }

  // t5 += timit();

  return nL;
}

int use_nbsortswap(Bucket *B, int P) {
  int nL = 0;
  t1 = t2 = t3 = t4 = t5 = 0.0;
  vector<pair<int,int> > L, R, aL, aR;

  REP(i,NBUCKETS){
    int nlow = 0;
    REP(j,9) nlow += B[i].arr[rand()%B[i].n] < P;
    if (nlow < 1 || 9 <= nlow){
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P))
- B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else
R.push_back(make_pair(pos, i));
    } else if (nlow >= 5){
      aL.push_back(make_pair(nlow,i));
    } else {
      aR.push_back(make_pair(nlow,i));
    }
  }

  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());

  // t1 += timit();

  while (!L.empty() && !R.empty()){
    pair<int,int> &pL = L.back(); L.pop_back();
    pair<int,int> &pR = R.back(); R.pop_back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[j--]);
    if (i<B[pL.second].n) L.push_back(pL); else nL += BSIZE;
    if (j>0) R.push_back(pR);
  }

  for (int k=0; k < (int)L.size(); ){
    if (k+1 == (int)L.size()){ aL.push_back(L.back()); break; }
    pair<int,int> &pL = L[k], &pR = L.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i>=B[pL.second].n){ k++; nL += BSIZE; }
    if (j==0) L.pop_back();
  }

  for (int k=0; k < (int)R.size(); ){
    if (k+1 == (int)R.size()){ aR.push_back(R.back()); break; }
    pair<int,int> &pL = R[k], &pR = R.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(i, B[pR.second].n - j);
    while (m--) swap(a[--i], b[j++]);
    if (j>=B[pR.second].n){ R.pop_back(); nL += BSIZE; }
    if (i==0){ k++; }
  }

  // t2 += timit();

  L = aL; R = aR;
  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());
  int Lb = -1, Rb = -1;
  int hi[BSIZE], lo[BSIZE], nhi=0, nlo=0;
  int *lo2 = lo+BSIZE/2, nlo2 = 0;
  int *hi2 = hi+BSIZE/2, nhi2 = 0;
  while ((!L.empty() || nhi) && (!R.empty() || nlo)){
    if (nhi == 0){
      Lb = L.back().second; L.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi, nhi, P);
      nL += BSIZE;
    }
    if (nlo == 0){
      Rb = R.back().second; R.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo, nlo, P);
    }
    if (nhi && nlo){
      int n = min(nhi, nlo); assert(n > 0);
      while (n--) swap(B[Rb].arr[lo[--nlo]], B[Lb].arr[hi[--nhi]]);
      // nhi -= n;
      // nlo -= n;
    }
  }
  // if (nhi){
  //  nL -= BSIZE;
  //  L.push_back(make_pair(-1, Lb));
  // }
  // if (nlo){
  //  R.push_back(make_pair(-1, Rb));
  // }

  assert(L.empty() || R.empty());
  // t3 += timit();

  // nlo = nhi = 0;
  for (int k=0; ; ){
    if (nlo == 0){
      if (k == (int)R.size()) break;
      Rb = R[k++].second;
      parsplit<1>(B[Rb].arr, B[Rb].n, lo, nlo, P);
    } else if (nhi == 0){
      if (k == (int)R.size()) break;
      Lb = R.back().second; R.pop_back();
      parsplit<0>(B[Lb].arr, B[Lb].n, hi, nhi, P);
      nL += BSIZE;
    } else {
      int n = min(nlo, nhi); assert(n > 0);
      while (n--) swap(B[Rb].arr[lo[--nlo]], B[Lb].arr[hi[--nhi]]);
    }
  }
  // if (nhi){ nL -= BSIZE; REP(k,B[Lb].n) nL += B[Lb].arr[k] < P; nhi = 0; }
  // if (nlo){ REP(k,B[Rb].n) nL += B[Rb].arr[k] < P; nlo = 0; }

  // t4 += timit();

  // nlo = nhi = 0;
  for (int k=0; ; ){
    if (nhi == 0){
      if (k == (int)L.size()) break;
      Lb = L[k++].second;
      parsplit<0>(B[Lb].arr, B[Lb].n, hi, nhi, P);
      // fprintf(stderr, "nhi = %d\n", nhi);
      nL += BSIZE;
    } else if (nlo == 0){
      if (k == (int)L.size()) break;
      Rb = L.back().second; L.pop_back();
      parsplit<1>(B[Rb].arr, B[Rb].n, lo, nlo, P);
      // fprintf(stderr, "nlo = %d\n", nlo);
    } else {
      int n = min(nhi, nlo); assert(n > 0);
      while (n--) swap(B[Rb].arr[lo[--nlo]], B[Lb].arr[hi[--nhi]]);
    }
  }
  if (nhi){ nL -= BSIZE; REP(k,B[Lb].n) nL += B[Lb].arr[k] < P; nhi = 0; }
  if (nlo){ REP(k,B[Rb].n) nL += B[Rb].arr[k] < P; nlo = 0; }

  // t5 += timit();

  return nL;
}

int use_parnb(Bucket *B, int P) {
  int nL = 0;
  t1 = t2 = t3 = t4 = t5 = 0.0;
  vector<pair<int,int> > L, R, nb;
  L.reserve(NBUCKETS);
  R.reserve(NBUCKETS);
  nb.reserve(NBUCKETS);

  REP(i,NBUCKETS){
    int nlow = 0;
    REP(j,9) nlow += B[i].arr[rand()%B[i].n] < P;
    if (nlow < 1 || 9 <= nlow){
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P))
- B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else
R.push_back(make_pair(pos, i));
    } else {
      nb.push_back(make_pair(nlow,i));
    }
  }

  // t1 += timit();

  sort(L.rbegin(), L.rend());
  sort(R.begin(), R.end());

  while (!L.empty() && !R.empty()){
    pair<int,int> &pL = L.back(); L.pop_back();
    pair<int,int> &pR = R.back(); R.pop_back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[j--]);
    if (i<B[pL.second].n) L.push_back(pL); else nL += BSIZE;
    if (j>0) R.push_back(pR);
  }

  // fprintf(stderr,"b");
  for (int k=0; k < L.size(); ){
    if (k+1 == L.size()){ nL += L.back().first; break; }
    pair<int,int> &pL = L[k], &pR = L.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i>=B[pL.second].n){ k++; nL += BSIZE; }
    if (j==0) L.pop_back();
  }

  // fprintf(stderr,"c");
  for (int k=0; k < R.size(); ){
    if (k+1 == R.size()){ nL += R.back().first; break; }
    pair<int,int> &pL = R[k], &pR = R.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(i, B[pR.second].n - j);
    while (m--) swap(a[--i], b[j++]);
    if (j>=B[pR.second].n){ R.pop_back(); nL += BSIZE; }
    if (i==0) k++;
  }

  // t2 += timit();

  int F1 = -1, L1 = NBUCKETS, R1 = L1+1;
  int F2 = -1, L2 = R1+1, R2 = L2+1, k = 0;
  B[L1].n = B[R1].n = B[L2].n = B[R2].n = 0;
  for ( ; k+1 < (int)nb.size(); k+=2){
    int i = nb[k].second, j = nb[k+1].second;
    while (B[i].n && B[j].n){
      if (L1==-1){ assert(F1!=-1); L1 = F1; F1 = -1; } else if (R1==-1){
assert(F1!=-1); R1 = F1; F1 = -1; }
      if (L2==-1){ assert(F2!=-1); L2 = F2; F2 = -1; } else if (R2==-1){
assert(F2!=-1); R2 = F2; F2 = -1; }
      int n1 = min(B[i].n, min(B[L1].slack(), B[R1].slack())); assert(n1 > 0);
      int n2 = min(B[j].n, min(B[L2].slack(), B[R2].slack())); assert(n2 > 0);
      int n = min(n1,n2); assert(n > 0);
      int *src1 = B[i].arr + B[i].n - 1;
      int *src2 = B[j].arr + B[j].n - 1;
      B[i].n -= n;
      B[j].n -= n;
      int cnt1[2] = { B[L1].n, B[R1].n }, *arr1[2] = { B[L1].arr, B[R1].arr };
      int cnt2[2] = { B[L2].n, B[R2].n }, *arr2[2] = { B[L2].arr, B[R2].arr };
      while (n--){
        char t1 = *src1 >= P;
        char t2 = *src2 >= P;
        arr1[t1][cnt1[t1]++] = *(src1--);
        arr2[t2][cnt2[t2]++] = *(src2--);
      }
      B[L1].n = cnt1[0]; B[R1].n = cnt1[1];
      B[L2].n = cnt2[0]; B[R2].n = cnt2[1];

      if (!B[i].n){ assert(F1==-1); F1 = i; }
      if (!B[j].n){ assert(F2==-1); F2 = j; }
      if (!B[L1].slack()) L1 = -1, nL += BSIZE;
      if (!B[L2].slack()) L2 = -1, nL += BSIZE;
      if (!B[R1].slack()) R1 = -1;
      if (!B[R2].slack()) R2 = -1;
    }
  }
  if (L1!=-1) nL += B[L1].n;
  if (L2!=-1) nL += B[L2].n;

  if (k < (int) nb.size()){
    Bucket &b = B[nb[k].second];
    nL += partition(b.arr, b.arr + b.n, bind2nd(less<int>(), P)) - b.arr;
  }

  // t3 += timit();

  return nL;
}


int use_partest(Bucket *B, int P) {
  int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
  int nL = 0, lo[BSIZE]={0};

  REP(i,NBUCKETS){
//    parsplit<1>(B[i].arr, B[i].n, lo, lo2, nlo, nlo2, P);
    register int j = 0;
    for (int k=0; k<B[i].n; k++){
      // lo[nlo2] = k+1;
      assert(j < BSIZE);
      lo[j] += k;
      j += (B[i].arr[k] < P);
      // nlo2 += (B[i].arr[k+1] < P) << 1;
    }
    // B[i].n = j;
  }

  return nL;
}

int swapit(int &nhi, int &nlo, int &L, int &R, int *hi, int *lo){
  int m = min(nhi, nlo); assert(m > 0);
  assert(L != -1);
  assert(R != -1);
  // while (m--) swap(B[L].arr[hi[--nhi]], B[R].arr[lo[--nlo]]);
  nhi -= m;
  nlo -= m;
}

int use_nbswap2(Bucket *B, int P) {
  int nL = 0, i = 0, L=-1, R=-1;
  int hia[BSIZE], loa[BSIZE];
  int *hi[2] = {hia,hia+1};
  int *lo[2] = {loa,loa+1};
  int nhi[2]={0}, nlo[2]={0};
  int t1=0,t2=0;
  while (true){
    if (nhi[0] && nlo[0]){
      swapit(nhi[0], nlo[0], L, R, hi[0], lo[0]);
      if (nhi[0]==0 && nhi[1]==0) L = -1;
      if (nlo[0]==0 && nlo[1]==0) R = -1;
    } else if (nhi[0] && nlo[1]){
      swapit(nhi[0], nlo[1], L, R, hi[0], lo[1]);
      if (nhi[0]==0 && nhi[1]==0) L = -1;
      if (nlo[0]==0 && nlo[1]==0) R = -1;
    } else if (nhi[1] && nlo[0]){
      swapit(nhi[1], nlo[0], L, R, hi[1], lo[0]);
      if (nhi[0]==0 && nhi[1]==0) L = -1;
      if (nlo[0]==0 && nlo[1]==0) R = -1;
    } else if (nhi[1] && nlo[1]){
      swapit(nhi[1], nlo[1], L, R, hi[1], lo[1]);
      if (nhi[0]==0 && nhi[1]==0) L = -1;
      if (nlo[0]==0 && nlo[1]==0) R = -1;
    } else if (L == -1){
      if (i == NBUCKETS) break;
      nL += BSIZE;
      L = i++;
    } else if (R == -1){
      if (i == NBUCKETS) break;
      R = i++;
    } else if (nhi[0] == 0 && nhi[1] == 0){
      assert(L!=-1);
      for (int j=0; j<B[L].n; j+=2){
        hi[0][nhi[0]] = j;
        hi[1][nhi[1]] = j+1;
        nhi[0] += B[L].arr[j] >= P;
        nhi[1] += B[L].arr[j+1] >= P;
      }
      if (nhi[0] == 0 && nhi[1] == 0) L = -1;
    } else if (nlo[0] == 0 && nlo[1] == 0){
      assert(R!=-1);
      for (int j=0; j<B[R].n; j+=2){
        lo[0][nlo[0]] = j;
        lo[1][nlo[1]] = j+1;
        nlo[0] += (B[R].arr[j] < P) << 1;
        nlo[1] += (B[R].arr[j+1] < P) << 1;
      }
      if (nlo[0] == 0 && nlo[1] == 0) R = -1;
    } else {
      assert(0);
    }
  }
  if (L!=-1){ nL -= BSIZE; REP(j,B[L].n) nL += B[L].arr[j] < P; }
  if (R!=-1) REP(j,B[R].n) nL += B[R].arr[j] < P;
  return nL;
}

int use_nbp(Bucket *B, int P) {
  int F1 = -1, L1 = NBUCKETS, R1 = L1+1;
  int F2 = -1, L2 = R1+1, R2 = L2+1, nL = 0;
  B[L1].n = B[R1].n = B[L2].n = B[R2].n = 0;
  for (int i=0; i<NBUCKETS; i+=2){
    while (B[i].n && B[i+1].n){
      if (L1==-1){ assert(F1!=-1); L1 = F1; F1 = -1; } else if (R1==-1){
assert(F1!=-1); R1 = F1; F1 = -1; }
      if (L2==-1){ assert(F2!=-1); L2 = F2; F2 = -1; } else if (R2==-1){
assert(F2!=-1); R2 = F2; F2 = -1; }
      int n1 = min(B[i].n, min(B[L1].slack(), B[R1].slack())); assert(n1 > 0);
      int n2 = min(B[i+1].n, min(B[L2].slack(), B[R2].slack())); assert(n2 > 0);
      int n = min(n1,n2); assert(n > 0);
      int *src1 = B[i].arr + B[i].n - 1;
      int *src2 = B[i+1].arr + B[i+1].n - 1;
      B[i].n -= n;
      B[i+1].n -= n;
      int cnt1[2] = { B[L1].n, B[R1].n }, *arr1[2] = { B[L1].arr, B[R1].arr };
      int cnt2[2] = { B[L2].n, B[R2].n }, *arr2[2] = { B[L2].arr, B[R2].arr };
      while (n--){
        char t1 = *src1 >= P;
        char t2 = *src2 >= P;
        arr1[t1][cnt1[t1]++] = *(src1--);
        arr2[t2][cnt2[t2]++] = *(src2--);
      }
      B[L1].n = cnt1[0]; B[R1].n = cnt1[1];
      B[L2].n = cnt2[0]; B[R2].n = cnt2[1];

      if (!B[i].n){ assert(F1==-1); F1 = i; }
      if (!B[i+1].n){ assert(F2==-1); F2 = i+1; }
      if (!B[L1].slack()) L1 = -1, nL += BSIZE;
      if (!B[L2].slack()) L2 = -1, nL += BSIZE;
      if (!B[R1].slack()) R1 = -1;
      if (!B[R2].slack()) R2 = -1;
    }
  }
  if (L1!=-1) nL += B[L1].n;
  if (L2!=-1) nL += B[L2].n;
  return nL;
}
*/

#define BSIZE (1 << 12)

static long long* allocate_array(int sz);

class Bucket {
 public:
  long long* arr;
  int n;

  Bucket(long long* a, int sz) : arr(a), n(sz) {}

  long long first() const { return arr[0]; }

  long long last() const { return arr[n - 1]; }

  Bucket split(long long p) {
    Bucket R(allocate_array(BSIZE), 0);
    long long* x = arr;
    long long* y = R.arr;
    long long* z = arr;
    for (int i = n; i--; z++) {
      int is_less = *z < p;
      *x = *z;
      x += is_less;
      *y = *z;
      y += !is_less;
    }
    n = int(x - arr);
    R.n = int(y - R.arr);
    return R;
  }

  bool is_less_than(long long p) const {
    for (int i = 0; i < n; i++) {
      if (arr[i] >= p) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) const {
    for (int i = 0; i < n; i++) {
      if (arr[i] < p) {
        return false;
      }
    }
    return true;
  }

  long long csum() const {
    long long ret = 0;
    for (int i = 0; i < n; i++) {
      ret += arr[i];
    }
    return ret;
  }
};

class Chain {
  vector<Bucket> buckets;
  int slack_index;
  int slack;

 public:
  Chain() : buckets(), slack_index(0), slack(0) {}

  void append(Bucket b) {
    slack += BSIZE - b.n;
    buckets.push_back(b);
  }

  vector<Bucket>& get_buckets() { return buckets; }

  int get_slack() const { return slack; }

  void compact() {
    sort(buckets.begin(), buckets.end(),
         [](const auto& a, const auto& b) { return a.n < b.n; });
    for (int i = 0; i + 1 < int(buckets.size());) {
      Bucket& L = buckets[i];
      Bucket& R = buckets.back();
      int m = min(BSIZE - L.n, R.n);
      memcpy(L.arr + L.n, R.arr + R.n - m, m * sizeof(long long));
      L.n += m;
      R.n -= m;
      if (L.n == BSIZE) i++;
      if (R.n == 0) buckets.pop_back();
    }
  }

  pair<long long*, int> next_slack() {
    while (slack_index < int(buckets.size())) {
      Bucket& b = buckets[slack_index];
      if (b.n < BSIZE) {
        return make_pair(b.arr + b.n, BSIZE - b.n);
      }
      slack_index++;
    }
    Bucket b(allocate_array(BSIZE), 0);
    append(b);
    return make_pair(b.arr, BSIZE);
  }

  void consume_slack(int amt) {
    assert(slack_index < int(buckets.size()));
    Bucket& b = buckets[slack_index];
    b.n += amt;
    slack -= amt;
    assert(b.n <= BSIZE);
  }

  bool is_less_than(long long p) const {
    for (Bucket b : buckets) {
      if (!b.is_less_than(p)) {
        return false;
      }
    }
    return true;
  }

  bool is_at_least(long long p) const {
    for (Bucket b : buckets) {
      if (!b.is_at_least(p)) {
        return false;
      }
    }
    return true;
  }

  int num_elements() const {
    int num = 0;
    for (Bucket b : buckets) {
      num += b.n;
    }
    return num;
  }

  long long csum() const {
    long long ret = 0;
    for (Bucket b : buckets) {
      ret += b.csum();
    }
    return ret;
  }
};

static pair<Chain, Chain> std_partition(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  for (Bucket b : chain.get_buckets()) {
    int n_flipped = 0;
    for (int i = 1; i < b.n && !n_flipped; i++) {
      n_flipped += b.arr[i - 1] > b.arr[i];
    }
    if (n_flipped == 0) {
      // is ascending.
      if (b.last() < p) {
        left_chain.append(b);
        continue;
      }
      if (b.first() >= p) {
        right_chain.append(b);
        continue;
      }
    }
    int i =
        int(partition(b.arr, b.arr + b.n, [p](auto const v) { return v < p; }) -
            b.arr);
    Bucket L(allocate_array(BSIZE), i);
    Bucket R(allocate_array(BSIZE), b.n - i);
    memcpy(L.arr, b.arr, L.n * sizeof(long long));
    memcpy(R.arr, b.arr + i, R.n * sizeof(long long));
    left_chain.append(L);
    right_chain.append(R);
  }
  return make_pair(left_chain, right_chain);
}

static pair<Chain, Chain> nobranch_partition(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  for (Bucket b : chain.get_buckets()) {
    int n_flipped = 0;
    for (int i = 1; i < b.n && !n_flipped; i++) {
      n_flipped += b.arr[i - 1] > b.arr[i];
    }
    if (n_flipped == 0) {
      // is ascending.
      if (b.last() < p) {
        left_chain.append(b);
        continue;
      }
      if (b.first() >= p) {
        right_chain.append(b);
        continue;
      }
    }

    Bucket L(allocate_array(BSIZE), 0);
    Bucket R(allocate_array(BSIZE), 0);

    long long* x = L.arr;
    long long* y = R.arr;
    long long* z = b.arr;

    for (long long* end = b.arr + b.n; z < end; z++) {
      int is_less = *z < p;
      *x = *z;
      x += is_less;
      *y = *z;
      y += !is_less;
    }
    L.n = int(x - L.arr);
    R.n = int(y - R.arr);

    left_chain.append(L);
    right_chain.append(R);
  }
  return make_pair(left_chain, right_chain);
}

static pair<Chain, Chain> nobranch_partition2(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  for (Bucket b : chain.get_buckets()) {
    int n_flipped = 0;
    for (int i = 1; i < b.n && !n_flipped; i++) {
      n_flipped += b.arr[i - 1] > b.arr[i];
    }
    if (n_flipped == 0) {
      // is ascending.
      if (b.last() < p) {
        left_chain.append(b);
        continue;
      }
      if (b.first() >= p) {
        right_chain.append(b);
        continue;
      }
    }
    right_chain.append(b.split(p));
    left_chain.append(b);
  }
  return make_pair(left_chain, right_chain);
}

static pair<Chain, Chain> nobranch_compact(Chain& chain, long long p) {
  auto res = nobranch_partition2(chain, p);
  res.first.compact();
  res.second.compact();
  return res;
}

static pair<Chain, Chain> nobranch_compact2(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  Random rng(1234);
  for (Bucket b : chain.get_buckets()) {
    int n_flipped = 0;
    for (int i = 1; i < b.n && !n_flipped; i++) {
      n_flipped += b.arr[i - 1] > b.arr[i];
    }
    if (n_flipped == 0) {
      // is ascending.
      if (b.last() < p) {
        left_chain.append(b);
        continue;
      }
      if (b.first() >= p) {
        right_chain.append(b);
        continue;
      }
    }

    int nhi = 0;
    for (int i = 0; i < 11; i++) {
      nhi += b.arr[rng.nextInt(b.n)] >= p;
    }

    // fprintf(stderr, "slack %d %d\n", left_chain.get_slack(),
    //         right_chain.get_slack());
    if (2 * nhi >= 11) {
      // fprintf(stderr, "left\n");
      long long* x = b.arr;
      long long* z = b.arr;
      while (z < b.arr + b.n) {
        auto s = right_chain.next_slack();
        int m = std::min(int(b.arr + b.n - z), s.second);
        long long* y = s.first;
        for (; m--; z++) {
          int is_less = *z < p;
          *x = *z;
          x += is_less;
          *y = *z;
          y += !is_less;
        }
        right_chain.consume_slack(int(y - s.first));
      }
      b.n = int(x - b.arr);
      left_chain.append(b);

      assert_dbg(left_chain.is_less_than(p));
      assert_dbg(right_chain.is_at_least(p));

    } else {
      // fprintf(stderr, "right\n");
      long long* x = b.arr;
      long long* z = b.arr;
      while (z < b.arr + b.n) {
        auto s = left_chain.next_slack();
        int m = std::min(int(b.arr + b.n - z), s.second);
        long long* y = s.first;
        for (; m--; z++) {
          int is_less = *z < p;
          *x = *z;
          x += !is_less;
          *y = *z;
          y += is_less;
        }
        left_chain.consume_slack(int(y - s.first));
      }
      b.n = int(x - b.arr);
      right_chain.append(b);

      assert_dbg(left_chain.is_less_than(p));
      assert_dbg(right_chain.is_at_least(p));
    }
  }
  return make_pair(left_chain, right_chain);
}

static void swap_offsets(long long* first, long long* last, unsigned char* hi,
                         unsigned char* lo, int num) {
  if (num > 0) {
    long long* l = first + hi[0];
    long long* r = last - lo[0];
    long long tmp(*l);
    *l = *r;
    for (int i = 1; i < num; ++i) {
      l = first + hi[i];
      *r = (*l);
      r = last - lo[i];
      *l = (*r);
    }
    *r = (tmp);
  }
}

#define block_size 64

long long* partition_right_branchless(long long* begin, long long* end,
                                      long long pivot) {
  long long* first = begin;
  long long* last = end;
  auto comp = std::less<long long>();
  while (first < last && comp(*first, pivot)) first++;
  while (first < last && !comp(*(last - 1), pivot)) last--;

  // The branchless partitioning is derived from "BlockQuicksort: How Branch
  // Mispredictions don’t affect Quicksort" by Stefan Edelkamp and Armin Weiss.
  unsigned char hi[block_size];
  unsigned char lo[block_size];
  int nhi = 0, nlo = 0, ihi = 0, ilo = 0;

  while (last - first > 2 * block_size) {
    if (nhi == 0) {
      ihi = 0;
      long long* it = first;
      for (unsigned char i = 0; i < block_size;) {
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
      }
    }
    if (nlo == 0) {
      ilo = 0;
      long long* it = last;
      for (unsigned char i = 0; i < block_size;) {
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
      }
    }

    // Swap elements and update block sizes and first/last boundaries.
    int num = std::min(nhi, nlo);
    swap_offsets(first, last, hi + ihi, lo + ilo, num);
    nhi -= num;
    nlo -= num;
    ihi += num;
    ilo += num;
    if (nhi == 0) first += block_size;
    if (nlo == 0) last -= block_size;
  }

  int l_size = 0, r_size = 0;
  int unknown_left = int((last - first) - ((nlo || nhi) ? block_size : 0));
  if (nlo) {
    // Handle leftover block by assigning the unknown elements to the other
    // block.
    l_size = unknown_left;
    r_size = block_size;
  } else if (nhi) {
    l_size = block_size;
    r_size = unknown_left;
  } else {
    // No leftover block, split the unknown elements in two blocks.
    l_size = unknown_left / 2;
    r_size = unknown_left - l_size;
  }

  // Fill offset buffers if needed.
  if (unknown_left && !nhi) {
    ihi = 0;
    long long* it = first;
    for (unsigned char i = 0; i < l_size;) {
      hi[nhi] = i++;
      nhi += !comp(*it++, pivot);
    }
  }
  if (unknown_left && !nlo) {
    ilo = 0;
    long long* it = last;
    for (unsigned char i = 0; i < r_size;) {
      lo[nlo] = ++i;
      nlo += comp(*--it, pivot);
    }
  }

  int num = std::min(nhi, nlo);
  swap_offsets(first, last, hi + ihi, lo + ilo, num);
  nhi -= num;
  nlo -= num;
  ihi += num;
  ilo += num;
  if (nhi == 0) first += l_size;
  if (nlo == 0) last -= r_size;

  // We have now fully identified [first, last)'s proper position. Swap the last
  // elements.
  if (nhi) {
    while (nhi--) std::iter_swap(first + hi[ihi + nhi], --last);
    first = last;
  }
  if (nlo) {
    while (nlo--) std::iter_swap(last - lo[ilo + nlo], first), ++first;
    last = first;
  }

  // Put the pivot in the right place.
  return first;
}

static pair<Chain, Chain> block_qsort_basic(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  for (Bucket b : chain.get_buckets()) {
    int n_flipped = 0;
    for (int i = 1; i < b.n && !n_flipped; i++) {
      n_flipped += b.arr[i - 1] > b.arr[i];
    }
    if (n_flipped == 0) {
      // is ascending.
      if (b.last() < p) {
        left_chain.append(b);
        continue;
      }
      if (b.first() >= p) {
        right_chain.append(b);
        continue;
      }
    }
    int i = int(partition_right_branchless(b.arr, b.arr + b.n, p) - b.arr);
    Bucket L(allocate_array(BSIZE), i);
    Bucket R(allocate_array(BSIZE), b.n - i);
    memcpy(L.arr, b.arr, L.n * sizeof(long long));
    memcpy(R.arr, b.arr + i, R.n * sizeof(long long));
    left_chain.append(L);
    right_chain.append(R);
  }
  return make_pair(left_chain, right_chain);
}

static pair<int, int> partition_branchless(long long* L, int nL, long long* R,
                                           int nR, long long pivot) {
  // fprintf(stderr, "n = %d %d\n", nL, nR);
  long long* first = L;
  long long* last = R + nR;
  auto comp = std::less<long long>();

  // TODO: Unroll.
  while (first < L + nL && comp(*first, pivot)) first++;
  while (R < last && !comp(*(last - 1), pivot)) last--;

  // The branchless partitioning is derived from "BlockQuicksort: How Branch
  // Mispredictions don’t affect Quicksort" by Stefan Edelkamp and Armin Weiss.
  unsigned char hi[block_size];
  unsigned char lo[block_size];
  int nhi = 0, nlo = 0, ihi = 0, ilo = 0;

  while (first + block_size <= L + nL && last - block_size >= R) {
    if (nhi == 0) {
      ihi = 0;
      long long* it = first;
      for (unsigned char i = 0; i < block_size;) {
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
        hi[nhi] = i++;
        nhi += !comp(*it++, pivot);
      }
    }
    if (nlo == 0) {
      ilo = 0;
      long long* it = last;
      for (unsigned char i = 0; i < block_size;) {
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
        lo[nlo] = ++i;
        nlo += comp(*--it, pivot);
      }
    }

    // Swap elements and update block sizes and first/last boundaries.
    int num = std::min(nhi, nlo);
    swap_offsets(first, last, hi + ihi, lo + ilo, num);
    nhi -= num;
    nlo -= num;
    ihi += num;
    ilo += num;
    if (nhi == 0) first += block_size;
    if (nlo == 0) last -= block_size;
  }

  int l_size = min(block_size, int(L + nL - first));
  if (first < L + nL && !nhi) {
    ihi = 0;
    long long* it = first;
    for (unsigned char i = 0; i < l_size;) {
      hi[nhi] = i++;
      nhi += !comp(*it++, pivot);
    }
  }

  int r_size = min(block_size, int(last - R));
  if (last > R && !nlo) {
    ilo = 0;
    long long* it = last;
    for (unsigned char i = 0; i < r_size;) {
      lo[nlo] = ++i;
      nlo += comp(*--it, pivot);
    }
  }

  // fprintf(stderr, "rem %d %d\n", l_size, r_size);
  int num = std::min(nhi, nlo);
  swap_offsets(first, last, hi + ihi, lo + ilo, num);
  nhi -= num;
  nlo -= num;
  ihi += num;
  ilo += num;
  if (nhi == 0) first += l_size;
  if (nlo == 0) last -= r_size;
  assert(first <= L + nL);
  assert(last >= R);

  return make_pair(first - L, R + nR - last);
}

static pair<Chain, Chain> block_qsort_opt(Chain& chain, long long p) {
  Chain left_chain, right_chain;
  vector<Bucket>& buckets = chain.get_buckets();
  int i = 0, j = int(buckets.size()) - 1;
  Bucket* L = nullptr;
  Bucket* R = nullptr;
  int ihi = 0, ilo = 0;
  while (i <= j) {
    // fprintf(stderr, "%d %d\n", i, j);
    if (L == nullptr) {
      L = &buckets[i++];
      ihi = 0;
    } else if (R == nullptr) {
      R = &buckets[j--];
      ilo = 0;
    } else {
      auto m =
          partition_branchless(L->arr + ihi, L->n - ihi, R->arr, R->n - ilo, p);
      // fprintf(stderr, "m %d %d\n", m.first, m.second);
      ihi += m.first;
      assert(ihi <= L->n);
      ilo += m.second;
      assert(ilo <= R->n);
      if (ihi == L->n) {
        assert_dbg(left_chain.is_less_than(p));
        left_chain.append(*L);
        assert_dbg(left_chain.is_less_than(p));
        L = nullptr;
      }
      if (ilo == R->n) {
        assert_dbg(right_chain.is_at_least(p));
        right_chain.append(*R);
        assert_dbg(right_chain.is_at_least(p));
        R = nullptr;
      }
    }
  }
  if (L) {
    right_chain.append(L->split(p));
    assert_dbg(right_chain.is_at_least(p));
    left_chain.append(*L);
    assert_dbg(left_chain.is_less_than(p));
  }
  if (R) {
    right_chain.append(R->split(p));
    assert_dbg(right_chain.is_at_least(p));
    left_chain.append(*R);
    assert_dbg(left_chain.is_less_than(p));
  }
  return make_pair(left_chain, right_chain);
}

#define MAXN (1 << 26)
#define NBUCKETS (1 << 10)

static long long arr[MAXN];
static int arr_sz;

static long long* allocate_array(int sz) {
  arr_sz += sz;
  return arr + arr_sz - sz;
}

static void varying_randomness(function<void(int)> cb) {
  for (int percent_e3 = 100000, m = 1;;) {
    percent_e3 = max(0, percent_e3);
    cb(percent_e3);
    if (percent_e3 <= 0) break;
    percent_e3 -= m;
    m *= 2;
  }
}

static void run_test(const char* name,
                     function<pair<Chain, Chain>(Chain&, long long)> algo) {
  fprintf(stderr, "%20s: ", name);
  Random rng(140384);
  double total_time = 0;
  long long chk = 0;
  long long nleft = 0, nright = 0;
  varying_randomness([&](int percent_e3) {
    Chain chain;
    arr_sz = 0;
    for (int i = 0; i < NBUCKETS; i++) {
      chain.append(Bucket(allocate_array(BSIZE), BSIZE));
    }

    percent_e3 = max(0, percent_e3);
    for (int i = 0; i < arr_sz; i++) {
      arr[i] = i;  // rng.nextLong();
      if (rng.nextInt(100000) < 100000 - percent_e3) {
        swap(arr[i], arr[rng.nextInt(i + 1)]);
      }
    }

    long long p = arr_sz / 2;  // Roughly in the middle.
    pair<Chain, Chain> res;
    double split_time = time_it([&] { res = algo(chain, p); });
    nleft += res.first.num_elements();
    nright += res.second.num_elements();
    chk = (chk * 13 + res.first.csum()) * 7 + res.second.csum();

    fprintf(stderr, "%6.2lf", split_time * 1000);
    total_time += split_time;
  });

  fprintf(stderr, " %10.3lf %5.2lf %20lld\n", total_time * 1000,
          100.0 * nleft / (nleft + nright), chk);
}

int main(int argc, char* argv[]) {
  memset(arr, 1, sizeof(long long) * MAXN);
  fprintf(stderr, "%20s  ", "");
  varying_randomness(
      [](int percent_e3) { fprintf(stderr, "%6.2lf", percent_e3 / 1000.0); });
  fprintf(stderr, " %10s ratio %20s\n", "total_time", "csum");

  run_test("block_qsort_opt", block_qsort_opt);
  run_test("std_partition", std_partition);
  run_test("nobranch_partition", nobranch_partition);
  run_test("nobranch_partition2", nobranch_partition2);
  run_test("nobranch_compact", nobranch_compact);
  run_test("nobranch_compact2", nobranch_compact2);
  run_test("block_qsort_basic", block_qsort_basic);
}
