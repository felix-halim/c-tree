#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <stack>
#include <set>
#include <set>
#include <vector>
#include <chrono>

#include "update.h"

using namespace std;
using namespace std::chrono;

#define REP(i, n) for (int i = 0, _n = n; i < _n; i++)

template<typename Func>
double time_it(Func f) {
  auto t0 = high_resolution_clock::now();
  f();
  auto t1 = high_resolution_clock::now();
  return duration_cast<microseconds>(t1 - t0).count() * 1e-6;
}

struct Bucket {
  int arr[BSIZE], n;
  Bucket *next;
  Bucket(): n(0), next(0) {}
  int slack() { return BSIZE - n; }
};

pair<int, int> count_lower(Bucket *B, int P) {
  int cnt = 0, sum = 0;
  while (B) {
    REP(i, B->n) if (B->arr[i] < P) cnt++, sum += B->arr[i];
    B = B->next;
  }
  return make_pair(cnt, sum);
}

// Count only using std::partition, without actually moving the data.
Bucket* hypothetical_partition(Bucket *B, int P) {
  for (; B; B = B->next) {
    int pos = partition(B->arr, B->arr + B->n, bind2nd(less<int>(), P)) - B->arr;
    assert(pos >= 0 && pos <= B->n);
  }
  return NULL;
}

// Partitions the buckets using if-else branch.
Bucket* use_if(Bucket *B, int P) {
  Bucket *F = NULL; // Free-ed bucket.
  Bucket *L = new Bucket(), *Lh = L;
  Bucket *R = new Bucket();
  while (B) {
    Bucket *next = B->next;
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
      if (!L->slack()) { assert(F && !F->n); L->next = F; L = F; F = NULL; L->next = NULL; }
      if (!R->slack()) { assert(F && !F->n); R->next = F; R = F; F = NULL; }
    }
    B = next;
  }
  return Lh;
}

/*
int use_par(Bucket *B, int P) {
  Bucket *L = new Bucket(BSIZE);
  Bucket *R = new Bucket(BSIZE);
  set<Bucket*> F;
  int nL = 0;
  while (B) {
    int pos = partition(B->arr, B->arr + B->n, bind2nd(less<int>(), P)) - B->arr;
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
      L = i;
      nL += BSIZE;
    }

    if (R->slack() >= BSIZE - pos){
      memmove(R->arr + R->n, B->arr + pos, sizeof(int)*(BSIZE-pos));
      R->n += BSIZE-pos;
      if (L != i){
        F.insert(i);
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
  nL += L->n;
  return nL;
}

int use_paro(Bucket *B, int P) {
  int t1=0,t2=0,t3=0,t4=0;
  vector<pair<int,int> > L, R;
  REP(i,NBUCKETS){
    int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P)) - B[i].arr;
    if (pos*2 > B[i].n) L.push_back(make_pair(pos, i));
    else R.push_back(make_pair(pos, i));
  }

  // t1 += timit();

  // fprintf(stderr, "%d %d <<\n", L.size(), R.size());
  int nL = 0;
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
  // fprintf(stderr, "%d %d <<\n", L.size(), R.size());

  // t2 += timit();

  for (int k=0; k < (int)L.size(); ){
    if (k+1 == (int)L.size()){ nL += L.back().first; break; }
    pair<int,int> &pL = L[k], &pR = L.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(B[pL.second].n - i, j);
    while (m--) swap(a[i++], b[--j]);
    if (i>=B[pL.second].n){ k++; nL += BSIZE; }
    if (j==0) L.pop_back();
  }

  // t3 += timit();

  for (int k=0; k < (int)R.size(); ){
    if (k+1 == (int)R.size()){ nL += R.back().first; break; }
    pair<int,int> &pL = R[k], &pR = R.back();
    int *a = B[pL.second].arr, &i = pL.first;
    int *b = B[pR.second].arr, &j = pR.first;
    int m = min(i, B[pR.second].n - j);
    // fprintf(stderr, ">> %d >> %d %d, m = %d; %d %d\n", R.size(),i,j,m,pL.second,pR.second);
    while (m--) swap(a[--i], b[j++]);
    // fprintf(stderr, ">>> %d %d %d, m = %d\n", i,j,B[pR.second].n,m);
    if (j>=B[pR.second].n){ R.pop_back(); nL += BSIZE; }
    if (i==0){ k++;  }
  }

  // t4 += timit();

  // fprintf(stderr, "%d == %d %.4lf %.4lf %.4lf %.4lf, %d %d\n", count_lower(P), nL, t1,t2,t3,t4,L.size(),R.size());
  return nL;
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
      // fprintf(stderr,"i = %d, j = %d/%d, R = %d, %d %d, %d\n",i,j,B[i].n,R,cnt[0],B[i].n,B[R].slack());
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
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) + (B[i].arr[j+3] < P);
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
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) + (B[i].arr[j+3] < P) +
          (B[i].arr[j+4] < P) + (B[i].arr[j+5] < P) + (B[i].arr[j+6] < P) + (B[i].arr[j+7] < P);
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
      nL += (B[i].arr[j] < P) + (B[i].arr[j+1] < P) + (B[i].arr[j+2] < P) + (B[i].arr[j+3] < P) +
          (B[i].arr[j+4] < P) + (B[i].arr[j+5] < P) + (B[i].arr[j+6] < P) + (B[i].arr[j+7] < P) +
          (B[i].arr[j+8] < P) + (B[i].arr[j+9] < P) + (B[i].arr[j+10] < P) + (B[i].arr[j+11] < P) +
          (B[i].arr[j+12] < P) + (B[i].arr[j+13] < P) + (B[i].arr[j+14] < P) + (B[i].arr[j+15] < P);
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
      char c1 = ((B[i].arr[j] < P) << 3) | ((B[i].arr[j+1] < P) << 2) | ((B[i].arr[j+2] < P) << 1) | ((B[i].arr[j+3] < P) << 0);
      char c2 = ((B[i].arr[j+4] < P) << 3) | ((B[i].arr[j+5] < P) << 2) | ((B[i].arr[j+6] < P) << 1) | ((B[i].arr[j+7] < P) << 0);
      char c3 = ((B[i].arr[j+8] < P) << 3) | ((B[i].arr[j+9] < P) << 2) | ((B[i].arr[j+10] < P) << 1) | ((B[i].arr[j+11] < P) << 0);
      char c4 = ((B[i].arr[j+12] < P) << 3) | ((B[i].arr[j+13] < P) << 2) | ((B[i].arr[j+14] < P) << 1) | ((B[i].arr[j+15] < P) << 0);
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
  int nL = 0, *a[4] = { B[NBUCKETS].arr, B[NBUCKETS+1].arr, B[NBUCKETS+2].arr, B[NBUCKETS+3].arr};
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
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P)) - B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else R.push_back(make_pair(pos, i));
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
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P)) - B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else R.push_back(make_pair(pos, i));
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
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P)) - B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else R.push_back(make_pair(pos, i));
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
      int pos = partition(B[i].arr, B[i].arr + B[i].n, bind2nd(less<int>(), P)) - B[i].arr;
      if (pos*2 > B[i].n) L.push_back(make_pair(pos, i)); else R.push_back(make_pair(pos, i));
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
      if (L1==-1){ assert(F1!=-1); L1 = F1; F1 = -1; } else if (R1==-1){ assert(F1!=-1); R1 = F1; F1 = -1; }
      if (L2==-1){ assert(F2!=-1); L2 = F2; F2 = -1; } else if (R2==-1){ assert(F2!=-1); R2 = F2; F2 = -1; }
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
      if (L1==-1){ assert(F1!=-1); L1 = F1; F1 = -1; } else if (R1==-1){ assert(F1!=-1); R1 = F1; F1 = -1; }
      if (L2==-1){ assert(F2!=-1); L2 = F2; F2 = -1; } else if (R2==-1){ assert(F2!=-1); R2 = F2; F2 = -1; }
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
const char **uw = update_workload;

int main(int argc, char *argv[]) {
  Update update(argv[1], 0);
  update.load();
  int N = update.size() / 2;

  vector<int> samples;
  int algo, ith, *arr = update.get_arr();
  REP(i, 10000) samples.push_back(arr[rand() % N]);
  sort(samples.begin(), samples.end());

  sscanf(argv[2],"%d",&algo);
  sscanf(argv[3],"%d",&ith);
  Bucket *B = NULL, *T = NULL;
  int NBUCKETS = 0;

  double insert_time = time_it([&] {
    int i = 0;
    for (; i + BSIZE <= N; NBUCKETS++) {
      Bucket *b = new Bucket();
      REP(j, BSIZE) b->arr[b->n++] = arr[i++];
      if (B) T->next = b; else B = T = b;
      T = b;
    }
    if (i < N) {
      Bucket *b = new Bucket();
      REP(j, N - i) b->arr[b->n++] = arr[i++];
      if (B) T->next = b; else B = T = b;
      T = b;
      NBUCKETS++;
    }
  });

  int P = samples[(ith * samples.size()) / 100];
  pair<int, int> cs1 = count_lower(B, P);
  Bucket *smaller = NULL;
  double split_time = time_it([&] {
    switch (algo) {
      case 1 : smaller = hypothetical_partition(B, P); break;
      case 2 : smaller = use_if(B, P); break;
      /*
      case 3 : smaller = use_par(B, P); break;
      case 4 : smaller = use_nb(B, P); break;
      case 5 : smaller = use_nbp(B, P); break;
      case 6 : smaller = use_nbo(B, P); break;
      case 8 : smaller = use_cnt4(B, P); break;
      case 9 : smaller = use_cpy(B, P); break;
      case 10 : smaller = use_cpyp2(B, P); break;
      case 11 : smaller = use_cpyp4(B, P); break;
      case 12 : smaller = use_cpyp8(B, P); break;
      case 13 : smaller = use_cpyp16(B, P); break;
      case 14 : smaller = use_memo(B, P); break;
      case 15 : smaller = use_cpy4(B, P); break;
      case 16 : smaller = use_noswap(B, P); break;
      case 17 : smaller = use_nbswap(B, P); break;
      case 18 : smaller = use_nbswap2(B, P); break;
      case 19 : smaller = use_nbsortswap(B, P); break;
      case 20 : smaller = use_nbsortswap2(B, P); break;
      case 21 : smaller = use_partest(B, P); break;
      case 22 : smaller = use_paro(B, P); break;
      case 23 : smaller = use_parnb(B, P); break;
      */
    }
  });
  int cnt = 0, sum = 0;
  while (smaller) {
    cnt += smaller->n;
    // assert(smaller->n > 0 || !smaller->next);
    while (smaller->n--) {
      sum += smaller->arr[smaller->n];
    }
    smaller = smaller->next;
  }
  assert(cs1.first == cnt);
  assert(cs1.second == sum);

  printf("%d,%d,%d,%d,%d,%.6lf,%.6lf,%d,%d\n", N, algo, BSIZE, NBUCKETS, ith, insert_time, split_time, cnt, sum);
}
