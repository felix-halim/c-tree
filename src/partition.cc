/*
To run, go to the root folder and execute:

make partition
*/

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <vector>

#include "random.h"
#include "time_it.h"

#define BSIZE (1 << 12)
#define block_size 64

static long long arr[BSIZE];

int std_partition(long long p) {
  return int(
      std::partition(arr, arr + BSIZE, [p](long long v) { return v < p; }) -
      arr);
}

template <typename RandomAccessIterator>
void swap_offsets(RandomAccessIterator first, RandomAccessIterator last,
                  unsigned char* hi, unsigned char* lo, int num,
                  bool use_swaps) {
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
  if (use_swaps) {
    // This case is needed for the descending distribution, where we need
    // to have proper swapping for pdqsort to remain O(n).
    for (int i = 0; i < num; ++i) {
      std::iter_swap(first + hi[i], last - lo[i]);
    }
  } else if (num > 0) {
    RandomAccessIterator l = first + hi[0];
    RandomAccessIterator r = last - lo[0];
    T tmp(*l);
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
    swap_offsets(first, last, hi + ihi, lo + ilo, num, nhi == nlo);
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
  swap_offsets(first, last, hi + ihi, lo + ilo, num, nhi == nlo);
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

int fusion(long long p) {
  return int(partition_right_branchless(arr, arr + BSIZE, p) - arr);
}

int lo[BSIZE], hi[BSIZE];
int marker(long long p) {
  int nhi = 0, nlo = 0;
  for (int i = 0; i < BSIZE; i++) {
    hi[nhi] = i;
    nhi += (arr[i] >= p);
  }

  for (int i = 0; i < BSIZE; i++) {
    lo[nlo] = i;
    nlo += (arr[i] < p);
  }

  int* hip = hi;
  int* lop = lo + nlo - 1;
 
  int m = std::min(nhi, nlo);
  while (m-- && *hip < *lop) {
    std::swap(arr[*(hip++)], arr[*(lop--)]);
  }

  return nlo;
}

int main() {
  Random rng(140384);

  double total_time = 0;
  long long chk = 0;
  for (int m = 0; m < 200000; m++) {
    for (int i = 0; i < BSIZE; i++) {
      arr[i] = rng.nextLong();
    }

    long long p = arr[rng.nextInt(BSIZE)];
    int pos = -1;
    total_time += time_it([&]() {
      // pos = std_partition(p);
      pos = fusion(p);
      // pos = marker(p);
    });

    // fprintf(stderr, "pos = %d\n", pos);
    for (int i = 0; i < BSIZE; i++) {
      // fprintf(stderr, "%4d %4d %20lld, %20lld\n", i, pos, arr[i], p);
      assert((i < pos) ? (arr[i] < p) : (arr[i] >= p));
    }

    chk = chk * 13 + pos;
  }

  fprintf(stderr, "tt = %.6lf, chk = %lld\n", total_time, chk);
  assert(chk == 2144819404491265167LL);
}
