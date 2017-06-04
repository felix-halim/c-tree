#ifndef CRACK_MEMORY_MANAGER_H_
#define CRACK_MEMORY_MANAGER_H_

#include <cassert>
#include <cstring>

#include <memory>
#include <functional>
#include <algorithm>

namespace trimmer {

template<typename T, template<class> class alloc>
class memory_manager {

  alloc<int> ints_allocator;
  alloc<T> arr_allocator;

  int* ints;
  T* arr;

  bool is_zero(T* p) {
    char* b = (char*) p;
    for (unsigned i = 0; i < sizeof(T); i++)
      if (b[i]) return false;
    return true;
  }

  int* free_indices() { return ints + 2; }
  int free_indices_len() const { return ints[0]; }
  void set_free_indices_len(int len) { ints[0] = len; }

  const int& arr_len() const { return ints[1]; }
  void set_arr_len(int len) { ints[1] = len; }

  void free_indices_push(int idx) {
    int* fidxs = free_indices();
    int len = free_indices_len();
    fidxs[len++] = idx;
    std::push_heap(fidxs, fidxs + len, std::greater<int>());
    set_free_indices_len(len);
  }
  int free_indices_pop() {
    int* fidxs = free_indices();
    int len = free_indices_len();
    std::pop_heap(fidxs, fidxs + len, std::greater<int>());
    set_free_indices_len(--len);
    return fidxs[len];
  }
  bool free_indices_non_allocated_is_zeroed() {
    int* fidxs = free_indices();
    int len = free_indices_len();
    for (int i = 0; i < len; i++)
      if (!is_zero(&arr[fidxs[i]])) {
        fprintf(stderr, "arr[%d] is non zero\n", fidxs[i]);
        return false;
      }
    return true;
  }


 public:

  memory_manager():
    ints(ints_allocator.allocate(2)),
    arr(nullptr) {
      set_free_indices_len(0);
      set_arr_len(0);
  }

  int capacity() const {
    return arr_len();
  }

  int size() const {
    return arr_len() - free_indices_len();
  }

  int index_of(T* p) const {
    int index = int(p - arr);
    assert(index >= 0 && index < arr_len());
    return index;
  }

  T* get(int index) const {
    assert(index >= 0 && index < arr_len());
    return &arr[index];
  }

  void reserve(int nlen) {
    int len = arr_len();
    ints = ints_allocator.reallocate(ints, 2 + len, 2 + nlen);
    arr = arr_allocator.reallocate(arr, len, nlen);
    for (int i = len; i < nlen; i++) free_indices_push(i);
    memset(arr + len, 0, sizeof(T) * (nlen - len));
    set_arr_len(nlen);
  }

  T* create() {
    // fprintf(stderr, "create(), free = %d\n", free_indices_len());
    if (free_indices_len() == 0) {
      // Grow arr.
      if (arr_len() == 0) {
        ints = ints_allocator.reallocate(ints, 2, 3);
        arr = arr_allocator.allocate(1);
        free_indices_push(0);
        set_arr_len(1);
        memset(arr, 0, sizeof(T));
      } else {
        fprintf(stderr, "doubling %d -> %d\n", arr_len(), 2 * arr_len());
        reserve(arr_len() * 2);
      }
    }
    assert(free_indices_len() > 0);
    int index = free_indices_pop();
    // fprintf(stderr, "created arr index = %d, %p\n", index, &arr[index]);
    assert(is_zero(&arr[index]));
    return &arr[index];
  }

  void destroy(int p) {
    destroy(get(p));
  }

  void destroy(T* p) {
    p->~T();
    memset(p, 0, sizeof(T));
    // fprintf(stderr, "destroying %p %d\n", p, index_of(p));
    free_indices_push(index_of(p));
  }

  long long size_in_bytes() {
    return ((long long) sizeof(T)) * arr_len() + sizeof(int) * (arr_len() + 2LL);
  }

  void print_stats() {
    fprintf(stderr, "len = %5d / %5d, size = %10lld\n", size(), arr_len(), size_in_bytes());
  }

  void check() {
    assert(free_indices_non_allocated_is_zeroed());
  }
};

};

#endif
