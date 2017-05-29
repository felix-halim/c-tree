#ifndef CRACK_UTIL_H_
#define CRACK_UTIL_H_

#include <random>

namespace trimmer {

using std::mt19937;
using std::uniform_int_distribution;

constexpr int INTERNAL_NODE_BIT = 1 << 30;
constexpr int EMPTY_NODE = -1;

static bool is_internal(int p) {
  return p & INTERNAL_NODE_BIT;
}

template<typename T, typename compare> bool eq(T const &a, T const &b, compare &lt) { return !lt(a, b) && !lt(b, a); }
template<typename T, typename compare> bool lt(T const &a, T const &b, compare &lt) { return lt(a, b); }
template<typename T, typename compare> bool lte(T const &a, T const &b, compare &lt) { return !lt(b, a); }
template<typename T, typename compare> bool gt(T const &a, T const &b, compare &lt) { return lt(b, a); }
template<typename T, typename compare> bool gte(T const &a, T const &b, compare &lt) { return !lt(a, b); }

template<typename T, typename compare>
static int sorted_lower_pos(T const data[], int n, T const &value, compare &lt) {
  int pos = 0;
  while (pos < n && lt(data[pos], value)) pos++;
  return pos;
}

class Random {
  mt19937 gen;
  uniform_int_distribution<> dis;

 public:
  Random() : gen(140384) {}
  Random(int seed) : gen(seed) {}
  int nextInt() { return dis(gen); }
  int nextInt(int N) { return dis(gen) % N; } // Poor I know.
};

template<typename internal_type, typename leaf_type, template<class> class alloc>
class node_manager {
  memory_manager<internal_type, alloc> internal_mm_;
  memory_manager<leaf_type, alloc> leaf_mm_;

 public:

  internal_type* create_internal() {
    return internal_mm_.create();
  }

  leaf_type* create_leaf() {
    return leaf_mm_.create();
  }

  int next_power_of_two(int v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
  }

  void reserve_leaf(int n) {
    assert(n > 0);
    int nleaves = (leaf_type::capacity - 1 + n) / leaf_type::capacity;
    leaf_mm_.reserve(next_power_of_two(nleaves));
  }

  void destroy_internal(int p) {
    assert(is_internal(p));
    internal_mm_.destroy(p ^ INTERNAL_NODE_BIT);
  }

  void destroy_leaf(int p) {
    assert(!is_internal(p));
    leaf_mm_.destroy(p);
  }

  int index_of(internal_type* p) const {
    return internal_mm_.index_of(p) | INTERNAL_NODE_BIT;
  }

  int index_of(leaf_type* p) const {
    return leaf_mm_.index_of(p);
  }

  internal_type* get_internal(int p) const {
    assert(is_internal(p));
    return internal_mm_.get(p ^ INTERNAL_NODE_BIT);
  }

  leaf_type* get_leaf(int p) const {
    assert(!is_internal(p));
    return leaf_mm_.get(p);
  }

  void print_stats() {
    fprintf(stderr, "internal: "); internal_mm_.print_stats();
    fprintf(stderr, "leaf    : "); leaf_mm_.print_stats();
  }
};

}

#endif
