#ifndef CRACK_MULTISET_H_
#define CRACK_MULTISET_H_

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "memory_manager.h"
#include "mallocator.h"
#include "multiset_util.h"
#include "multiset_node.h"
#include "random.h"

namespace trimmer {

using std::less;
using std::pair;

template <typename T, typename compare, template<class> class alloc = std::allocator>
class multiset {
 public:
  typedef T       key_type;
  typedef T       value_type;
  typedef compare key_compare;
  typedef compare value_compare;
  // typedef alloc   allocator_type;
  typedef int     size_type;

 private:
  typedef internal_node<T, compare> internal_type;
  typedef leaf_node<T, compare> leaf_type;

  node_manager<internal_type, leaf_type, alloc> nmgr;
  value_compare lt_;
  int root_;

  internal_type* get_internal(int p) const {
    assert(is_internal(p));
    return nmgr.get_internal(p);
  }

  leaf_type* get_leaf(int p) const {
    assert(!is_internal(p));
    return nmgr.get_leaf(p);
  }

  int rec(int b, std::function<int(int)> cnt) const {
    if (b == EMPTY_NODE) return 0;
    if (is_internal(b)) {
      internal_type* p = get_internal(b);
      int n = p->size();
      int ret = cnt(b);
      for (int i = 0; i <= n; i++) {
        ret += rec(p->child(i), cnt);
      }
      return ret;
    }
    int ret = 0;
    while (b != EMPTY_NODE) {
      leaf_type* p = get_leaf(b);
      ret += cnt(b);
      b = p->next();
    }
    return ret;
  }

 public:

  multiset(): root_(EMPTY_NODE) {}

  multiset(value_compare &lt): root_(EMPTY_NODE), lt_(lt) {}

  ~multiset() {}

  multiset& operator= (const multiset& x);

  class iterator {
   public:
    multiset<T, compare, alloc> *m_;
    int b_;      // Bucket index.
    int pos_;    // Position in the bucket.

    iterator(multiset<T, compare, alloc> *m, int b, int pos):
      m_(m),
      b_(b),
      pos_(pos) {}

    const value_type& operator*() {
      // TODO: remove.
      assert(b_ != EMPTY_NODE);
      if (is_internal(b_)) {
        return m_->get_internal(b_)->data(pos_);
      }
      return m_->get_leaf(b_)->data(pos_);
    }

    const value_type* operator->() {
      return &operator*();
    }

    iterator operator++(int amt) {
      assert(0);
    }

    iterator operator++() {
      assert(0);
    }

    iterator operator--(int amt);
    iterator operator--();

    bool operator==(iterator const &that) const {
      if (m_ != that.m_ || b_ != that.b_) return false;
      if (b_ == EMPTY_NODE) return true;
      return pos_ == that.pos_;
    }
  };

  class const_iterator {
    
  };

  iterator begin();
  const_iterator begin() const;

  iterator end() {
    return iterator(this, EMPTY_NODE, EMPTY_NODE);
  }

  const_iterator end() const;

  bool empty() const;

  size_type size() const {
    return rec(root_, [&](int b) { return is_internal(b) ? this->get_internal(b)->size() : this->get_leaf(b)->size(); });
  }

  size_type max_size() const;


  void print_stats() {
    fprintf(stderr, "\nstats:\n");
    int islack = 0;
    int isize = rec(root_, [&](int b){
      if (!is_internal(b)) return 0;
      int s = this->get_internal(b)->size();
      islack += internal_type::capacity - s;
      return s;
    });
    fprintf(stderr, "internal size = %6d, slack = %d\n", isize, islack);

    int lslack = 0;
    int lsize = rec(root_, [&](int b){
      if (is_internal(b)) return 0;
      int s = this->get_leaf(b)->size();
      lslack += leaf_type::capacity - s;
      return s;
    });
    fprintf(stderr, "leaf     size = %6d, slack = %d\n", lsize, lslack);

    fprintf(stderr, "memory_manager:\n");
    nmgr.print_stats();
  }

  // Inserts {value} to node {b}, updating the pointer to head via callback {set_head}.
  void insert(int b, const value_type& value, std::function<void(int)> set_head) {
    if (b == EMPTY_NODE) {
      b = new_leaf_node(EMPTY_NODE);
      get_leaf(b)->set_next(EMPTY_NODE);
      get_leaf(b)->append(value);
      set_head(b);
    } else if (!get_leaf(b)->is_full()) {
      get_leaf(b)->append(value);
    } else {
      int next = get_leaf(b)->next();
      if (next != EMPTY_NODE && !get_leaf(next)->is_full()) {
        get_leaf(next)->append(value);
      } else {
        int nb = new_leaf_node(get_leaf(b)->parent());
        get_leaf(nb)->append(value);
        get_leaf(nb)->set_next(next);
        get_leaf(b)->set_next(nb);
      }
    }
  }

  // iterator insert(const value_type& value) {
  void insert(const value_type& value) {
    // fprintf(stderr, "ins %d\n", value.value);
    if (root_ != EMPTY_NODE && is_internal(root_)) {
      int b = get_internal(root_)->next();
      insert(b, value, [&](int head) { get_internal(root_)->set_next(head); });
    } else {
      insert(root_, value, [&](int root) { root_ = root; });
    }
    // fprintf(stderr, "done ins %d\n", value.value);
  }

  iterator insert(iterator position, const value_type& val);

  void insert(T* first, T* last) {
    // fprintf(stderr, "batch %d\n", N);
    int i = 0, N = last - first;
    nmgr.reserve_leaf(N);
    while (i + leaf_type::capacity <= N) {
      leaf_type* p = nmgr.create_leaf();
      p->init(EMPTY_NODE);
      p->append(first + i, leaf_type::capacity);
      p->set_next(root_);
      root_ = nmgr.index_of(p);
      i += leaf_type::capacity;
      // fprintf(stderr, "chain %d, %d\n", i, b->size());
    }
    // fprintf(stderr, "done %d %d\n", i, N);
    while (i < N) {
      insert(first[i++]);
    }
    // fprintf(stderr, "inserted %d elements\n", size());
  }

  void try_destroy_internal_node(int b) {
    if (root_ == b) return;
    internal_type* p = get_internal(b);
    int parent = p->parent();
    internal_type* ip = get_internal(parent);
    int cpos = ip->child_pos(b);
    if (cpos > 0 && (p->size() + get_internal(ip->child(cpos - 1))->size() + 2) <= internal_type::capacity) {
      // fprintf(stderr, "destroy_internal1 %d\n", b);
      int nb = ip->child(cpos - 1);
      internal_type* np = get_internal(nb);
      update_child_parents(b, nb);
      p->move_to_lesser(np, ip->data(cpos - 1), lt_);
      ip->erase_pos(cpos - 1, 1);
      nmgr.destroy_internal(b);
      try_destroy_internal_node(parent);
    } else if (cpos < ip->size() && (p->size() + get_internal(ip->child(cpos + 1))->size() + 2) <= internal_type::capacity) {
      // fprintf(stderr, "destroy_internal2 %d\n", b);
      int nb = ip->child(cpos + 1);
      internal_type* np = get_internal(nb);
      update_child_parents(nb, b);
      np->move_to_lesser(p, ip->data(cpos), lt_);
      ip->erase_pos(cpos, 1);
      nmgr.destroy_internal(nb);
      try_destroy_internal_node(parent);
    } else if (root_ == parent && !ip->size()) {
      // fprintf(stderr, "destroy_internal_root %d\n", b);
      nmgr.destroy_internal(parent);
      root_ = b;
      get_internal(root_)->set_parent(EMPTY_NODE);
    }
  }

  void destroy_leaf_node(int b) {
    assert(b != EMPTY_NODE);
    assert(!is_internal(b));
    assert(!get_leaf(b)->size());

    if (b == root_) {
      nmgr.destroy_leaf(b);
      root_ = EMPTY_NODE;
      return;
    }

    leaf_type* p = get_leaf(b);
    int parent = p->parent();
    assert(parent != EMPTY_NODE);

    internal_type* ip = get_internal(parent);
    int cpos = ip->child_pos(b);
    if (cpos > 0) {
      // fprintf(stderr, "destroy_leaf %d, cpos = %d / %d\n", b, cpos, ip->size());
      leaf_type* np = get_leaf(ip->child(cpos - 1));
      if (np->is_full()) {
        // fprintf(stderr, "full %d\n", cpos - 1);
        p->append(ip->data(cpos - 1));
        ip->set_data(cpos - 1, np->remove_median(lt_));
        np->move_data_at_least(ip->data(cpos - 1), p, lt_);
      } else {
        // fprintf(stderr, "has room %d\n", cpos - 1);
        np->append(ip->data(cpos - 1));
        ip->erase_pos(cpos - 1, 1);
        nmgr.destroy_leaf(b);
        try_destroy_internal_node(parent);
      }
    } else if (cpos < ip->size()) {
      // fprintf(stderr, "destroy_leaf2 %d, cpos = %d / %d\n", b, cpos, ip->size());
      leaf_type* np = get_leaf(ip->child(cpos + 1));
      if (np->is_full()) {
        p->append(ip->data(cpos));
        ip->set_data(cpos, np->remove_median(lt_));
        np->move_data_lt_to(ip->data(cpos), p, lt_);
      } else {
        np->append(ip->data(cpos));
        ip->erase_pos(cpos);
        nmgr.destroy_leaf(b);
        try_destroy_internal_node(parent);
      }
    } else {
      assert(0);
    }
  }

  void erase(iterator it) {
    // assert(check());
    assert(it.m_ == this);
    assert(it.b_ != EMPTY_NODE);
    // fprintf(stderr, "erase it = %d, %d, val = %d\n", it.b_, it.pos_, it->value);
    if (is_internal(it.b_)) {
      iterator next_largest = find_bucket(*it, false);
      assert(next_largest.b_ != EMPTY_NODE);
      // fprintf(stderr, "next largest = %d, %d %d\n", next_largest->value, next_largest.b_, next_largest.pos_);
      assert(!is_internal(next_largest.b_));
      assert(lte(*next_largest, *it, lt_));
      get_internal(it.b_)->set_data(it.pos_, *next_largest);
      erase(next_largest);
    } else {
      // fprintf(stderr, "erasing a leaf %d at %d\n", it.b_, it.pos_);
      get_leaf(it.b_)->erase_pos(it.pos_);
      if (!get_leaf(it.b_)->size()) {
        // fprintf(stderr, "blank leaf %d at %d\n", it.b_, it.pos_);
        destroy_leaf_node(it.b_);
      }
    }
    // assert(check());
  }

  size_type erase(const value_type& val) {
    iterator it = find_bucket(val, true);
    if (it == end()) return 0;
    if (eq(*it, val, lt_)) {
      erase(it);
      return 1;
    }
    return 0;
  }

  void erase(iterator first, iterator last);

  void swap(multiset& x);
  void clear();

  key_compare key_comp() const;
  value_compare value_comp() const;

  iterator find(const value_type& val) const;
  size_type count(const value_type& val) const;

  void set_parent(int child, int parent) {
    if (is_internal(child)) {
      get_internal(child)->set_parent(parent);
    } else {
      get_leaf(child)->set_parent(parent);
    }
  }

  int new_internal_node(int parent, int left_child) {
    internal_type* np = nmgr.create_internal();
    int new_ib = nmgr.index_of(np);
    np->init(parent, left_child);
    np->set_next(EMPTY_NODE);
    set_parent(left_child, new_ib);
    return new_ib;
  }

  int new_leaf_node(int parent) {
    leaf_type* p = nmgr.create_leaf();
    p->init(parent);
    p->set_next(EMPTY_NODE);
    return nmgr.index_of(p);
  }

  void update_child_parents(int b, int parent) {
    internal_type* p = get_internal(b);
    for (int i = 0; i <= p->size(); i++) {
      set_parent(p->child(i), parent);
    }
  }

  int internal_split(int ib) {
    internal_type* p = get_internal(ib);
    int new_ib = new_internal_node(p->parent(), p->mid_child());
    internal_type* np = get_internal(new_ib);
    get_internal(ib)->move_half_to(np);
    update_child_parents(new_ib, new_ib);
    return new_ib;
  }

  void internal_insert(int internalb, T const &value, int newb, int left = 0) {
    // fprintf(stderr, "gi %d\n", internalb);
    get_internal(internalb)->insert(value, newb, left, lt_);
    set_parent(newb, internalb);
    // assert(check());
  }

  void leaf_split_one_chain(int b, std::function<void(T, int)> cb) {
    assert(!is_internal(b));
    int new_b = get_leaf(b)->detach_and_get_next();
    assert(get_leaf(new_b)->next() == EMPTY_NODE);
    // assert(get_leaf(new_b)->size() == leaf_type::capacity);

    T D[leaf_type::capacity * 2];
    int N = get_leaf(b)->copy_data_to(D);
    N += get_leaf(new_b)->copy_data_to(D + N);
    assert(N <= leaf_type::capacity * 2);
    int H = N / 2;
    std::nth_element(D, D + H, D + N, lt_);
    get_leaf(b)->copy_data_from(D, H);
    get_leaf(new_b)->copy_data_from(D + H + 1, N - H - 1);
    cb(D[H], new_b);
  }

  // Reservoir sampling (http://en.wikipedia.org/wiki/Reservoir_sampling).
  T pick_random_pivot(int b) {
    leaf_type* p = get_leaf(b);
    assert(leaf_type::capacity == p->size());

    // TODO: use randomized seed.
    Random rng(140384);

    // Replace R with the next buckets in the chain using reservoir sampling.
    for (int i = 1, nb = p->next(); nb != EMPTY_NODE; i++) {
      int j = rng.nextInt(i);
      leaf_type* np = get_leaf(nb);
      assert(np->size() > 0);
      if (j < p->size()) {
        assert(np->size() > 0);
        int k = rng.nextInt(np->size());
        np->swap_random_data_at_with(k, p, j);
      }
      nb = np->next();
    }
    return p->remove_median(lt_);
  }

  void add_chain(int head, int next) {
    assert(!is_internal(head));
    leaf_type* p = get_leaf(head);
    leaf_type* np = get_leaf(next);
    np->set_next(p->next());
    p->set_next(next);
  }

  void distribute_values(int b, T pivot, int chain[2]) {
    while (get_leaf(b)->size()) {
      T value = get_leaf(b)->remove_last();
      int i = gte(value, pivot, lt_);
      insert(chain[i], value, [&](int next) { assert(0); /* should not change head. */ });
    }
  }

  // TODO: optimize on near sorted input.
  void leaf_split_long_chain(int b, std::function<void(T, int)> cb) {
    // fprintf(stderr, "split long chain = %d\n", get_leaf(b)->size());
    T pivot = pick_random_pivot(b);
    // fprintf(stderr, "pivot = %d\n", pivot.value);
    int nb = new_leaf_node(get_leaf(b)->parent());
    leaf_type* p = get_leaf(b);
    p->move_data_at_least(pivot, get_leaf(nb), lt_);
    int chain[2] { b, nb };

    // fprintf(stderr, "fusion %d\n", b);
    int Nb = p->detach_and_get_next();
    int Lb = EMPTY_NODE, Rb = EMPTY_NODE;
    // TODO: optimize locality.
    int hi[leaf_type::capacity], nhi = 0;
    int lo[leaf_type::capacity], nlo = 0;
    while (true) {
      // fprintf(stderr, "%d %d\n", Lb, Rb);
      if (nhi && nlo) {
        assert(Lb != EMPTY_NODE && Rb != EMPTY_NODE);
        get_leaf(Lb)->fusion(get_leaf(Rb), hi, lo, nhi, nlo);
        if (!nhi) { add_chain(chain[0], Lb); Lb = EMPTY_NODE; }
        if (!nlo) { add_chain(chain[1], Rb); Rb = EMPTY_NODE; }
      } else if (Lb == EMPTY_NODE) {
        if (Nb == EMPTY_NODE) break;
        Lb = Nb;
        // fprintf(stderr, "%d size = %d / %d\n", Lb, get_leaf(Lb)->size(), leaf_type::capacity);
        // assert(get_leaf(Lb)->is_full());
        Nb = get_leaf(Nb)->detach_and_get_next();
      } else if (!nhi) {
        assert(Lb != EMPTY_NODE);
        get_leaf(Lb)->mark_hi(pivot, hi, nhi, lt_);
        if (!nhi){ add_chain(chain[0], Lb); Lb = EMPTY_NODE; }
      } else if (Rb == EMPTY_NODE) {
        if (Nb == EMPTY_NODE) break;
        Rb = Nb;
        // assert(get_leaf(Rb)->is_full());
        Nb = get_leaf(Nb)->detach_and_get_next();
      } else if (!nlo) {
        assert(Rb != EMPTY_NODE);
        get_leaf(Rb)->mark_lo(pivot, lo, nlo, lt_);
        if (!nlo){ add_chain(chain[1], Rb); Rb = EMPTY_NODE; }
      } else {
        assert(0);
      }
    }
    assert(Nb == EMPTY_NODE);

    // fprintf(stderr, "fusioned\n");
    if (Lb != EMPTY_NODE) distribute_values(Lb, pivot, chain), nmgr.destroy_leaf(Lb);
    if (Rb != EMPTY_NODE) distribute_values(Rb, pivot, chain), nmgr.destroy_leaf(Rb);
    // fprintf(stderr, "done split long chain nb = %d, %d\n", nb, chain[1]);
    assert(b == chain[0]);
    assert(nb == chain[1]);
    // assert(leaf_check());
    cb(pivot, nb);
  }

  // TODO: pararellize.
  void leaf_split(int b, std::function<void(T, int)> cb) {
    assert(!is_internal(b));
    assert(get_leaf(b)->next() != EMPTY_NODE);
    // fprintf(stderr, "leaf_split\n");

    if (get_leaf(get_leaf(b)->next())->next() == EMPTY_NODE) {
      leaf_split_one_chain(b, cb);
    } else {
      leaf_split_long_chain(b, cb);
    }
    // fprintf(stderr, "leaf_split2\n");
  }

  void compact_head_chain(int b) {
    leaf_type* p = get_leaf(b);
    while (!p->is_full() && p->next() != EMPTY_NODE) {
      int next = p->next();
      leaf_type* np = get_leaf(next);
      np->move_data_to(p);
      if (!np->size()) {
        p->set_next(np->next());
        nmgr.destroy_leaf(next);
      }
    }
  }

  bool split_chain(int b) {
    // assert(check());
    assert(!is_internal(b));

    compact_head_chain(b);

    if (get_leaf(b)->next() == EMPTY_NODE) return false;

    // fprintf(stderr, "\nsplit_chain %d\n", b);
    // assert(check());
    leaf_split(b, [&](T promotedValue, int nb) {
      assert(!is_internal(b));
      assert(!is_internal(nb));
      int parent = get_leaf(b)->parent();
      // fprintf(stderr, "leaf split parent = %d, is_internal = %d\n", parent, is_internal(parent));
      while (parent != EMPTY_NODE && nb != EMPTY_NODE) {
        if (get_internal(parent)->is_full()) {
          // fprintf(stderr, "parful\n");
          // Optional optimization: transfer_one_to_left_or_right();
          int inb = internal_split(parent);
          T promotedValueInternal = get_internal(parent)->remove_last();
          if (lt(promotedValue, promotedValueInternal, lt_)) {
            internal_insert(parent, promotedValue, nb);
          } else {
            internal_insert(inb, promotedValue, nb);
          }
          promotedValue = promotedValueInternal;
          nb = inb;
          parent = get_internal(nb)->parent();
        } else {
          // fprintf(stderr, "internal\n");
          internal_insert(parent, promotedValue, nb);
          nb = EMPTY_NODE;
        }
      }
      // assert(check());
      // fprintf(stderr, "nb = %d\n", nb);
      if (nb != EMPTY_NODE) {
        // fprintf(stderr, "OLD ROOT_ %d\n", root_);
        assert(parent == EMPTY_NODE);
        if (is_internal(root_)) {
          assert(get_internal(root_)->parent() == EMPTY_NODE);
        } else {
          assert(get_leaf(root_)->parent() == EMPTY_NODE);
        }
        root_ = new_internal_node(EMPTY_NODE, root_);
        internal_insert(root_, promotedValue, nb);
        // fprintf(stderr, "NEW ROOT_ %d\n", root_);
      }
    });
    // fprintf(stderr, "done split_chain %d\n", b);
    // assert(check());
    // fprintf(stderr, "promotedValue = %d\n", promotedValue);

    // fprintf(stderr, "\n\n\ndone split %d\n", b);
    // debug();
    // assert(check());
    return true;
  }

  void flush_pending_inserts(int ib) {
    // TODO: optimize
    while (get_internal(ib)->next() != EMPTY_NODE) {
      int next = get_internal(ib)->next();
      if (!get_leaf(next)->size()) {
        get_internal(ib)->set_next(get_leaf(next)->next());
        nmgr.destroy_leaf(next);
      } else {
        T value = get_leaf(next)->remove_last();
        // fprintf(stderr, "flush %d != %d, size = %d\n", next, ib, get_leaf(next)->size());
        int pos = get_internal(ib)->lower_pos(value, lt_);
        int c = get_internal(ib)->child(pos);
        // fprintf(stderr, "pos = %d / %d, child = %d\n", pos, get_internal(ib)->size(), c);
        if (is_internal(c)) {
          // fprintf(stderr, "internal next = %d\n", get_internal(c)->next());
          insert(get_internal(c)->next(), value, [&](int head) { get_internal(c)->set_next(head); });
        } else {
          // fprintf(stderr, "is leaf %d\n", c);
          insert(c, value, [&](int head) { assert(0); /* head shouldn't change */ });
        }
      }
    }
  }

  // Returns an iterator that points to th elower_bound of {value}.
  iterator find_bucket(value_type const &value, bool stop_at_internal) {
    // fprintf(stderr, "find_bucket %d\n", b);
    for (int b = root_; ; ) {
      // fprintf(stderr, "find bucket loop %d, stop_at_internal = %d\n", b, stop_at_internal);
      assert(b != EMPTY_NODE);
      if (is_internal(b)) {
        internal_type* p = get_internal(b);
        // fprintf(stderr, "find_bucket is_internal %d, size = %d, next = %d\n", b, p->size(), p->next());
        if (p->next() == EMPTY_NODE) {
          int pos = p->lower_pos(value, lt_);
          if (stop_at_internal && p->equal(pos, value, lt_)) {
            // fprintf(stderr, "find_bucket_internal %d\n", b);
            return iterator(this, b, pos); // Found in the internal bucket.
          }
          b = p->child(pos);    // Search the child.
          // fprintf(stderr, "find_bucket child %d, pos = %d / %d\n", b, pos, p->size());
        } else {
          flush_pending_inserts(b);
        }
      } else if (split_chain(b)) {
        // fprintf(stderr, "splited_chain %d\n", b);
        b = get_leaf(b)->parent();
      } else {
        leaf_type* p = get_leaf(b);
        int pos = p->lower_pos(value, lt_);
        // fprintf(stderr, "other %d, %d / %d, root_ par = %d\n", b, pos, p->size(), get_internal(root_)->parent());
        if (pos < p->size()) {
          // fprintf(stderr, "find_bucket_leaf %d\n", b);
          return iterator(this, b, pos);
        }
        if (!stop_at_internal) {
          // fprintf(stderr, "find_bucket_!stop_at_internal %d\n", b);
          assert(pos > 0);
          return iterator(this, b, pos - 1);
        }
        // Find the first next element after {b, pos}.
        for (b = p->parent(); b != EMPTY_NODE; ) {
          // fprintf(stderr, "b = %d, root_ = %d\n", b, root_);
          internal_type* ip = get_internal(b);
          pos = ip->lower_pos(value, lt_);
          // fprintf(stderr, "pos = %d / %d\n", pos, ip->size());
          if (pos < ip->size()) break;
          b = ip->parent();
        }
        // fprintf(stderr, "find_bucket_internal2 %d\n", b);
        return iterator(this, b, pos);
      }
    }
  }

  iterator lower_bound(const value_type& value) {
    // TODO: compact leaf
    // assert(check());
    // fprintf(stderr, "lower_bound %d\n", value);
    return find_bucket(value, true);
  }

  iterator upper_bound(const value_type& val) const;
  pair<iterator,iterator> equal_range(const value_type& val) const;

  bool check() {
    return rec(root_, [&](int b) {
      if (is_internal(b)) {
        internal_type* p = get_internal(b);
        for (int i = 0; i <= p->size(); i++) {
          int nb = p->child(i);
          if (is_internal(nb)) {
            assert(b == get_internal(nb)->parent());
          } else {
            assert(b == get_leaf(nb)->parent());
          }
        }
      }
      return 1;
    });
  }
};





/*
  void leaf_sort() {
    assert(is_valid());
    if (P) {
      assert(!locked);
      sort(D, D + N);
      P = 0;
    }
  }

  T leaf_promote_first() {
    assert(is_valid());
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

  T internal_promote_first(int *C) {
    assert(is_valid());
    T ret = D[0];
    N--;
    for (int i = 0; i < N; i++) {
      D[i] = D[i + 1];
      C[i] = C[i + 1];
    }
    C[N] = C[N + 1];
    return ret;
  }
};

template <typename T>
class CTree {

  bool optimize(int b = 0) {
    if (b == 0) {
      if (root_ == 0) return false;
      while (optimize(root_)) {
        if (!is_leaf(root_) && INTERNAL_BUCKET(root_)->size() == 0) {
          // fprintf(stderr, "PROMOTE ROOT_\n");
          int c = INTERNAL_BUCKET(root_)->child(0);
          delete_leaf_bucket(root_);
          root_ = c;
          set_parent(root_, 0);
        }
        // fprintf(stderr, "optimize\n");
      }
      return false;
    }

    if (is_leaf(b)) {
      bool splitted = split_chain(b);
      while (split_chain(b));
      assert(is_leaf(b));
      LEAF_BUCKET(b)->leaf_sort();
      return splitted;
    }

    bool changed = false;
    for (int i = 0; i <= INTERNAL_BUCKET(b)->size(); i++) {
      if (optimize(INTERNAL_BUCKET(b)->child(i))) {
        i = -1;
        changed = true;
      }
    }

    // OPTIONAL:
    // for (int i = 0; i < BUCKET(b)->size(); i++) {
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

  bool check(int b = 0, T lo = -2147483648) {
    if (b == 0) b = root_;
    // fprintf(stderr, "check %d, leaf = %d\n", b, BUCKET(b)->is_leaf());
    if (is_leaf(b)) {
      if (b == root_) assert(LEAF_BUCKET(b)->parent() == 0);
      return LEAF_BUCKET(b)->leaf_check(lo, true, 0, false);
    }
    IBucket<T> *ib = INTERNAL_BUCKET(b);
    if (parent_of(ib->child(0)) != b) {
      fprintf(stderr, "parent mismatch %d != %d\n", parent_of(ib->child(0)), b);
      return false;
    }
    if (ib->size() && !check(ib->child(0), lo, ib->data(0))) return false;
    for (int i = 0; i < ib->size(); i++) {
      assert(i == 0 || ib->data(i - 1) <= ib->data(i));
      if (parent_of(ib->child(i + 1)) != b) {
        fprintf(stderr, "parent mismatch internal %d != %d\n", parent_of(ib->child(i + 1)), b);
        return false;
      }
      if (!check(ib->child(i + 1), ib->data(i), (i + 1 < ib->size()) ? ib->data(i + 1) : 2147483647)) return false;
    }
    return true;
  }
};


class CTree {

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
      ib->set_data(rpos, M->remove_last());
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
      ib->set_data(rpos, M->remove_last());
    }
    assert(!R->is_full());
    R->internal_insert(ib->data(rpos), M->child(M->size()), -1);
    ib->internal_erase(rpos, 0);
    delete M;
    return true;
  }

  int leaf_shift_left(int b, int pos) {
    int L = child(b, pos);
    int R = child(b, pos + 1);
    assert(BUCKET(L)->is_leaf());
    assert(BUCKET(R)->is_leaf());
    if (BUCKET(L)->next != 0) return false;
    if (BUCKET(R)->next != 0) return false;

    // Move from M to L as many as possible.
    int changed = 0;
    while (!BUCKET(L)->is_full() && BUCKET(R)->size()) {
      leaf_insert(L, BUCKET(b)->D[pos]);
      BUCKET(b)->D[pos] = BUCKET(R)->leaf_promote_first();
      changed = 1;
    }
    if (!BUCKET(L)->is_full() && !BUCKET(R)->size()) {
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
    assert(BUCKET(L)->next == 0);
    assert(BUCKET(R)->next == 0);

    // Move from R to L as many as possible.
    bool changed = false;
    while (!BUCKET(L)->is_full() && BUCKET(R)->size() && numMove-- > 0) {
      internal_insert(L, BUCKET(b)->D[pos], child(R, 0));
      BUCKET(b)->D[pos] = BUCKET(R)->internal_promote_first(CHILDREN(R));
      changed = true;
    }
    // if (!BUCKET(L)->is_full() && !BUCKET(R)->size()) {
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
    assert(BUCKET(L)->next == 0);
    assert(BUCKET(R)->next == 0);

    // Move from L to R as many as possible.
    bool changed = false;
    while (!BUCKET(R)->is_full() && BUCKET(L)->size() && numMove-- > 0) {
      internal_insert(R, BUCKET(b)->D[pos], child(L, BUCKET(L)->size()), -1);
      BUCKET(b)->D[pos] = BUCKET(L)->remove_last();
      changed = true;
    }
    // if (!BUCKET(R)->is_full());
    // R->internal_insert(ib->data(rpos), L->child(L->size()), -1);
    // ib->internal_erase(rpos, 0);
    // delete L;
    return changed;
  }


  int internal_find_child_pos(int b, int c) {
    int *C = CHILDREN(b);
    for (int i = 0; i <= BUCKET(b)->size(); i++) {
      if (C[i] == c) return i;
    }
    assert(0);
    return 0;
  }


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
    if (B->next == 0) slack += B->slack(), start = 0;
    for (int i = 1; i <= BUCKET(b)->size(); i++) {
      B = BUCKET(CHILDREN(b)[i]);
      if (B->next != 0) {
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

  transfer_one_to_left_or_right() {
    assert(!BUCKET(parent)->is_leaf());
    int pp = BUCKET(parent)->parent;
    if (pp != 0) {
      assert(!BUCKET(pp)->is_leaf());
      int pos = internal_find_child_pos(pp, parent);
      // assert(check());

      if (pos > 0 && BUCKET(parent)->D[0] < promotedValue && internal_shift_left(pp, pos - 1, 1)) {
        // fprintf(stderr, "shift left %d, p = %d, b = %d, root_ = %d, pos = %d / %d\n", pp, parent, b, root_, pos, BUCKET(pp)->size());
        assert(!BUCKET(parent)->is_full());
        // assert(check());
        internal_insert(parent, promotedValue, nb);
        nb = 0;
        // assert(check());
        break;
      } else if (pos < BUCKET(pp)->size() && promotedValue < BUCKET(parent)->D[BUCKET(parent)->size() - 1] && internal_shift_right(pp, pos, 1)) {
        // assert(check());
        assert(!BUCKET(parent)->is_full());
        internal_insert(parent, promotedValue, nb);
        nb = 0;
        // assert(check());
        break;
      }
    }
  }
*/
}

#endif
