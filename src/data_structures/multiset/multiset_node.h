#ifndef CRACK_NODE_H_
#define CRACK_NODE_H_

#include <vector>

namespace trimmer {

using std::vector;

template<typename T, typename compare>
class node {
 protected:
  int size_;          // Size of this node.
  int next_;          // Points to the head of leaf_node.
  int parent_;        // Parent node id.

  bool valid() const { return size_ >= 0; }

 public:

  node(): size_(0), parent_(0) {}

  ~node() {
    assert(valid());
    size_ = EMPTY_NODE;
  }

  int size() const {
    assert(valid());
    return size_;
  }

  int parent() const {
    assert(valid());
    return parent_;
  }

  int next() const {
    assert(valid());
    return next_;
  }

  void set_parent(int parent) {
    assert(valid());
    assert(parent >= -1);
    parent_ = parent;
  }

  void set_next(int next) {
    assert(valid());
    next_ = next;
  }
};

// B-Tree internal node.
template<typename T, typename compare>
class internal_node : public node<T, compare> {
 public:
  const static int capacity = 32;

  void init(int parent, int left_child_id) {
    this->parent_ = parent;
    this->size_ = 0;
    this->child_[0] = left_child_id;
  }

  bool is_full() const {
    assert(this->valid());
    assert(this->size_ <= capacity);
    return this->size_ == capacity;
  }

  T& data(int i) {
    assert(this->valid());
    assert(i >= 0 && i < this->size_);
    return this->data_[i];
  }

  int lower_pos(T const &value, compare &lt) const {
    assert(this->valid());
    return sorted_lower_pos(data_, this->size_, value, lt);
  }

  int mid_child() const {
    assert(this->valid());
    assert(is_full());
    return child_[this->size_ / 2];
  }

  void move_half_to(internal_node* that) {
    assert(this->valid());
    assert(is_full());
    assert(this->size_ % 2 == 0);
    int H = this->size_ / 2;
    for (int i = H, j = 0; i < this->size_; i++) {
      that->data_[j++] = data_[i];
      that->child_[j] = child_[i + 1];
    }
    this->size_ = that->size_ = H;
  }

  T remove_last() {
    assert(this->valid());
    assert(this->size_ > 0);
    return data_[--this->size_];
  }

  void insert(T const &value, int nb, int left, compare &lt) {
    assert(this->valid());
    assert(!is_full());
    // fprintf(stderr, "ii %d\n", b);
    int i = this->size_ - 1;
    while (i >= 0 && lt(value, data_[i])) {
      data_[i + 1] = data_[i];
      child_[i + 2] = child_[i + 1];
      i--;
    }
    data_[i + 1] = value;
    if (left == -1) {
      child_[i + 2] = child_[i + 1];
      child_[i + 1] = nb;
    } else {
      child_[i + 2] = nb;
    }
    this->size_++;
  }

  bool equal(int pos, T const &value, compare &lt) {
    assert(this->valid());
    return pos >= 0 && pos < this->size_ && eq(data_[pos], value, lt);
  }

  void set_data(int i, T const &value) {
    assert(this->valid());
    assert(i >= 0 && i < this->size_);
    this->data_[i] = value;
  }

  int child(int i) const {
    assert(this->valid());
    assert(i >= 0 && i <= this->size_);
    return this->child_[i];
  }

  int child_pos(int c) const {
    assert(this->valid());
    for (int i = 0; i <= this->size_; i++)
      if (this->child_[i] == c) return i;
    assert(0);
    return this->size_;
  }

  // void set_child(int i, int node_id) {
  //   assert(this->valid());
  //   assert(i >= 0 && i <= this->size_);
  //   this->child_[i] = node_id;
  // }

  void erase_pos(int pos, int stride = 0) {
    assert(this->valid());
    assert(pos >= 0 && pos < this->size_);
    this->size_--;
    while (pos < this->size_) {
      data_[pos] = data_[pos + 1];
      child_[pos + stride] = child_[pos + stride + 1];
      pos++;
    }
    if (!stride) child_[pos] = child_[pos + 1];
  }

  void move_to_lesser(internal_node* that, T value, compare &lt) {
    assert(this->valid());
    assert(that->size_ > 0);
    assert(lt(that->data_[that->size_ - 1], value));
    that->data_[that->size_] = value;
    that->child_[++that->size_] = child_[0];
    for (int i = 0; i < this->size_; ) {
      that->data_[that->size_] = data_[i];
      that->child_[++that->size_] = child_[++i];
    }
  }


 private:
  int child_[capacity + 1]; // Children node ids.
  T data_[capacity];        // Data values.
};


// B-Tree leaf node.
template<typename T, typename compare>
class leaf_node : public node<T, compare> {
 public:
  const static int capacity = 256;

  void init(int parent) {
    assert(this->valid());
    this->set_parent(parent);
  }

  const T& data(int i) const {
    assert(this->valid());
    assert(i >= 0 && i < this->size_);
    return this->data_[i];
  }

  bool is_full() const {
    assert(this->valid());
    assert(this->size_ <= capacity);
    return this->size_ == capacity;
  }

  T remove_last() {
    assert(this->valid());
    assert(this->size_ > 0);
    return data_[--this->size_];
  }

  void mark_hi(T P, int *hi, int &nhi, compare &lt) const {
    assert(this->valid());
    for (int i = 0; i < this->size_; i++) {
      hi[nhi] = i;
      nhi += gte(data_[i], P, lt);
    }
  }

  void mark_lo(T P, int *lo, int &nlo, compare &lt) const {
    assert(this->valid());
    for (int i = 0; i < this->size_; i++) {
      lo[nlo] = i;
      nlo += lt(data_[i], P);
    }
  }

  // Only swaps as necessary.
  void fusion(leaf_node *that, int *hi, int *lo, int &nhi, int &nlo) {
    assert(this->valid());
    T *Lp = data_, *Rp = that->data_;
    int m = std::min(nhi, nlo); assert(m > 0);
    int *hip = hi + nhi - 1, *lop = lo + nlo - 1;
    nhi -= m; nlo -= m;
    while (m--) std::swap(Lp[*(hip--)], Rp[*(lop--)]);
  }

  void move_data_to(leaf_node *that) {
    assert(this->valid());
    while (!that->is_full() && this->size_) {
      that->data_[that->size_++] = data_[--this->size_];
    }
    that->sorted_ = false;
  }

  void move_data_at_least(T const &value, leaf_node* that, compare &lt) {
    assert(this->valid());
    for (int i = 0; i < this->size_; i++) {
      if (gte(data_[i], value, lt)) {
        that->append(data_[i]);
        data_[i--] = data_[--this->size_];
      }
    }
  }

  void move_data_lt_to(T const &value, leaf_node* that, compare &lt) {
    assert(this->valid());
    for (int i = 0; i < this->size_; i++) {
      if (lt(data_[i], value)) {
        that->append(data_[i]);
        data_[i--] = data_[--this->size_];
      }
    }
  }

  void swap_random_data_at_with(int pos, leaf_node* that, int that_pos) {
    assert(this->valid());
    assert(this->size_ > 0 && pos < this->size_);
    assert(that->valid());
    assert(that->size_ > 0 && that_pos < that->size_);
    std::swap(data_[pos], that->data_[that_pos]);
  }

  const T& remove_median(compare &lt) {
    assert(this->valid());
    assert(this->size_ > 0);
    int H = this->size_ / 2;
    std::nth_element(data_, data_ + H, data_ + this->size_, lt);
    std::swap(data_[H], data_[--this->size_]);
    sorted_ = false;
    return data_[this->size_];
  }

  void copy_data_from(T *from, int cnt) {
    assert(this->valid());
    for (int i = 0; i < cnt; i++) {
      data_[i] = from[i];
    }
    this->size_ = cnt;
    sorted_ = false;
  }

  int copy_data_to(T *to) const {
    assert(this->valid());
    assert(this->size_ <= capacity);
    for (int i = 0; i < this->size_; i++) {
      to[i] = data_[i];
    }
    return this->size_;
  }

  int detach_and_get_next() {
    assert(this->valid());
    int ret = this->next_;
    this->next_ = EMPTY_NODE;
    sorted_ = false;
    return ret;
  }

  int lower_pos(T const &value, compare& lt) {
    assert(this->valid());
    assert(this->next_ == EMPTY_NODE);
    if (!sorted_) {
      std::sort(data_, data_ + this->size_, lt);
      sorted_ = true;
    }
    return sorted_lower_pos(data_, this->size_, value, lt);
  }

  // Append the value to this leaf node, possibly creating a linked chain of leaf nodes.
  void append(T const &value) {
    assert(this->valid());
    assert(this->size_ < capacity);
    data_[this->size_++] = value;
    sorted_ = false;
  }

  // Bulk append.
  void append(T const *values, int cnt) {
    assert(this->valid());
    assert(!this->size_);
    assert(capacity == cnt);
    memcpy(data_, values, sizeof(T) * capacity);
    this->size_ = capacity;
    sorted_ = false;
  }

  void erase_pos(int pos) {
    assert(this->valid());
    assert(pos >= 0 && pos < this->size_);
    data_[pos] = data_[--this->size_];
    sorted_ = false;
  }

 private:
  int sorted_;        // Whether the data_ array is in sorted state.
  T data_[capacity];  // Data values.
};

}

#endif
