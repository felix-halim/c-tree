/*
TODO:
- test fusion performance
- test with payload
- measure index size
- streaming sort, like merge sort! may be good for disk based.
- more compact than btree
- explain evolution from cracking
*/

var CELL_W = 25;
var CELL_H = 20;
var LEAF_CAPACITY = 11;
var NEXT_BUCKET_ID = 1;

function assert(pred) { if (!pred) throw new Error(); }

function NewBucket() {
  var id = NEXT_BUCKET_ID++;
  var data = [];
  var index = []; // Mini indexes inside this bucket.
  var next = null;

  function append(num) {
    assert(typeof num == 'number');
    var ret = this;
    if (data.length < LEAF_CAPACITY) {
      data.push(num);
    } else {
      ret = new NewBucket();
      ret.append(num);
      ret.set_next(this);
    }
    return ret;
  }

  function traverse(cb) {
    var b = this;
    for (var i = 0; b; i++) {
      cb('leaf_bucket', i, b.id(), b.data(), b.index(), b.next());
      b = b.next();
    }
  }

  function n_chains() {
    var ret = 0, b = this;
    while (b) {
      ret++;
      b = b.next();
    }
    return ret;
  }

  function add_index(i) {
    for (var j = 0; j < index.length; j++) {
      if (index[j].pos > i) {
        index.splice(j, 0, { pos: i });
        break;
      }
    }
  }

  function swap(i, j) {
    var t = data[i];
    data[i] = data[j];
    data[j] = t;
  }

  function sort(L, R) {
    for (var i = L; i < R; i++)
      for (var j = i + 1; j < R; j++)
        if (data[i] > data[j]) swap(i, j);
  }

  function crack(L, R) {
    var i = (L + R) >> 1; // Naive.
    swap(i, L);
    var pivot = data[L];
    for (var x = L + 1, y = R - 1; x <= y; ) {
      while (x <= y && data[x] < pivot) x++;
      while (x <= y && data[y] > pivot) y--;
      if (x < y) swap(x, y);
    }
    swap(L, x - 1);
    return x - 1;
  }

  function flush_pending_inserts() {
    assert(index.length);
    var j = index.length - 1;
    if (j == 0) {
      index[j].pos = data.length;
      index[j].sorted = false;
    } else {
      while (index[j].pos < data.length) {
        var k = index[j].pos++;
        for (var i = j - 1; i >= 0; i--) {
          index[i].sorted = false;
          if (data[k] >= data[index[i].pos]) break;
          swap(k, index[i].pos + 1);
          swap(index[i].pos, index[i].pos + 1);
          k = index[i].pos++;
        }
      }
    }
  }

  function shuffle_out(p) {
    assert(index.length);
    for (var i = 0, j = 0; i < data.length; i++) {
      if (i >= index[j].pos) j++;
      assert(j < index.length);
      if (data[i] == p.value) {
        // Shuffle it to the next index position.
        if (index[j].pos == data.length) {
          swap(i, data.length - 1);
          data.pop();
          index[j].sorted = false;
          index[j].pos--;
          if (index[j].pos == 0 || (j > 0 && index[j].pos == index[j - 1].pos)) index.splice(j, 1);
          p.type = 'fusion';
          p.message = { text: p.value + ' is shuffled out', color: 'green' };
          return p.cb(p, function () {
            p.type = 'delete_end';
            p.cb(p);
            // shuffle_out(p);
          });
        } else {
          if (j > 0 && index[j - 1].pos == i) {
            swap(i, i + 1);
            i++;
          }
          swap(i, index[j].pos - 1);
          swap(index[j].pos - 1, index[j].pos);
          index[j].pos--;
          index[j].sorted = false;
          if (index[j].pos == 0 || (j > 0 && index[j].pos == index[j - 1].pos)) index.splice(j, 1);
          p.type = 'fusion';
          p.message = { text: p.value + ' is shuffled to pos = ' + (index[j].pos + 1), color: 'green' };
          return p.cb(p, function () {
            shuffle_out(p);
          });
        }
      }
    }
    assert(false);
  }

  function lower_bound(p) {
    if (!index.length) index.push({ pos: data.length });

    if (index[index.length - 1].pos < data.length) {
      flush_pending_inserts();
      return p.cb(p, function () {
        lower_bound(p);
      });
    }

    var L = 0, R = data.length;
    // console.log(JSON.stringify(data) + ' ' + JSON.stringify(index));
    assert(index[index.length - 1].pos == data.length);
    for (var i = 0; i < index.length; i++) {
      R = index[i].pos;
      if (p.value < data[R]) break;
      if (i + 1 < index.length) L = R; else break;
    }
    if (R - L > 3) {
      var i = 0;
      while (i == 0) i = crack(L + (L == 0 ? 0 : 1), R);
      p.message = { text: 'Crack the bucket range = [' + L + ', ' + R + '), index = ' + i, color: 'red' };
      add_index(i);
      if (p.value < data[i]) {
        R = i;
      } else {
        L = i;
      }

      p.cb(p, function () {
        lower_bound(p);
      });
    } else {
      if (!index[i].sorted) {
        index[i].sorted = true;
        sort(L, R);
      }
      var found = false;
      for (; L < R; L++) if (data[L] == p.value) found = true;
      if (found) {
        p.message = { text: 'Found value ' + p.value, color: 'green' };
      } else {
        p.message = { text: 'Value ' + p.value + ' not found', color: 'red' };
      }
      p.cb(p, function () {
        p.type = 'end';
        p.result = found;
        p.cb(p, function (delete_it) {
          if (!delete_it) return;
          // console.log(L + ' ' + R + ' ' + JSON.stringify(data) + ' ' + JSON.stringify(index) + ' v = ' + p.value);
          assert(found);
          shuffle_out(p);
        });
      });
    }
  }

  function clear_index() {
    index = [];
  }

  return {
    append: append,
    traverse: traverse,
    id: function () { return id; },
    size: function () { return data.length; },
    data: function () { return data; },
    index: function () { return index; },
    next: function () { return next; },
    set_next: function (n) { next = n; clear_index(); if (n) n.clear_index(); },
    remove: function (i) { return data.splice(i, 1)[0]; },
    n_chains: n_chains,
    lower_bound: lower_bound,
    clear_index: clear_index,
  };
}


function NewTrimmer() {
  var root_index = []; // Can be implemented using BTree, ART, etc.
  var root_pointer = [ NewBucket() ];
  var inserted_values = {};

  function root_index_position(value) {
    // BTree or ART is much faster that plain array here.
    // This is just for illustration of root index maintains a sorted values (and pointers).
    var i = 0;
    for (; i < root_index.length && root_index[i] < value; i++);
    return i;
  }

  function insert(value) {
    // This demo only supports unique values. Special care must be taken to support unique values.
    if (inserted_values[value]) return false;
    inserted_values[value] = true;

    var i = root_index_position(value);
    root_pointer[i] = root_pointer[i].append(value);

    return true;
  }

  function traverse(cb) {
    cb('root_index', root_index, root_pointer);
  }

  function pick_median(ith, leafb) {
    var arr = [];
    for (var b = leafb; b; b = b.next())
      for (var i = 0; i < b.size(); i++)
        arr.push(b.data()[i]);
    arr.sort(function (a, b) { return a - b; });
    var pivot = arr[arr.length >> 1];
    for (var b = leafb; b; b = b.next())
      for (var i = 0; i < b.size(); i++)
        if (b.data()[i] == pivot) {
          var t = b.remove(i);
          assert(!((t < pivot) || (t > pivot)));
          if (!b.size()) {
            assert(b == leafb);
            root_pointer[ith] = b.next();
          }
          return pivot;
        }
    assert(false);
  }

  function fusion_split(b1, b2, pivot) {
    var x = 0, y = 0;
    // Faster fusion is in the paper. This is for illustration only.
    while (true) {
      while (x < b1.size() && b1.data()[x] < pivot) x++;
      while (y < b2.size() && b2.data()[y] > pivot) y++;
      if (x == b1.size() || y == b2.size()) break;
      var t = b1.data()[x];
      b1.data()[x] = b2.data()[y];
      b2.data()[y] = t;
    }
    return { i1 : x, i2 : y };
  }

  function attach(b, i) {
    b.clear_index();
    if (!root_pointer[i]) {
      b.set_next(null);
      root_pointer[i] = b;
    } else if (root_pointer[i].size() == LEAF_CAPACITY) {
      b.set_next(root_pointer[i]);
      root_pointer[i] = b;
    } else if (b.size() == LEAF_CAPACITY) {
      b.set_next(root_pointer[i].next());
      root_pointer[i].set_next(b);
    } else {
      while (b.size() > 0 && root_pointer[i].size() < LEAF_CAPACITY) {
        root_pointer[i].append(b.remove(0));
      }
      if (b.size()) {
        b.set_next(root_pointer[i]);
        root_pointer[i] = b;
      }
    }
  }

  function split_chain(p) {
    var b = p.chain;
    if (!b) {
      p.root_i += (p.value < p.pivot) ? 0 : 1;
      return fusion(p);
    }

    var b1 = b, b2 = b.next();
    if (p.marked) {
      var pos = fusion_split(b1, b2, p.pivot);
      if (pos.i1 == b1.size()) {
        attach(b1, p.root_i);
        p.chain = b2;
        p.message = { text: 'Red bucket move to the left chain', color: 'red' };
        p.cb(p, function () {
          if (pos.i2 == b2.size()) {
            var b = b2.next();
            attach(b2, p.root_i + 1);
            p.chain = b;
            p.message = { text: 'Blue bucket move to the right chain', color: 'blue' };
            p.cb(p, function () {
              p.marked = false;
              split_chain(p);
            });
          } else {
            p.chain = b2;
            p.marked = false;
            split_chain(p);
          }
        });
      } else if (pos.i2 == b2.size()) {
        b1.set_next(b2.next());
        attach(b2, p.root_i + 1);
        p.chain = b1;
        p.message = { text: 'Blue bucket move to the right chain', color: 'blue' };
        p.cb(p, function () {
          p.chain = b1;
          p.marked = false;
          split_chain(p);
        });
      }
    } else if (b2) {
      fusion_split(b1, b2, p.pivot);
      p.chain = b1;
      p.message = { text: 'Fusion technique to separate red and blue elements', color: 'orange' };
      p.cb(p, function () {
        p.marked = true;
        split_chain(p);
      });
    } else if (b1) {
      b1.set_next(b2 = NewBucket());
      for (var j = 0; j < b1.size(); j++) {
        if (b1.data()[j] > p.pivot) {
          b2.append(b1.remove(j));
          j--;
        }
      }
      p.chain = b1;
      p.message = { text: 'Partition the red and blue elements to its own bucket', color: 'brown' };
      p.cb(p, function () {
        p.marked = true;
        split_chain(p);
      });
    } else {
      assert(false);
    }
  }

  function fusion(p) {
    if (root_pointer[p.root_i].next()) {
      // Perfect median finding, in practice, rough median finding is faster.
      p.chain = root_pointer[p.root_i];
      p.pivot = pick_median(p.root_i, p.chain);
      root_index.splice(p.root_i, 0, p.pivot);
      root_pointer.splice(p.root_i + 1, 0, null);
      root_pointer[p.root_i] = null;
      p.message = { text: 'Performing stochastic split_chain with pivot = ' + p.pivot, color: 'blue' };
      p.cb(p, function () {
        p.marked = false;
        split_chain(p);
      });
    } else {
      root_pointer[p.root_i].lower_bound(p);
    }
  }

  function erase(v, cb) {
    var result = false;
    lower_bound(v, function (action, resume) {
      if (action.type == 'end') {
        if (!action.result) {
          // Check the root.
          var i = root_index_position(v);
          if (root_index[i] == v) {
            result = true;
            return cb(action, function () {
              delete_root(i);
              cb({ type: 'fusion'}, function () {
                cb({ type: 'delete_end'});
              });
            });
          }
          return cb({ type: 'delete_end'});
        }
        result = true;
      }
      cb(action, resume);
    });
    if (!inserted_values[v]) inserted_values[v] = false;
    assert(inserted_values[v] == result);
    inserted_values[v] = false;
  }

  function delete_root(i) {
    // There are several options:
    // - Attach to the right chain to the left (this demo).
    // - Break chains and get the biggest element to replace.
    var x = root_index.splice(i, 1)[0];
    // console.log('deleted root = ' + x);
    var b = root_pointer.splice(i + 1, 1)[0];
    var last = root_pointer[i];
    while (last.next()) last = last.next();
    last.set_next(b);
    last.clear_index();
    b.clear_index();
  }

  function lower_bound(value, cb) {
    var i = root_index_position(value);
    if (i < root_index.length && root_index[i] == value) {
      // console.log('delete root ' + value);
      return cb({ type: 'end', result: true }, function (delete_it) {
        if (!delete_it) return;
        delete_root(i);
        cb({ type: 'fusion' }, function () {
          cb({ type: 'delete_end' });
        });
      });
    }
    //console.log('delete normal ' + value);
    fusion({ type: 'fusion', value: value, root_i: i, cb: cb });
  }

  function has_value(v) { return inserted_values[v]; }

  return {
    traverse: traverse,
    insert: insert,
    lower_bound: lower_bound,
    erase: erase,
    has_value: has_value,
  };
}
