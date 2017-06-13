/*
TODO:
- test fusion performance
- test with payload
- drop/forget indexes
- measure index size
- btree from eager only 10 secs
- streaming sort, like merge sort! may be good for disk based.
- cracker indexes as lightweight as btree 
- more compact than btree
- explain evolution from cracking
*/
var CELL_W = 25;
var CELL_H = 25;
var LEAF_CAPACITY = 4;
var INTERNAL_CAPACITY = 3;
var NEXT_BUCKET_ID = 1;

var LEAF_GAP = CELL_W * 2;
var LEVEL_GAP = CELL_H * 4;

function assert(pred) { if (!pred) throw new Error(); }

function NewBucket(p, first_child) {
  var id = NEXT_BUCKET_ID++;
  var data = [];
  var children = null;
  var next = null;
  var tail = null;
  var parent = p;
  var sorted = false;
  var pos_x = 0;
  var pos_y = 0;
  var thisRef = {
    append: append,
    is_full: is_full,
    is_leaf: is_leaf,
    equal: equal,
    internal_lower_pos: internal_lower_pos,
    move_data_at_least: move_data_at_least,
    detach_and_get_next: detach_and_get_next,
    leaf_promote_last: leaf_promote_last,
    mark_hi: mark_hi,
    mark_lo: mark_lo,
    fusion: fusion,
    internal_insert: internal_insert,
    leaf_lower_pos: leaf_lower_pos,
    min_width: min_width,
    move_half_to: move_half_to,
    print_chain: print_chain,
    snapshot: snapshot,
    rec_set_position: rec_set_position,
    x: function () { return pos_x; },
    y: function () { return pos_y; },
    id: function () { return id; },
    data: function () { return data; },
    tail: function () { return tail; },
    next: function () { return next; },
    parent: function () { return parent; },
    sorted: function () { return sorted; },
    child: function (i) { return children[i]; },
    child_size: function (i) { return children.length; },
    size: function (i) { return data.length; },
    set_tail: function (t) { tail = t; },
    set_next: function (n) { next = n; },
    set_parent: function (p) { parent = p; },
    set_data: function (i, v) { data[i] = v; },
    set_child: function (i, v) { children[i] = v; },
    set_data_length: function (len) { data.length = len; if (!is_leaf()) children.length = len + 1; },
    set_position: function (x, y) { pos_x = x; pos_y = y; },
    remove_last: function () { assert(data.length > 0); var t = data.pop(); thisRef.set_data_length(data.length); return t; },
    last_data_is_at_least : function (value) { return (data[data.length - 1] >= value) ? 1 : 0; },
    mid_child : function () { assert(!is_leaf()); return children[(data.length + 1) >> 1]; },
  };

  if (first_child) {
    children = [ first_child ];
    first_child.set_parent(thisRef);
  }

  function print_chain() {
    console.log('chain:');
    for (var b = thisRef; b; b = b.next()) {
      console.log(b.data());
    }
    console.log('chain done');
  }

  function is_leaf() { return children == null; }

  function is_full() { return data.length == (is_leaf() ? LEAF_CAPACITY : INTERNAL_CAPACITY); }

  function append(num) {
    data.push(num);
    assert(data.length <= LEAF_CAPACITY);
  };

  function move_half_to(that) {
    var H = (data.length + 1) >> 1;
    for (var i = H, j = 0; i < data.length; i++) {
      that.set_data(j++, data[i]);
      that.set_child(j, children[i + 1]);
    }
    that.set_data_length(data.length - H);
    thisRef.set_data_length(H);
  }

  function leaf_sort() {
    if (!sorted) {
      data.sort(function (a, b) { return a - b; });
      sorted = true;
    }
  }

  function leaf_promote_last() {
    assert(data.length > 0);
    return data.pop();
  }

  function leaf_lower_pos(value) {
    leaf_sort();
    return internal_lower_pos(value);
  }

  function leaf_promote_first() {
    sorted = false;
    var smallest_pos = 0;
    var pos = 1;
    while (pos < N) {
      if (data[pos] < data[smallest_pos]) smallest_pos = pos;
      pos++;
    }
    var ret = data[smallest_pos];
    var last = data.pop();
    if (smallest_pos < data.length) data[smallest_pos] = last;
    return ret;
  }

  function detach_and_get_next() {
    var ret = next;
    next = tail = null;
    return ret;
  }

  function mark_hi(P, hi, nhi) {
    for (var i = 0; i < data.length; i++) {
      hi[nhi] = i;
      nhi += data[i] >= P;
    }
    return nhi;
  }

  function mark_lo(P, lo, nlo) {
    for (var i = 0; i < data.length; i++) {
      lo[nlo] = i;
      nlo += data[i] < P;
    }
    return nlo;
  }

  // Only swaps as necessary.
  function fusion(that, hi, lo, nhi, nlo) {
    var Lp = thisRef.data(), Rp = that.data();
    var m = Math.min(nhi, nlo); assert(m > 0);
    while (m--) {
      var t = Lp[hi[--nhi]];
      Lp[hi[nhi]] = Rp[lo[--nlo]];
      Rp[lo[nlo]] = t;
    }
    return [nhi, nlo];
  }

  function leaf_erase_pos(pos) {
    assert(pos >= 0 && pos < data.length);
    var ret = data[pos];
    var last = data.pop();
    if (pos < data.length) data[pos] = last;
    sorted = false;
    return ret;
  }

  function leaf_largest_pos() {
    assert(!next);
    var largest_pos = 0, pos = 1;
    while (pos < data.length) {
      if (data[pos] >= data[largest_pos])
        largest_pos = pos;
      pos++;
    }
    return largest_pos;
  }

  function leaf_erase(v) {
    for (var i = 0; i < data.length; i++) {
      if (data[i] == v) {
        var last = data.pop();
        if (i < data.length) data[i] = last;
        sorted = false;
        return true;
      }
    }
    return false;
  }

  function internal_lower_pos(value) {
    var pos = 0;
    while (pos < data.length && data[pos] < value) pos++;
    return pos;
  }

  function move_data_at_least(value, that) {
    for (var i = 0; i < data.length; i++) {
      if (data[i] >= value) {
        that.append(data[i]);
        var last = data.pop();
        if (i < data.length) data[i--] = last;
      }
    }
  }

  function equal(pos, value) {
    return pos >= 0 && pos < data.length && data[pos] == value;
  }

  function internal_insert(value, nb, left) {
    assert(!is_full());
    var i = data.length - 1;
    while (i >= 0 && data[i] > value) {
      data[i + 1] = data[i];
      children[i + 2] = children[i + 1];
      i--;
    }
    data[i + 1] = value;
    if (left == -1) {
      children[i + 2] = children[i + 1];
      children[i + 1] = nb;
    } else {
      children[i + 2] = nb;
    }
  }

  function min_width() {
    if (is_leaf()) {
      var n_leaves = 1;
      for (var b = next; b; b = b.next()) n_leaves++;
      return n_leaves * LEAF_GAP + CELL_W;
    }
    var ret = 0;
    for (var i = 0; i < children.length; i++) {
      ret += children[i].min_width();
    }
    return Math.max(CELL_W * INTERNAL_CAPACITY, ret);
  }

  function rec_set_position(sx) {
    if (is_leaf()) {
      sx += LEAF_GAP;
      for (var b = next; b; b = b.next()) {
        b.set_position(sx, pos_y);
        sx += LEAF_GAP;
      }
    } else {
      var w = 0;
      for (var i = 0; i < children.length; i++) {
        children[i].set_position(sx + w, pos_y + LEVEL_GAP);
        children[i].rec_set_position(sx + w);
        w += children[i].min_width();
      }
      pos_x = sx + (w - CELL_W * INTERNAL_CAPACITY) / 2;
    }
  }

  function snapshot(s) {
    if (is_leaf()) {
      s.leaves = s.leaves || [];
      s.leaves.push({ x : pos_x, y : pos_y, has_next : !!next, id : id });
      for (var b = next; b; b = b.next()) {
        s.leaves.push({ x : b.x(), y : b.y(), has_next : !!b.next(), id : b.id() });
      }

      s.leaf_connectors = s.leaf_connectors || [];
      if (next) s.leaf_connectors.push({ x : pos_x, y : pos_y, id : id });
      for (var b = next; b && b.next(); b = b.next()) {
        s.leaf_connectors.push({ x : b.x(), y : b.y(), id : b.id() });
      }

      s.cell_values = s.cell_values || [];
      for (var i = 0; i < data.length; i++) s.cell_values.push({ x : pos_x, y : pos_y + CELL_H * i, sorted: thisRef.sorted(), value: data[i] });
      for (var b = next; b; b = b.next()) {
        var d = b.data();
        for (var i = 0; i < d.length; i++) s.cell_values.push({ x : b.x(), y : b.y() + CELL_H * i, sorted: b.sorted(), value: d[i] });
      }
    } else {
      s.internals = s.internals || [];
      s.internals.push({ x : pos_x, y : pos_y, id : id });
      for (var i = 0; i < children.length; i++) children[i].snapshot(s);

      s.cell_values = s.cell_values || [];
      for (var i = 0; i < data.length; i++) s.cell_values.push({ x : pos_x + CELL_W * i, y : pos_y, value: data[i], sorted: true });

      s.children_links = s.children_links || [];
      var sumw = 0;
      for (var i = 0; i < children.length; i++) {
        var xx = 0;
        if (children[i].is_leaf()) {
          xx = children[i].x() + 0.5 * CELL_W;
        } else {
          xx = children[i].x() + (i ? 0 : 1) * CELL_W * INTERNAL_CAPACITY;
        }
        s.children_links.push({ x1: pos_x + i * CELL_W, y1: pos_y + CELL_H, x2: xx,  y2: pos_y + LEVEL_GAP, child_id : children[i].id() });
        sumw += children[i].min_width();
      }
    }
    return s;
  }

  return thisRef;
}


function NewCTree(svg) {
  var root = NewBucket(null);
  root.set_position(CELL_W, CELL_H);

  function add(num) {
    leaf_insert(root, num);
  }

  function distribute_values(leafb, pivot, chain) {
    while (leafb.size()) {
      var i = leafb.last_data_is_at_least(pivot);
      leaf_insert(chain[i], leafb.leaf_promote_last());
    }
  }

  function leaf_insert(B, value) {
    if (!B.is_full()) return B.append(value);
    var tail = B.tail();
    assert(!tail || !tail.next());
    if (!tail || tail.is_full()) {
      var px = tail ? tail.x() : B.x();
      var py = tail ? tail.y() : B.y();
      add_chain(B, tail = NewBucket(null));
      tail.set_position(px + LEAF_GAP, py);
    }
    tail.append(value);
  }

  function add_chain(B, next) {
    if (!B.next()) {
      B.set_next(next);
      B.set_tail(next);
    } else {
      B.tail().set_next(next);
      B.set_tail(next);
    }
  }

  function pick_median(leafb) {
    var arr = [], last = 0, prev_last = 0;
    for (var b = leafb; b; prev_last = last, last = b, b = b.next())
      for (var i = 0; i < b.size(); i++)
        arr.push(b.data()[i]);
    arr.sort(function (a, b) { return a - b; });
    var pivot = arr[arr.length >> 1];
    for (var b = leafb; b; b = b.next())
      for (var i = 0; i < b.size(); i++)
        if (b.data()[i] == pivot) {
          if (b == last) {
            var t = b.remove_last();
            if (i == b.size()) {
              assert(t == pivot);
            } else {
              b.data()[i] = t;
            }
            if (last.size() == 0) {
              // free last
              assert(prev_last);
              if (leafb.next() == last) {
                leafb.set_next(null);
                leafb.set_tail(null);
              } else {
                prev_last.set_next(null);
                leafb.set_tail(prev_last);
              }
            }
          } else {
            b.data()[i] = last.remove_last();
          }
          return pivot;
        }
    assert(false);
  }

  function leaf_split(leafb) {
    // leafb.print_chain();

    var promotedValue, new_leafb = null;
    assert(leafb.next());

    var pivot = promotedValue = pick_median(leafb);
    new_leafb = NewBucket(leafb.parent());
    leafb.move_data_at_least(pivot, new_leafb);
    var chain = [ leafb, new_leafb ];
    // console.log('before split: ' + JSON.stringify(leafb.data()) + ' <> ' + JSON.stringify(new_leafb.data()) + ' pivot = ' + pivot);

    var Nb = leafb.detach_and_get_next();
    var Lb = null, Rb = null;
    var hi = [], nhi = 0;
    var lo = [], nlo = 0;
    while (true) {
      if (nhi && nlo) {
        assert(Lb && Rb);
        // console.log('fubefore ' + JSON.stringify(Lb.data()) + ' ' + JSON.stringify(Rb.data()) + ' ' + hi + " <> " + lo + ' -> ' + nhi + ' <> ' + nlo);
        var arr = Lb.fusion(Rb, hi, lo, nhi, nlo);
        // console.log('fuafter ' + JSON.stringify(Lb.data()) + ' ' + JSON.stringify(Rb.data()));
        nhi = arr[0]; nlo = arr[1];
        if (!nhi) { add_chain(chain[0], Lb); Lb = null; }
        if (!nlo) { add_chain(chain[1], Rb); Rb = null; }
      } else if (!Lb) {
        if (!Nb) break;
        Lb = Nb;
        Nb = Nb.detach_and_get_next();
        if (!Lb.is_full()) break;
      } else if (!nhi) {
        assert(Lb);
        nhi = Lb.mark_hi(pivot, hi, nhi);
        if (!nhi){ add_chain(chain[0], Lb); Lb = null; }
      } else if (!Rb) {
        if (!Nb) break;
        Rb = Nb;
        Nb = Nb.detach_and_get_next();
        if (!Rb.is_full()) break;
      } else if (!nlo) {
        assert(Rb);
        nlo = Rb.mark_lo(pivot, lo, nlo);
        if (!nlo){ add_chain(chain[1], Rb); Rb = null; }
      } else {
        assert(0);
      }
    }
    assert(!Nb);

    if (Lb) distribute_values(Lb, pivot, chain);
    if (Rb) distribute_values(Rb, pivot, chain);

    // console.log('split1: ' + JSON.stringify(leafb.data()) + ' <> ' + JSON.stringify(new_leafb.data()));

    consolidate(leafb);
    consolidate(new_leafb);

    // console.log('split2: ' + JSON.stringify(leafb.data()) + ' <> ' + JSON.stringify(new_leafb.data()));
    // leafb.print_chain();
    // new_leafb.print_chain();

    return [promotedValue, new_leafb];
  }

  function consolidate(leafb) {
    while (!leafb.is_full()) {
      if (leafb.tail()) {
        var x = leafb.tail().remove_last();
        leafb.append(x);
        if (!leafb.tail().size()) {
          var prev_last = 0, last = leafb;
          for (; last.next(); prev_last = last, last = last.next()) {}
          assert(last == leafb.tail());
          if (leafb.next() == last) {
            leafb.set_next(null);
            leafb.set_tail(null);
          } else {
            prev_last.set_next(null);
            leafb.set_tail(prev_last);
          }
        }
      } else {
        break;
      }
    }
  }

  function internal_insert(internalb, value, newb, left) {
    internalb.internal_insert(value, newb, left);
    newb.set_parent(internalb);
  }

  function internal_split(internalb) {
    var new_internalb = NewBucket(internalb.parent(), internalb.mid_child());
    internalb.move_half_to(new_internalb);
    for (var i = 0; i <= new_internalb.size(); i++) {
      new_internalb.child(i).set_parent(new_internalb);
    }
    return new_internalb;
  }

  function split_chain(leafb) {
    if (!leafb.next()) return false;

    var arr = leaf_split(leafb, promotedValue, nb);
    var promotedValue = arr[0];
    var nb = arr[1];
    var parent = leafb.parent();
    while (parent && nb) {
      if (parent.is_full()) {
        var inb = internal_split(parent);
        var promotedValueInternal = parent.remove_last();
        if (promotedValue >= promotedValueInternal) {
          internal_insert(inb, promotedValue, nb);
        } else {
          internal_insert(parent, promotedValue, nb);
        }
        promotedValue = promotedValueInternal;
        nb = inb;
        parent = nb.parent();
      } else {
        internal_insert(parent, promotedValue, nb);
        nb = false;
      }
    }
    if (nb) {
      root = NewBucket(null, root);
      internal_insert(root, promotedValue, nb);
      console.log("NEW ROOT " + promotedValue);
    }
    root.rec_set_position(CELL_W);
    return true;
  }

  // Returns <bucket, pos> if found in internal node,
  // Otherwise returns <bucket, splitted> for leaf node.
  function find_bucket(value, include_internal) {
    var b = root, splitted = 0;
    while (true) {
      if (b.is_leaf()) {
        if (!split_chain(b)) break;
        b = b.parent();
        assert(b);
        splitted = 1;
      } else {
        var pos = b.internal_lower_pos(value);
        if (include_internal && b.equal(pos, value)) {
          return [b, pos]; // Found in the internal bucket.
        }
        b = b.child(pos);    // Search the child.
      }
    }
    return [b, splitted];
  }

  function lower_bound(value) {
    var p = find_bucket(value, true);

    var ret = [false, 0];
    if (!p[0].is_leaf()) {
      // Found in internal bucket.
      ret = [true, value];
    } else {
      var pos = p[0].leaf_lower_pos(value);
      if (pos < p[0].size()) {
        ret = [true, p[0].data(pos)];
      } else {
        var b = p[0].parent();
        while (b) {
          pos = b.internal_lower_pos(value);
          if (pos < b.size()) {
            ret = [true, b.data(pos)];
            break;
          }
          b = b.parent();
        }
      }
    }
    return ret;
  }

  function snapshot() {
    return root.snapshot({});
  }

  return {
    add: add,
    lower_bound: lower_bound,
    snapshot: snapshot,
  };
}


/*
data = d3.csv.parse(data);

var base = false;
var fmt = d3.format(".3f");
data.forEach(function(d) {
  if (d.algorithm == 'comb_noup') base = d;
  d.insert_time = parseFloat(d.insert_time);
  d.query_time = parseFloat(d.query_time);
  d.total_time = d.insert_time + d.query_time;
});

data.sort(function (a, b) {
  return (a.algorithm < b.algorithm) ? -1 :
         (a.algorithm > b.algorithm) ?  1 : 0;
});

function build_table(N, Q) {
  var arr = [];
  data.forEach(function(d) {
    if (d.N != N || d.Q != Q) return;
    arr.push(d);
  });

  if (arr.length == 0) return '';

  arr.sort(function (a, b) { return a.total_time - b.total_time; });

  var s = '\
    <table border=1>\
      <caption>N = ' + N + ', Q = ' + Q + '</caption>\
      <thead>\
        <tr>\
          <th>Run Date</th>\
          <th>Hostname</th>\
          <th>Algorithm</th>\
          <th>Insert Time</th>\
          <th>Query Time</th>\
          <th>Total Time</th>\
          <th>Checksum</th>\
          <th>N Leaves</th>\
          <th>N Capacity</th>\
          <th>N Internals</th>\
          <th>Tree Depth</th>\
          <th>N Slacks</th>\
          <th>Internal Allocator</th>\
          <th>Leaf Allocator</th>\
          <th>Note</th>\
        </tr>\
      </thead>\
      <tbody>';
  for (var i = 0; i < arr.length; i++) {
    var d = arr[i];
    s += '\
    <tr>\
    <td>'+ d.timestamp +'</td>\
    <td>'+ d.hostname +'</td>\
    <td>'+ d.algorithm +'</td>\
    <td>'+ fmt(d.insert_time) +'</td>\
    <td>'+ fmt(d.query_time) +'</td>\
    <td>'+ fmt(d.total_time) +'</td>\
    <td>'+ d.checksum +'</td>\
    <td>'+ (d.nLeaves || 0) +'</td>\
    <td>'+ (d.nCap || 0) +'</td>\
    <td>'+ (d.nInternals || 0) +'</td>\
    <td>'+ (d.max_depth || 0) +'</td>\
    <td>'+ (d.slack || 0) +'</td>\
    <td>'+ (d.ia_free ? (d.ia_free + ' / ' + d.ia_size) : 0) +'</td>\
    <td>'+ (d.la_free ? (d.la_free + ' / ' + d.la_size) : 0) +'</td>\
    <td>'+ (d.version || '') +'</td>\
    </tr>\
    ';
  }
  return  s + '</tbody></table>'
}

var s = '';
for (var N = 100000000; N <= 100000000; N *= 10) {
  for (var Q = 1; Q <= 1e9; Q *= 10) {
    s += build_table(N, Q);
  }
}

document.write(s);

*/
