# Crack Tree (C-Tree) [![Build Status](https://travis-ci.org/felix-halim/c-tree.svg?branch=master)](https://travis-ci.org/felix-halim/c-tree)

Crack Tree (C-Tree) is a container data structure inspired from Database Cracking philosophy (i.e., always do just enough). Like Cracking, C-Tree incrementally sorts the elements as a side effect of query processing. Unlike Cracking, C-Tree first laid out the elements in buckets that are chained together like a linked list and then later transitioned into a B-Tree (or other trees) like structure as the elements are touched by queries. See [this demo](https://felix-halim.github.io/c-tree/docs/trimmer.html).

Advantages of C-Tree over the original single array Cracking:
- Updates are faster since it avoids the need to shift/move existing elements to make space.
- Buckets can be partitioned more effectively, chained to other buckets, minimizing data movements.

Advantages of C-Tree over other container data structures (e.g., B-Tree, AVL-Tree, ART)
- C-Tree amortized cost is always lower under any number of query.
The amortized cost is calculated from when the elements are inserted to the container
until the results of the queries are produced.

C-Tree leverages modern CPU optimizations such as SIMD for the partitioning algorithm.
C-Tree partitioning algorithm takes two buckets L and R,
marks elements in L to be swapped with R and vice versa,
then it swaps only the marked elements.
With this optimization, C-Tree is able to sort 10^8 random 64 bit signed integer
about 20-30 percent faster than the standard STL sort.
See [ctree_sort.cc](https://github.com/felix-halim/c-tree/blob/master/src/ctree_sort.cc).

See [the comparisons](https://felix-halim.github.io/c-tree/docs/graphs.html)
of C-Tree vs other data structures on point and range sum queries,
with interleaved updates, varying selectivity and skewness.

Note that while C-Tree has the lowest amortized cost on any number of query,
the first query may take significantly longer.
This issue can be avoided by issuing several "warming up" queries to the C-Tree
before the first query.
