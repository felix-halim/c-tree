# Benchmark tool for Indexing Data Structures [![Build Status](https://travis-ci.org/felix-halim/indexing-benchmark.svg?branch=master)](https://travis-ci.org/felix-halim/indexing-benchmark)

Tool for benchmarking data structures that supports range-sum query and updates. The benchmark first inserts random 10^8 64 signed integers and then performs range-sum queries that may be interleaved with updates. The queries are varied by selectivity and skewness.
