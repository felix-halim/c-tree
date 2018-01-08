#include "experiments/ctree.h"
#include "experiments/ska_sort.h"
#include "experiments/vergesort.h"
#include "random.h"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <vector>

static std::vector<long long> arr;

static void init_arr() {
  if (arr.size()) return;
  Random random(14);
  for (int i = 0; i < 100000; i++) {
    arr.push_back(random.nextLong());
  }
}

static void CTree_Sort(benchmark::State& state) {
  init_arr();
  for (auto _ : state) {
    state.PauseTiming();
    std::vector<long long> a = arr;
    state.ResumeTiming();
    ctreesort<1024>(a.data(), a.size());
  }
}
BENCHMARK(CTree_Sort);

static void Std_Sort(benchmark::State& state) {
  init_arr();
  for (auto _ : state) {
    state.PauseTiming();
    std::vector<long long> a = arr;
    state.ResumeTiming();
    std::sort(a.begin(), a.end());
  }
}
BENCHMARK(Std_Sort);

static void Verge_Sort(benchmark::State& state) {
  init_arr();
  for (auto _ : state) {
    state.PauseTiming();
    std::vector<long long> a = arr;
    state.ResumeTiming();
    vergesort::vergesort(a.begin(), a.begin() + a.size());
  }
}
BENCHMARK(Verge_Sort);

static void Ska_Sort(benchmark::State& state) {
  init_arr();
  for (auto _ : state) {
    state.PauseTiming();
    std::vector<long long> a = arr;
    state.ResumeTiming();
    ska_sort(a.begin(), a.begin() + a.size());
  }
}
BENCHMARK(Ska_Sort);

BENCHMARK_MAIN();