#include <algorithm>
#include <thread>

void parallel_sort(long long arr[], long long tmp[], int N) {
  int q = N / 4;
  std::thread t0([&]() { std::sort(arr, arr + q); });
  std::thread t1([&]() { std::sort(arr + q, arr + q * 2); });
  std::thread t2([&]() { std::sort(arr + q * 2, arr + q * 3); });
  std::thread t3([&]() { std::sort(arr + q * 3, arr + N); });
  t0.join();
  t1.join();
  t2.join();
  t3.join();

  std::thread t4(
      [&]() { std::merge(arr, arr + q, arr + q, arr + 2 * q, tmp); });
  std::thread t5([&]() {
    std::merge(arr + 2 * q, arr + 3 * q, arr + 3 * q, arr + N, tmp + 2 * q);
  });
  t4.join();
  t5.join();

  std::merge(tmp, tmp + 2 * q, tmp + 2 * q, tmp + N, arr);
}
