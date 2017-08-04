#include <algorithm>
#include <thread>

static void parallel_sort(long long arr[], int N) {
  int q = N / 4;
  long long *q0 = arr;
  long long *q1 = arr + q;
  long long *q2 = arr + 2 * q;
  long long *q3 = arr + 3 * q;
  long long *q4 = arr + N;

  std::thread t0([&]() { std::sort(q0, q1); });
  std::thread t1([&]() { std::sort(q1, q2); });
  std::thread t2([&]() { std::sort(q2, q3); });
  std::thread t3([&]() { std::sort(q3, q4); });
  t0.join();
  t1.join();
  t2.join();
  t3.join();

  std::thread t4([&]() { std::inplace_merge(q0, q1, q2); });
  std::thread t5([&]() { std::inplace_merge(q2, q3, q4); });
  t4.join();
  t5.join();

  std::inplace_merge(q0, q2, q4);
}

static void parallel_sort8(long long arr[], int N) {
  int q = N / 8;
  long long *q0 = arr;
  long long *q1 = arr + q;
  long long *q2 = arr + 2 * q;
  long long *q3 = arr + 3 * q;
  long long *q4 = arr + 4 * q;
  long long *q5 = arr + 5 * q;
  long long *q6 = arr + 6 * q;
  long long *q7 = arr + 7 * q;
  long long *q8 = arr + N;

  std::thread t0([&]() { std::sort(q0, q1); });
  std::thread t1([&]() { std::sort(q1, q2); });
  std::thread t2([&]() { std::sort(q2, q3); });
  std::thread t3([&]() { std::sort(q3, q4); });
  std::thread t4([&]() { std::sort(q4, q5); });
  std::thread t5([&]() { std::sort(q5, q6); });
  std::thread t6([&]() { std::sort(q6, q7); });
  std::thread t7([&]() { std::sort(q7, q8); });
  t0.join();
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  t5.join();
  t6.join();
  t7.join();

  std::thread t8([&]() { std::inplace_merge(q0, q1, q2); });
  std::thread t9([&]() { std::inplace_merge(q2, q3, q4); });
  std::thread t10([&]() { std::inplace_merge(q4, q5, q6); });
  std::thread t11([&]() { std::inplace_merge(q6, q7, q8); });
  t8.join();
  t9.join();
  t10.join();
  t11.join();

  std::thread t12([&]() { std::inplace_merge(q0, q2, q4); });
  std::thread t13([&]() { std::inplace_merge(q4, q6, q8); });
  t12.join();
  t13.join();

  std::inplace_merge(q0, q4, q8);
}
