#include <cstdio>
#include <set>

#include "tester.h"

std::multiset<long long> s;

// op = 1: inserts the value.
void insert(long long value) {
  s.insert(value);
}

// op = 2: deletes the value. The value guaranteed to exists.
bool erase(long long value) {
  auto it = s.lower_bound(value);
  if (it != s.end() && *it == value) {
    s.erase(it);
    return true;
  }
  return false;
}

// op = 3: count values in range [a, b).
long long count(long long a, long long b) {
  auto it = s.lower_bound(a);
  long long cnt = 0;
  while (it != s.end() && *it < b) {
    cnt++;
    it++;
  }
  return cnt;
}

// op = 4: sum values in range [a, b).
long long sum(long long a, long long b) {
  auto it = s.lower_bound(a);
  long long sum = 0;
  while (it != s.end() && *it < b) {
    sum += *(it++);
  }
  return sum;
}
