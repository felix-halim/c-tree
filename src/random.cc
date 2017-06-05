#include <cstdio>
#include <cassert>

#include "time_it.h"
#include "random.h"

int main() {
  int n = 100000000;

  fprintf(stderr, "%.3lf\n", time_it([&]() {
    Random r(140384);
    long long csum = 0;
    for (int i = 0; i < n; i++) {
      uint64_t num = r.nextLong();
      // printf("%lu\n", num);
      csum += num;
    }
    printf("%20lld\n", csum);
  }));
}
