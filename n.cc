#include <stdio.h>

int main() {
  int m = 1238981;
  for (long long n = 2; n <= m; n++) {
    if ((n * n) % m == n) {
      fprintf(stderr, "n = %lld\n", n);
      break;
    }
  }
}
