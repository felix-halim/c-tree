#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../common.h"
#include "../random.h"

// Generates a sorted input, perform a random sum query.
void gen(int n, int q) {
  assert(argc == 3);

  for (int i = 0; i < n; i++) {
    printf("%d %d\n", INSERT, i + 1);
  }

  Random r(140384);
  for (int i = 0; i < q; i++) {
    int j = r.nextInt(n);
    printf("%d %d %d\n", SUM, j, j + 1);
  }
}
