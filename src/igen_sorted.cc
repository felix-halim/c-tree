#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.h"

using namespace std;

/**
 * Generates a sorted input, perform a random sum query.
 */
int main(int argc, char *argv[]) {
  assert(argc == 3);

  int n = atoi(argv[1]);
  for (int i = 0; i < n; i++) {
    printf("%d %d\n", INSERT, i + 1);
  }

  int q = atoi(argv[2]);
  for (int i = 0; i < q; i++) {
    int j = rand() % n;
    printf("%d %d %d\n", SUM, j, j + 1);
  }
}
