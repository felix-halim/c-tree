#include <cassert>
#include <cstdio>

#include <algorithm>
#include <random>

using namespace std;

int main(int argc, char *argv[]) {
  int N;
  mt19937 gen(28284);
  uniform_int_distribution<> dis;
  sscanf(argv[1], "%d", &N);

  // Produces unique sparse positive 2 * N integers.
  int gap = 2047483647 / N;
  int *arr = new int[N * 2];
  arr[0] = dis(gen) % gap + 1;
  for (int i = 1; i < N * 2; i++)
    arr[i] = arr[i - 1] + dis(gen) % gap + 1;

  shuffle(arr, arr + N * 2, gen);

  char fn[100];
  sprintf(fn, "%d.data", N);
  FILE *out = fopen(fn, "wb");
  int nw = fwrite(arr, sizeof(int), N * 2, out);
  if (nw != N * 2) fprintf(stderr, "Written %d, expected: %d\n", nw, N * 2);
  fclose(out);
}
