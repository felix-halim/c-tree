#include <cassert>
#include <cstdio>
#include <cstring>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp1[BATCH];
int tmp2[BATCH];

int main() {
  FILE *in1 = fopen("data/skyserver.udata", "rb");
  FILE *in2 = fopen("data/skyserver.udata2", "rb");
  assert(in1);
  assert(in2);
  for (int nth = 0; !feof(in1) && !feof(in2); nth++) {
    fprintf(stderr, "nth = %d\n", nth);
    int N1 = fread(tmp1, sizeof(int), BATCH, in1);
    int N2 = fread(tmp2, sizeof(int), BATCH, in2);
    assert(N1 == N2);
    for (int i = 0; i < N1; i++) {
      if (tmp1[i] != tmp2[i])
        fprintf(stderr, "%d. %d %d\n", i, tmp1[i], tmp2[i]);
    }
    int res = memcmp(tmp1, tmp2, sizeof(int) * N1);
    assert(res == 0);
  }
}
