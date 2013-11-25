#include <cassert>
#include <cstdio>
#include <cstring>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp1[BATCH];

int main() {
  FILE *in1 = fopen("data/skyserver.data", "rb");
  assert(in1);
  int mx = 0;
  for (int nth = 0; !feof(in1); nth++) {
    int N1 = fread(tmp1, sizeof(int), BATCH, in1);
    for (int i = 0; i < N1; i++) {
      mx = max(mx, tmp1[i]);
    }
    fprintf(stderr, "nth = %d, %d\n", nth, mx);
  }
}
