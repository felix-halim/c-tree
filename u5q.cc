#include <cassert>
#include <cstdio>

#include <set>
#include <map>
#include <algorithm>

#include "src/comb.h"

using namespace std;

#define BATCH 10000000

int tmp[BATCH] = { 0 };
unsigned wtmp[BATCH];

int main() {
  fprintf(stderr, "%d\n", sizeof(unsigned));
  FILE *in = fopen("data/skyserver.queries", "rb");
  FILE *out = fopen("data/skyserver.u5queries", "wb");
  assert(in);
  assert(out);
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int ntmp = 0;
    for (int i = 0; i < N; i++, nth++) {
      unsigned t = ((unsigned) tmp[i]) * 10;
      wtmp[ntmp++] = t;
    }
    fwrite(wtmp, sizeof(unsigned), ntmp, out);
    fflush(out);
    fprintf(stderr, "%d. size = %lu.\n", nth, cnt.size());
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
