#include <cassert>
#include <cstdio>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp[BATCH] = { 0 };
int wtmp[BATCH];
map<int, int> cnt;

int main() {
  FILE *in = fopen("data/skyserver.data", "rb");
  FILE *out = fopen("data/skyserver.u5data", "wb");
  assert(in);
  assert(out);
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int ntmp = 0;
    for (int i = 0; i < N; i++, nth++) {
      int t = tmp[i] * 5;
      int ninc = 0;
      while (cnt.count(t)) {
        t++;
        ninc++;
      }
      if (ninc > 10000) {
        fprintf(stderr, "?");
      } else if (ninc > 1000) {
        fprintf(stderr, ";");
      } else if (ninc > 100) {
        fprintf(stderr, ".");
      }
      cnt[t] = 1;
      wtmp[ntmp++] = t;
    }
    fwrite(wtmp, sizeof(int), ntmp, out);
    fflush(out);
    fprintf(stderr, "%d. size = %lu.\n", nth, cnt.size());
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
