#include <cassert>
#include <cstdio>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp[BATCH] = { 0 };
long long wtmp[BATCH];
map<int, int> cnt;

int main() {
  FILE *in = fopen("data/skyserver.data", "rb");
  FILE *out = fopen("data/skyserver.udata", "wb");
  assert(in);
  assert(out);
  int ndup = 0;
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int ntmp = 0;
    for (int i = 0; i < N; i++, nth++) {
      long long t = tmp[i];
      assert(cnt[t] < 100);
      t = t * 100 + cnt[t]++;
      wtmp[ntmp++] = t;
    }
    fwrite(wtmp, sizeof(long long), ntmp, out);
    fflush(out);
    fprintf(stderr, "%d. size = %lu.\n", nth, cnt.size());
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
