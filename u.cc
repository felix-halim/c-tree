#include <cassert>
#include <cstdio>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp[BATCH] = { 0 };
int wtmp[BATCH];
map<int, int> s;

int main() {
  FILE *in = fopen("data/skyserver.udata", "rb");
  FILE *out = fopen("data/skyserver.udata2", "wb");
  assert(in);
  assert(out);
  int ndup = 0, dup_max = 0;
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int nt = 0;
    for (int i = 0; i < N; i++, nth++) {
      int t = tmp[i];
      // printf("%d ", tmp[i]);
      if (s.count(t)) {
        s[t]++;
        dup_max = max(dup_max, s[t]);
        ndup++;

        int ninc = 0;
        while (s.count(t)) t++, ninc++;
        if (ninc > 100) {
          fprintf(stderr, ".");
        } else if (ninc > 1000) {
          fprintf(stderr, "x");
        } else if (ninc > 10000) {
          fprintf(stderr, "z");
        }
      }
      assert(!s.count(t) && t >= 0);
      s[t] = 1;
      wtmp[nt++] = t;
    }
    assert(nt == BATCH);
    fwrite(wtmp, sizeof(int), nt, out);
    fflush(out);
    fprintf(stderr, "%d. ndup = %d, size = %lu, dup_max = %d.\n", nth, ndup, s.size(), dup_max);
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
