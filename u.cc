#include <cassert>
#include <cstdio>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

int tmp[BATCH] = { 0 };
int wtmp[BATCH];
map<int, int> next;

int main() {
  FILE *in = fopen("data/skyserver.data", "rb");
  FILE *out = fopen("data/skyserver.udata", "wb");
  assert(in);
  assert(out);
  int ndup = 0;
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int nt = 0;
    for (int i = 0; i < N; i++, nth++) {
      int t = tmp[i];
      if (next.count(t)) {
        t = next[t];
        ndup++;

        int ninc = 0;
        while (next.count(t)) t = next[t], ninc++;
        if (ninc > 1) {
          fprintf(stderr, ".");
        } else if (ninc > 10) {
          fprintf(stderr, "x");
        } else if (ninc > 100) {
          fprintf(stderr, "z");
        }
      }
      assert(!next.count(t) && t >= 0);
      next[tmp[i]] = t + 1;
      wtmp[nt++] = t;
    }
    fwrite(wtmp, sizeof(int), nt, out);
    fflush(out);
    fprintf(stderr, "%d. ndup = %d, size = %lu.\n", nth, ndup, next.size());
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
