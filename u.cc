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
map<int, int> root;

int get_root(int t) {
  if (!root.count(t)) return root[t] = t;
  int nr = get_root(root[t]);
  return root[t] = nr;
}

int main() {
  FILE *in = fopen("data/skyserver.data", "rb");
  FILE *out = fopen("data/skyserver.udata2", "wb");
  assert(in);
  assert(out);
  int ndup = 0;
  for (int nth = 0; !feof(in); ) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int ntmp = 0;
    for (int i = 0; i < N; i++, nth++) {
      int r = get_root(tmp[i]);
      fprintf(stderr, "i = %d, r = %d\n", nth, r);
      int t = next.count(r) ? next[r] : r;
      fprintf(stderr, "i1 = %d, r = %d\n", nth, r);
      assert(t >= 0);
      int nr = get_root(t + 1);
      fprintf(stderr, "i2 = %d, r = %d\n", nth, r);
      int nt = next.count(nr) ? next[nr] : nr;
      fprintf(stderr, "i3 = %d, r = %d\n", nth, r);
      root[nr] = r;
      next[r] = nt;
      wtmp[ntmp++] = t;
      fprintf(stderr, "i4 = %d, r = %d\n", nth, r);
    }
    fwrite(wtmp, sizeof(int), ntmp, out);
    fflush(out);
    fprintf(stderr, "%d. ndup = %d, size = %lu.\n", nth, ndup, next.size());
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
