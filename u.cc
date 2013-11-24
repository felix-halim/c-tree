#include <cstdio>

#include <set>

using namespace std;

#define BATCH 1000000

int tmp[BATCH] = { 0 };
int wtmp[BATCH];
set<int> s;

int main() {
  FILE *in = fopen("rb", "data/skyserver.data");
  FILE *out = fopen("wb", "data/skyserver.udata");
  assert(in);
  assert(out);
  int ndup = 0;
  while (!feof(in)) {
    int N = fread(tmp, sizeof(int), BATCH, in);
    int nt = 0;
    for (int i = 0; i < N; i++) {
      if (s.count(tmp[i])) {
        ndup++;
      } else {
        s.insert(tmp[i]);
        wtmp[nt++] = tmp[i];
      }
    }
    if (nt > 0) {
      fwrite(wtmp, sizeof(int), nt, out);
    }
    fprintf(stderr, "ndup = %d\n", ndup);
  }
  fclose(out);
  if (ferror(in)) { fprintf(stderr,"Error reading file!\n"); }
  if (feof(in)) { fclose(in); in = NULL; }
}
