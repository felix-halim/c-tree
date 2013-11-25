#include <cassert>
#include <cstdio>
#include <cstring>

#include <set>
#include <map>
#include <algorithm>

using namespace std;

#define BATCH 10000000

unsigned int tmp1[BATCH];
set<unsigned int> s;

int main() {
  FILE *in1 = fopen("data/skyserver.u5data", "rb");
  assert(in1);
  for (int nth = 0; !feof(in1); nth++) {
    fprintf(stderr, "nth = %d\n", nth);
    int N1 = fread(tmp1, sizeof(unsigned int), BATCH, in1);
    for (int i = 0; i < N1; i++) {
      assert(!s.count(tmp1[i]));
      s.insert(tmp1[i]);
    }
  }
}
