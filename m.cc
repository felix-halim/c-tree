#include <cstdio>
#include <cstdlib>
#include <map>

using namespace std;

int main() {
  map<int, int> m;
  for (int i = 0; i < 100000; i++) {
    m[rand()] = rand();
  }

  int cnt = 0;
  for (int i = 0; i < 10000000; i++) {
    auto it = m.lower_bound(rand());
    if (it != m.end()) cnt += it->second;
  }
  printf("%d\n", cnt);
}
