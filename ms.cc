#include <cstdio>
#include <memory>
#include <vector>

using namespace std;

class A {
 public:
  int x;

  A() {
    fprintf(stderr, "A created\n");
  }

  ~A() {
    fprintf(stderr, "A destroyed\n");
  }
};

void sptr() {
  vector<shared_ptr<A>> v;
  v.push_back(make_shared<A>());

  fprintf(stderr, "cnt = %ld\n", v[0].use_count());

  v.pop_back();
  fprintf(stderr, "cnt\n");
}

int main() {
  // sptr();

  unique_ptr<A> u { new A() };

  unique_ptr<A> x = move(u);
}
