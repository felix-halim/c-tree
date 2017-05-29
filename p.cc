#include <cstdio>
#include <cstdlib>

using namespace std;

template <typename T>
class A {
public:
  T X;
  void sx() { X = 10; }
};

template <typename T>
class B : A<T> {
public:
  B():A<T>() {
    A<T>::sx();
  }

  T x() {
    return A<T>::X;
  }
};

int main () {
  B<int> b;
  return b.x();
}
