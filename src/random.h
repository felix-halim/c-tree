#ifndef RANDOM_H
#define RANDOM_H

#include <cassert>

// Extracted from java.util.Random
class Random {
  static const long long multiplier = 0x5DEECE66DLL;
  static const long long addend = 0xBLL;
  static const long long mask = (1LL << 48) - 1;
  unsigned long long seed;

public:
  Random(long long s){ setSeed(s); }

  void setSeed(long long s){ seed = (s ^ multiplier) & mask; }

  /**
   * This is a linear congruential pseudorandom number generator, as
   * defined by D. H. Lehmer and described by Donald E. Knuth in
   * <i>The Art of Computer Programming,</i> Volume 3:
   * <i>Seminumerical Algorithms</i>, section 3.2.1.
   */
  int next(int bits) {
    seed = (seed * multiplier + addend) & mask;
    return (int)(seed >> (48 - bits));
  }

  /**
   * @param bound the upper bound (exclusive).  Must be positive.
   * @return the next pseudorandom, uniformly distributed {@code int}
   *         value between zero (inclusive) and {@code bound} (exclusive)
   *         from this random number generator's sequence
   */
  int nextInt(int bound) {
    assert(bound > 0);
    int r = next(31);
    int m = bound - 1;
    if ((bound & m) == 0)  // i.e., bound is a power of 2
        r = (int)((bound * (long)r) >> 31);
    else {
        for (int u = r;
             u - (r = u % bound) + m < 0;
             u = next(31))
            ;
    }
    return r;
  }

  int nextInt(){
    return next(32);
  }

  long long nextLong(){
    return ((long long)(next(32)) << 32) + next(32);
  }

  bool nextBoolean(){ return next(1) != 0; }

  float nextFloat(){ return next(24) / ((float)(1 << 24)); }

  double nextDouble(){
    return (((long long)(next(26)) << 27) + next(27)) / (double)(1LL << 53);
  }
};

#endif
