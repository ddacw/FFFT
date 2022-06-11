#include "util.h"

cmplx GetW(int N) {
  double angle = 2. * M_PI / double(N);
  return exp(cmplx(0, -1) * angle);
  // return cmplx(cos(angle), sin(angle));
}

int CeilBit(int n) {
  int pow2 = 1, count = 0;
  while (pow2 < n) {
    pow2 <<= 1;
    count++;
  }
  return count;
}

int GetBit(int mask, int n) { return (mask >> n) & 1; }

int RevBit(int mask, int count) {
  for (int i = 0, j = count - 1; i < j; ++i, --j) {
    if (GetBit(mask, i) != GetBit(mask, j)) {
      mask ^= (1 << i);
      mask ^= (1 << j);
    }
  }
  return mask;
}

int PadZero(jarray& X) {
  int n_bit = CeilBit(X.size());
  int n = (1 << n_bit);
  while (X.size() < size_t(n)) {
    X.emplace_back(0);
  }
  return n_bit;
}
