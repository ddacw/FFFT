#include <complex>
#include <iostream>
#include <utility>
#include <vector>

typedef std::complex<double> cmplx;

int CeilBit(int n) {
  int pow2 = 1, count = 0;
  while (pow2 < n) {
    pow2 <<= 1;
    count ++;
  }
  return count;
}

int GetBit(int mask, int n) {
  return (mask >> n) & 1;
}

int RevBit(int mask, int count) {
  for (int i = 0, j = count - 1; i < j; ++i, --j) {
    if (GetBit(mask, i) != GetBit(mask, j)) {
      mask ^= (1 << i);
      mask ^= (1 << j);
    }
  }
  return mask;
}

void RefFFT(std::vector<cmplx>& X, bool invert) {
  if (X.size() == 1) {
    return;
  } 

  int n_bit = CeilBit(X.size());
  int n = (1 << n_bit);
  while (X.size() < size_t(n)) {
    X.emplace_back(0);
  }

  for (int i = 0; i < n; ++i) {
    int j = RevBit(i, n_bit);
    if (i < j) {
      std::swap(X[i], X[j]);
    }
  }

  std::vector<cmplx> next(n);
  for (int step = 1; step < n; step <<= 1) {
    double angle = M_PI / step;
    if (invert) {
      angle = -angle;
    }
    cmplx w(1), wn(cos(angle), sin(angle));
    std::vector<cmplx> W;
    for (int i = 0; i < step; ++i) {
      W.push_back(w);
      w *= wn;
    }

    int start_even = 0;
    int start_odd = start_even + step;

    while (start_even < n) {
      for(int i = 0; i < step; ++i) {
        next[start_even + i] = X[start_even + i] + W[i] * X[start_odd + i];
        next[start_odd + i]  = X[start_even + i] - W[i] * X[start_odd + i];
      }
      start_even += 2*step;
      start_odd = start_even + step;
    }
    for (int i = 0; i < n; ++i) {
      X[i] = next[i];
    }
  }
  if (invert) {
    for(int i = 0; i < n; ++i) {
      X[i] /= n;
    }
  }
}

int main(int argc, char* argv[]) {
  std::vector<cmplx> X;

  for(int i = 0; i < 8; ++i) {
    X.emplace_back(double(i));
  }

  RefFFT(X, 0);
  // RefFFT(X, 1);
  for(int i = 0; i < 8; ++i) {
    std::cout << X[i] << ' ';
  }
  std::cout << std::endl;
}
