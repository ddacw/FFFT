#include <iostream>

#include "fft.h"
#include "mfft.h"
// #include "util.h"

int main(int argc, char* argv[]) {
  jarray X;
  size_t num_threads = 2;
  for (int i = 0; i < 64; ++i) {
    X.emplace_back(double(i % 8) + 0.1 * (i % 8));
  }
  std::cerr << std::endl;

  jarray XP(X);

  MFFT::Transform(X, false, 3);
  MFFT::Transform(X, true, 3);
  FFT::FFTSeq(XP, 0);

  std::cout << std::endl;
  for (int i = 0; i < 64; ++i) {
    std::cout << X[i] << ' ' << XP[i] << std::endl;
  }
}
