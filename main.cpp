#include <iostream>

#include "fft.h"
#include "mfft.h"
// #include "util.h"

int main(int argc, char* argv[]) {
  jarray X;
  size_t num_threads = 2;
  for (int i = 0; i < 64; ++i) {
    X.emplace_back(double(i % 8));
    // std::cerr << X.back() << ' ';
  }
  std::cerr << std::endl;

  jarray XP(X);

  MFFT mfft(X, 3);
  mfft.Transform();

  FFT::FFTSeq(XP, 0);

  std::cout << std::endl;
  for (int i = 0; i < 64; ++i) {
    std::cout << mfft.X[i] << ' ' << XP[i] << std::endl;
  }
}
