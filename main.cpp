#include <iostream>

#include "fft.h"
#include "mfft.h"

int main(int argc, char* argv[]) {
  jarray X;
  for (int i = 0; i < 32; ++i) {
    X.emplace_back(double(i) + 0.1 * (i % 8));
  }

  jarray XP(X);

  MFFT::Transform(X, false, 4);
  // MFFT::Transform(X, true, 2);
  // MFFT::Transform(X, true, 2);
  // MFFT::Transform(X, true, 4);
  FFT::FFTParallel(XP, false, 2);
  // FFT::DFT(XP, true);
  std::cout << std::endl;
  for (int i = 0; i < 8; ++i) {
    std::cout << X[i] << ' ' << XP[i] << std::endl;
  }
}
