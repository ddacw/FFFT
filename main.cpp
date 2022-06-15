#include <iostream>

#include "fft.h"
#include "mfft.h"

int main(int argc, char* argv[]) {
  jarray X;
  for (int i = 0; i < 25; ++i) {
    X.emplace_back(double(i) + 0.1 * (i % 8));
  }
  

  jarray XP(X);

  MFFT::Transform(X, false, 5);
  // MFFT::Transform(X, true, 2);
  // MFFT::Transform(X, true, 4);
  FFT::FFTSlow(XP);

  std::cout << std::endl;
  for (int i = 0; i < 8; ++i) {
    std::cout << X[i] << ' ' << XP[i] << std::endl;
  }
}
