#include "fft.h"

// #include "gfft.h"
// #include "util.h"

int main(int argc, char* argv[]) {
  jarray X;
  size_t num_threads = 2;
  for (int i = 0; i < (1 << 24); ++i) {
    X.emplace_back(double(i % 8));
  }
  FFT::FFTParallel(X, 0, num_threads);

  for (int i = 0; i < 10; ++i) {
    std::cout << X[i] << ' ';
  }
  std::cout << std::endl;
  FFT::FFTParallel(X, 1, num_threads);
  for (int i = 0; i < 10; ++i) {
    std::cout << X[i] << ' ';
  }
  // GeneralizedFFT pfft(X, n_threads);
  // pfft.Print();
}
