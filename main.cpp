#include <cmath>
#include <iostream>

#include "fft.h"
#include "mfft.h"

const double eps = 1e-5;
typedef std::pair<int, cmplx> ii;
typedef std::vector<double> farray;

std::vector<ii> Compress(const std::vector<double>& D) {
  jarray X;
  for (double d : D) {
    X.emplace_back(d);
  }
  FFT::FFTParallel(X, false, 4);

  std::vector<ii> zip;
  zip.emplace_back(X.size(), 0);
  for (size_t i = 0; i < X.size(); ++i) {
    if (abs(X[i]) > eps) {
      zip.push_back({i, X[i]});
    }
  }
  std::cout << "Compression factor: " << zip.size() * 3. / D.size()
            << std::endl;
  return zip;
}

std::vector<double> Extract(const std::vector<ii>& zip) {
  size_t sz = zip[0].first;
  jarray X(sz);
  for (size_t i = 1; i < zip.size(); ++i) {
    X[zip[i].first] = zip[i].second;
  }
  FFT::FFTParallel(X, true, 4);
  std::vector<double> D;
  for (cmplx x : X) {
    D.push_back(x.real());
  }
  return D;
}

farray PolyMult(const farray& A, const farray& B) {
  jarray AX, BX;
  for (auto x : A) {
    AX.emplace_back(x);
  }
  for (auto x : B) {
    BX.emplace_back(x);
  }
  int n_bit = CeilBit(A.size() + B.size());
  size_t common_size = (1 << n_bit);
  while (AX.size() < common_size) {
    AX.emplace_back(0);
  }
  while (BX.size() < common_size) {
    BX.emplace_back(0);
  }
  MFFT::Transform(AX, false, 4);
  MFFT::Transform(BX, false, 4);
  for (size_t i = 0; i < BX.size(); ++i) {
    BX[i] *= AX[i];
  }
  MFFT::Transform(BX, true, 4);
  farray ret;
  for (cmplx x : BX) {
    ret.push_back(x.real());
  }
  while (ret.size() > 1 && ret.back() < 1e-9) {
    ret.pop_back();
  }
  return ret;
}

int main(int argc, char* argv[]) {
  // Test Compress
  // Proof of concept, since the performance is abysmal with the weather
  // data.
  std::vector<double> D;
  for (size_t i = 0; i < 128; ++i) {
    D.push_back(cos(i % 8));
  }
  auto zip = Compress(D);
  auto unzip = Extract(zip);

  // Test PolyMult.
  farray A = {1, 3, 5, 7, 9};
  farray B = {2, 4, 6, 8};
  farray AB = PolyMult(A, B);
  for (size_t i = 0; i < AB.size(); ++i) {
    std::cerr << AB[i] << ' ';
  }
}
