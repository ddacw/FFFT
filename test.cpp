#include <fstream>

#include "fft.h"

int main(int argc, char* argv[]) {
  // std::ifstream infile;
  // infile.open("data/temp.txt");
  // jarray X_orig;
  // double x;
  // while (infile >> x) {
  //   X_orig.push_back(x);
  // }

  // std::ofstream rec;
  // jarray X(X_orig);

  // rec.open("fft_cpp_rec.txt");
  // FFT::FFTRec(X);
  // for (size_t i = 0; i < X.size(); ++i) {
  //   rec << X[i] << std::endl;
  // }
  // rec.close();

  // X = X_orig;
  // rec.open("fft_cpp_seq.txt");
  // FFT::FFTSeq(X, 0);
  // for (size_t i = 0; i < X.size(); ++i) {
  //   rec << X[i] << std::endl;
  // }
  // rec.close();

  // size_t num_threads = 8;
  // X = X_orig;
  // rec.open("fft_cpp_par.txt");
  // FFT::FFTParallel(X, false, num_threads);
  // for (size_t i = 0; i < X.size(); ++i) {
  //   rec << X[i] << std::endl;
  // }
  // rec.close();
}
