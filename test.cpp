#include <chrono>
#include <cmath>
#include <fstream>
#include <random>

static std::uniform_int_distribution<> uniform(-1, RAND_MAX);

#include "fft.h"
#include "mfft.h"

const double eps_small = 1e-8;

// template <typename T>
// int benchmark_multi_thread(T& bst, int count, int num_threads) {
//   auto start = std::chrono::steady_clock::now();
//   size_t block_size = count / num_threads;
//   std::vector<std::thread> workers(num_threads - 1);
//   int remain = count;
//   for (int i = 0; i < num_threads - 1; ++i) {
//     workers[i] = std::thread(&add_random<T>, std::ref(bst), block_size, i);
//     remain -= block_size;
//   }
//   add_random(bst, remain, num_threads - 1);
//   for (int i = 0; i < num_threads - 1; ++i) {
//     workers[i].join();
//   }
//   auto finish = std::chrono::steady_clock::now();
//   auto elapsed =
//       std::chrono::duration_cast<std::chrono::microseconds>(finish - start)
//           .count();
//   return elapsed;
// }

bool CorrectnessAux(std::function<void(jarray&, bool, size_t)> fft, int N,
                    int num_threads) {
  auto engine = std::default_random_engine(N);
  jarray X_orig(N);
  for (int i = 0; i < N; ++i) {
    X_orig[i] = cos(double(i)) + cos(uniform(engine)) + sin(i + 1.) +
                sin(uniform(engine) + 1);
  }
  jarray X(X_orig);

  fft(X, false, num_threads);
  fft(X, true, num_threads);

  double error = 0;
  for (int i = 0; i < N; ++i) {
    error += abs(X[i] - X_orig[i]);
  }

  return error < eps_small;
}

void TestCorrectness(std::function<void(jarray&, bool, size_t)> fft,
                     const std::string& fft_name) {
  bool verdict = true;
  std::cerr << "0 ";
  verdict &= CorrectnessAux(fft, 512, 1);
  std::cerr << "1 ";
  verdict &= CorrectnessAux(fft, 128, 4);
  std::cerr << "2 ";
  verdict &= CorrectnessAux(fft, 511, 2);
  std::cerr << "3 ";
  verdict &= CorrectnessAux(fft, 305, 4);
  std::cerr << "4 ";
  verdict &= CorrectnessAux(fft, 277, 8);
  if (!verdict) {
    std::cerr << fft_name << " is incorrect !!!" << std::endl;
  } else {
    std::cerr << fft_name << " passed." << std::endl;
  }
}

int main(int argc, char* argv[]) {
  // Correctness.
  TestCorrectness(FFT::DFT, "FFT::DFT");
  TestCorrectness(FFT::FFTRec, "FFT::FFTRec");
  TestCorrectness(FFT::FFTSeq, "FFT::FFTSeq");
  TestCorrectness(FFT::FFTParallel, "FFT::FFTParallel");
  TestCorrectness(MFFT::Transform, "MFFT::Transform");

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
