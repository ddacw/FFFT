#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>

static std::uniform_int_distribution<> uniform(-1, RAND_MAX);

#include "fft.h"
#include "mfft.h"

const double eps_small = 1e-8;

typedef std::function<void(jarray&, bool, size_t)> func_fft;

jarray Generate(int N, int seed, double scale, double shift) {
  auto engine = std::default_random_engine();
  jarray X(N);
  for (int i = 0; i < N; ++i) {
    X[i] = cos(double(i)) + cos(uniform(engine)) + sin(i + 1.) +
           sin(uniform(engine) + 1);
    X[i] *= scale;
    X[i] += shift;
  }
  return X;
}

bool CorrectnessAux(func_fft fft, int N, int num_threads) {
  jarray X_orig = Generate(N, 0, 30.5, 3.05);
  jarray X(X_orig);

  // FFT.
  fft(X, false, num_threads);
  // IFFT.
  fft(X, true, num_threads);

  double error = 0;
  for (int i = 0; i < N; ++i) {
    error += abs(X[i] - X_orig[i]);
  }
  return error < eps_small;
}

void TestCorrectness(func_fft fft, const std::string& fft_name) {
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
  std::cerr << "5 ";
  verdict &= CorrectnessAux(fft, 1, 8);
  std::cerr << "6 ";
  verdict &= CorrectnessAux(fft, 7, 1);

  if (!verdict) {
    std::cerr << fft_name << " is incorrect !!!" << std::endl;
  } else {
    std::cerr << fft_name << " passed." << std::endl;
  }
}

int Benchmark(func_fft fft, int N, size_t num_threads) {
  jarray X = Generate(N, 0, 10, 2.3);
  auto start = std::chrono::steady_clock::now();
  fft(X, false, num_threads);
  auto finish = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::microseconds>(finish - start)
          .count();
  return elapsed;
}

const int Ns[8] = {5000,   10000,   50000,    100000,
                   500000, 1000000, 10000000, 50000000};

void TestExecutionTime(func_fft fft, std::string fft_name,
                       size_t max_num_threads) {
  std::cerr.imbue(std::locale(""));
  std::cerr << "\nBenchmark " << fft_name << ":" << std::endl;
  for (size_t num_threads = 1; num_threads <= max_num_threads;
       num_threads <<= 1) {
    std::cerr << std::setw(2) << num_threads << " thread(s): ";
    for (int N : Ns) {
      int elapsed = Benchmark(fft, N, num_threads);
      std::cerr << '|' << elapsed;
    }
    std::cerr << '|' << std::endl;
  }
}

int main(int argc, char* argv[]) {
  // Correctness.
  TestCorrectness(FFT::DFT, "FFT::DFT");
  TestCorrectness(FFT::FFTRec, "FFT::FFTRec");
  TestCorrectness(FFT::FFTSeq, "FFT::FFTSeq");
  TestCorrectness(FFT::FFTParallel, "FFT::FFTParallel");
  TestCorrectness(MFFT::Transform, "MFFT::Transform");

  TestExecutionTime(FFT::FFTRec, "FFT::FFTRec", 1);
  TestExecutionTime(FFT::FFTSeq, "FFT::FFTSeq", 1);
  TestExecutionTime(FFT::FFTParallel, "FFT::FFTParallel", 32);
  TestExecutionTime(MFFT::Transform, "MFFT::Transform", 32);

  // // Optional: test with the weather dataset.
  // std::ifstream infile;
  // infile.open("data/temp.txt");
  // jarray X_orig;
  // double x;
  // while (infile >> x) {
  //   X_orig.push_back(x);
  // }
}
