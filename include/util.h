#pragma once

#include <complex>
#include <exception>
#include <execution>
#include <functional>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>

typedef std::complex<double> cmplx;
typedef std::vector<size_t> dims;
typedef std::vector<cmplx> jarray;
const size_t DATA_SIZE = 1 << 18;

class Array {
 public:
  Array(jarray X, int num_threads, bool invert, bool init_dim, bool pad_zero);

  void update();

  bool invert;
  jarray X, next;
  jarray twiddle, ws; // roots of unity
  dims N;
  int n_bit, num_threads, n, step;
  int Lj, Lprev, Nj;
};

void Parallelize(std::function<void(Array&, int, int, int)> f, Array& arr,
                 int inc);

void Parallelize2(std::function<void(Array&, int, int, int, int)> f, Array& arr,
                  int N, int M);

bool InitN(dims& N, size_t n, size_t num_threads);

cmplx GetW(int N, bool invert, double k);

int CeilBit(int n);

int GetBit(int mask, int n);

int RevBit(int mask, int count);

int PadZero(jarray& X);
