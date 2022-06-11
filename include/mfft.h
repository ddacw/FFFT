#pragma once

#include <cassert>
#include <cmath>

#include "fft.h"
#include "util.h"

class MFFT {
 public:
  MFFT(const jarray& data, size_t p);

  void Transform();
  void Compress();
  void ReverseIndexAux(std::vector<size_t>& indices);
  void ReverseIndex();
  void Expand();
  void InitN();
  void Print();

  size_t GetIndex(const std::vector<size_t>& indices);

  size_t GetRevIndex(const std::vector<size_t>& indices);

  static cmplx GetW(int N) {
    double angle = 2. * M_PI / double(N);
    // return exp(cmplx(0, -1) * angle);
    return cmplx(cos(angle), sin(angle));
  }

  static void Fourier(jarray& x, const jarray& twiddle, int start, int N) {
    jarray next(N, 0.);
    cmplx wn = GetW(N);
    cmplx w(1.);
    cmplx w_twd = twiddle[start];
    for (int k = 0; k < N; ++k) {
      cmplx w_tmp(1.);
      for (int i = 0; i < N; ++i) {
        next[k] += x[start + i] * w_tmp;
        w_tmp *= w;
        w_tmp *= w_twd;
      }
      w *= wn;
    }
    for (int i = 0; i < N; ++i) {
      x[start + i] = next[i];
    }
  }

  size_t n;  // size of data, power of 2
  size_t p;  // number of threads such that the p-th root of n is an integer
  std::vector<size_t> N;
  jarray data;
  jarray X;

 private:
};

void MFFT::Transform() {
  int Lj = 1;
  int Lprev;
  cmplx wt(1.);

  for (int jn = 0; jn < N.size(); ++jn) {
    int Nj = N[jn];
    Lprev = Lj;
    Lj *= Nj;

    // twiddle
    jarray twiddle(n);
    cmplx wn = GetW(Lj);
    cmplx w(1.);

    for (int k = 0; k < Lprev; ++k) {
      for (int i = 0; i < n / Lprev; i += Nj) {
        twiddle[k * (n / Lprev) + i] = w;
      }
      w *= wn;
    }

    // transform
    for (int start = 0; start < n; start += Nj) {
      Fourier(X, twiddle, start, Nj);
    }

    // shift
    jarray next(n);
    for (int k = 0; k < Lj; ++k) {
      for (int i = 0; i < n / Lj; ++i) {
        // int Lprev = Lj / Nj;
        int kj = k / Lprev;
        int k1 = k % Lprev;
        next[k * (n / Lj) + i] = X[k1 * (n / Lprev) + i * Nj + kj];
      }
    }
    X = next;

    Print();
  }
}

MFFT::MFFT(const jarray& data, size_t p) : p(p) {
  for (cmplx x : data) {
    this->data.push_back(x);
  }
  // padding the array
  int n_bit = CeilBit(this->data.size());
  n = (1 << n_bit);
  while (this->data.size() < n) {
    this->data.emplace_back(0);
  }

  InitN();
  ReverseIndex();
}

void MFFT::InitN() {  // TODO: redo
  int root = std::round(std::pow(double(n), 1. / double(p)));
  size_t cur = 1;
  while (cur * root <= n) {
    N.push_back(root);
    cur *= root;
  }
  if (cur < n) {
    N.push_back(n / cur);
  }
}

size_t MFFT::GetIndex(const std::vector<size_t>& indices) {
  // assert(indices.size() == p);
  size_t index = 0;
  for (size_t i = 0; i < p; ++i) {
    index = index * N[i] + indices[i];
  }
  return index;
}

size_t MFFT::GetRevIndex(const std::vector<size_t>& indices) {
  assert(indices.size() == p);
  size_t index = 0;
  for (int i = p - 1; i >= 0; --i) {
    index = index * N[i] + indices[i];
  }
  return index;
}

void MFFT::ReverseIndexAux(std::vector<size_t>& indices) {
  if (indices.size() == p) {
    X[GetIndex(indices)] = data[GetRevIndex(indices)];
    return;
  }
  for (size_t i = 0; i < N[indices.size()]; ++i) {
    indices.push_back(i);
    ReverseIndexAux(indices);
    indices.pop_back();
  }
}

// x0[np,...,n1] =  x'[n1,...,np]
// Generates index-reversed X;
void MFFT::ReverseIndex() {
  X = jarray(data);
  std::vector<size_t> indices;
  ReverseIndexAux(indices);
}

void MFFT::Print() {
  // std::cerr << "data:";
  // for (size_t i = 0; i < n; ++i) {
  //   if (i % 10 == 0) {
  //     std::cerr << std::endl;
  //   }
  //   std::cerr << data[i] << ' ';
  // }
  std::cerr << std::endl << "X:";
  for (size_t i = 0; i < n; ++i) {
    if (i % 8 == 0) {
      std::cerr << std::endl;
    }
    std::cerr << X[i] << ' ';
  }
  std::cerr << std::endl << "N:";
  for (auto x : N) {
    std::cerr << x << ' ';
  }
}
