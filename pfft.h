#pragma once

#include <cmath>
#include "cassert"
#include "main.h"
#include "refft.h"

class GeneralizedFFT {
 public:
  GeneralizedFFT(const jarray& data, size_t p)
      : p(p) {
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
    IndexReverse();
  }

  void Transform();
  void Compress();
  void IndexReverseAux(std::vector<size_t>& indices);
  void IndexReverse();
  void Expand();
  void InitN();
  void Print();

  size_t GetIndex(const std::vector<size_t>& indices) {
    assert(indices.size() == p);
    size_t index = 0;
    for (size_t i = 0; i < p; ++i) {
      index = index * N[i] + indices[i];
    }
    return index;
  }

  size_t GetRevIndex(const std::vector<size_t>& indices) {
    assert(indices.size() == p);
    size_t index = 0;
    for (int i = p - 1; i >= 0; --i) {
      index = index * N[i] + indices[i];
    }
    return index;
  }

  static cmplx GetW(int N) {
    double angle = M_PI / double(N);
    return exp(cmplx(0, -1) * angle);
  }

 private:
  size_t n;  // size of data, power of 2
  size_t p;  // number of threads such that the p-th root of n is an integer
  std::vector<size_t> N;
  jarray data;
  jarray data0;
  jarray X;
  jarray Xp;
};

void GeneralizedFFT::InitN() {  // TODO: redo
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

void GeneralizedFFT::IndexReverseAux(std::vector<size_t>& indices) {
  if (indices.size() == p) {
    data0[GetIndex(indices)] = data[GetRevIndex(indices)];
    return;
  }
  for (size_t i = 0; i < N[indices.size()]; ++i) {
    indices.push_back(i);
    IndexReverseAux(indices);
    indices.pop_back();
  }
}

// x0[np,...,n1] =  x'[n1,...,np]
// Generates index-reversed data0;
void GeneralizedFFT::IndexReverse() {
  data0 = jarray(data);
  std::vector<size_t> indices;
  IndexReverseAux(indices);
}

void GeneralizedFFT::Print() {
  std::cerr << "data:";
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 0) {
      std::cerr << std::endl;
    }
    std::cerr << data[i] << ' ';
  }
  std::cerr << std::endl
            << "data0:";
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 0) {
      std::cerr << std::endl;
    }
    std::cerr << data0[i] << ' ';
  }
  std::cerr << std::endl
            << "N:";
  for (auto x : N) {
    std::cerr << x << ' ';
  }
}
