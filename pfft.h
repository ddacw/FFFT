#pragma once

#include <cmath>
#include "main.h"
#include "refft.h"

class PFFT {
 public:
  PFFT(const jarray& data, size_t p)
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

    // set p to be power of 2, optional
    int p_bit = CeilBit(p);
    if (this->p < (1 << p_bit)) {
      this->p = (1 << (p_bit - 1));
    }
    InitN();
  }

  void Transform();
  void Compress();
  void IndexReverse();
  void Expand();
  void InitN();
  void Print();

 private:
  size_t n;  // size of data, power of 2
  size_t p;  // number of threads, also power of 2 for simplicity
  std::vector<size_t> N;
  jarray data;
  jarray transformed;
};

void PFFT::InitN() {  // TODO: redo
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

void PFFT::Print() {
  std::cerr << "data:";
  for (size_t i = 0; i < n; ++i) {
    if (i % 10 == 0) {
      std::cerr << std::endl;
    }
    std::cerr << data[i] << ' ';
  }
  std::cerr << std::endl
            << "N:";
  for (auto x : N) {
    std::cerr << x << ' ';
  }
}
