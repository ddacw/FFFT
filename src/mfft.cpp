#include "mfft.h"

namespace MFFT {

typedef std::vector<size_t> dims;

void Fourier(jarray& x, const jarray& twiddle, int start, int N) {
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

bool InitN(dims& N, size_t n, size_t num_threads) {  // TODO: redo
  int root = std::round(std::pow(double(n), 1. / double(num_threads)));
  size_t cur = 1;
  while (cur * root <= n) {
    N.push_back(root);
    cur *= root;
  }
  if (cur < n) {
    N.push_back(n / cur);
  }

  size_t n_check = 1;
  for (size_t x : N) {
    n_check *= x;
  }
  return n == n_check;
}

dims GetIndices(int index, const dims& N) {
  dims indices;
  for (int i = N.size() - 1; i >= 0; --i) {
    indices.push_back(index % N[i]);
    index /= N[i];
  }
  return indices;
}

void ReverseIndex(jarray& X, jarray& next, const dims& N) {
  for (size_t i = 0; i < X.size(); ++i) {
    dims indices = GetIndices(i, N);
    int ri = 0;
    for (size_t j = 0; j < N.size(); ++j) {
      ri *= N[j];
      ri += indices[j];
    }
    next[ri] = X[i];
  }
}

void Transform(jarray& X, size_t num_threads) {
  int n_bit = PadZero(X);
  int n = (1 << n_bit);
  dims N;
  if (!InitN(N, n, num_threads)) {
    std::cerr << "Invalid sizes." << std::endl;
    return;
  }

  jarray next(n);
  ReverseIndex(X, next, N);
  X = next;

  int Lj = 1;
  int Lprev;

  for (size_t jn = 0; jn < N.size(); ++jn) {
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
    for (int k = 0; k < Lj; ++k) {
      for (int i = 0; i < n / Lj; ++i) {
        // int Lprev = Lj / Nj;
        int kj = k / Lprev;
        int k1 = k % Lprev;
        next[k * (n / Lj) + i] = X[k1 * (n / Lprev) + i * Nj + kj];
      }
    }
    X = next;
  }
}

}  // namespace MFFT
