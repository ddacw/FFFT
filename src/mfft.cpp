#include "mfft.h"

namespace MFFT {

void Fourier(jarray& x, bool invert, const jarray& twiddle, int start, int N) {
  jarray next(N, 0.);
  cmplx wn = GetW(N, invert);
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

void ReverseIndex(Array& arr, int start, int end) {
  for (size_t i = start; i < end; ++i) {
    dims indices = GetIndices(i, arr.N);
    int ri = 0;
    for (size_t j = 0; j < arr.N.size(); ++j) {
      ri *= arr.N[j];
      ri += indices[j];
    }
    arr.next[ri] = arr.X[i];
  }
}

void Transform(jarray& X_orig, bool invert, size_t num_threads) {
  Array arr(X_orig, num_threads, invert, true);

  const dims& N = arr.N;
  jarray& X = arr.X;
  jarray& next = arr.next;
  int n = arr.n;
  int Lj = 1;
  int Lprev;

  Parallelize(ReverseIndex, arr);
  arr.update();

  for (size_t jn = 0; jn < N.size(); ++jn) {
    int Nj = N[jn];
    Lprev = Lj;
    Lj *= Nj;

    // twiddle
    jarray twiddle(n);
    cmplx wn = GetW(Lj, invert);
    cmplx w(1.);

    for (int k = 0; k < Lprev; ++k) {
      for (int i = 0; i < n / Lprev; i += Nj) {
        twiddle[k * (n / Lprev) + i] = w;
      }
      w *= wn;
    }

    // transform
    for (int start = 0; start < n; start += Nj) {
      Fourier(X, invert, twiddle, start, Nj);
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
    arr.update();
  }
  if (invert) {
    for (int i = 0; i < n; ++i) {
      X[i] /= n;
    }
  }
  X_orig = X;
}

}  // namespace MFFT
