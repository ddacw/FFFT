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

dims GetIndices(int index, const dims& N) {
  dims indices;
  for (int i = N.size() - 1; i >= 0; --i) {
    indices.push_back(index % N[i]);
    index /= N[i];
  }
  return indices;
}

void ReverseIndex(Array& arr, int start, int end, int inc) {
  for (int i = start; i < end; i += inc) {
    dims indices = GetIndices(i, arr.N);
    int ri = 0;
    for (size_t j = 0; j < arr.N.size(); ++j) {
      ri *= arr.N[j];
      ri += indices[j];
    }
    if (i < ri) {
      std::swap(arr.X[i], arr.X[ri]);
    }
  }
}

void SeqTransform(Array& arr, int start, int end, int Nj) {
  for (int i = start; i < end; i += Nj) {
    Fourier(arr.X, arr.invert, arr.twiddle, i, Nj);
  }
}

void Transpose(Array& arr, int startk, int endk, int starti, int endi) {
  for (int k = startk; k < endk; ++k) {
    for (int i = starti; i < endi; ++i) {
      int kj = k / arr.Lprev;
      int k1 = k % arr.Lprev;
      arr.next[k * (arr.n / arr.Lj) + i] =
          arr.X[k1 * (arr.n / arr.Lprev) + i * arr.Nj + kj];
    }
  }
}

void Transform(jarray& X_orig, bool invert, size_t num_threads) {
  Array arr(X_orig, num_threads, invert, true);

  const dims& N = arr.N;
  jarray& X = arr.X;
  jarray& next = arr.next;
  jarray& twiddle = arr.twiddle;
  int& n = arr.n;
  int& Lj = arr.Lj;
  int& Lprev = arr.Lprev;
  int& Nj = arr.Nj;
  Lj = 1;

  Parallelize(ReverseIndex, arr, 1);

  for (size_t jn = 0; jn < N.size(); ++jn) {
    Nj = N[jn];
    Lprev = Lj;
    Lj *= Nj;

    // twiddle
    cmplx wn = GetW(Lj, invert);
    cmplx w(1.);

    for (int k = 0; k < Lprev; ++k) {
      for (int i = 0; i < n / Lprev; i += Nj) {
        twiddle[k * (n / Lprev) + i] = w;
      }
      w *= wn;
    }

    // fourier
    Parallelize(SeqTransform, arr, Nj);

    // transpose
    Parallelize2(Transpose, arr, Lj, n / Lj);
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
