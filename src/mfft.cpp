#include "mfft.h"

namespace MFFT {

dims GetIndices(int index, const dims& N) {
  dims indices;
  for (size_t i = 0; i < N.size(); ++i) {
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
    arr.next[ri] = arr.X[i];
  }
}

void ComputeW1(Array& arr, int start, int end, int inc) {
  cmplx w(1.);
  for (int k = start; k < end; k += inc) {
    arr.ws[k] = w;
    w *= arr.wln;
  }
}

void ComputeW2(Array& arr, int start, int end, int inc) {
  int block_size = (arr.n / arr.num_threads);
  for (int k = start; k < end; k += inc) {
    arr.ws[k] *= arr.wp[k / block_size];
  }
}

void ComputeWL(Array& arr) {
  Parallelize(ComputeW1, arr, 1);

  cmplx tmp = arr.ws[arr.n / arr.num_threads - 1] * arr.wln;
  arr.wp[0] = cmplx(1.);

  for (int i = 1; i < arr.num_threads; ++i) {
    arr.wp[i] = arr.wp[i - 1] * tmp;
  }

  Parallelize(ComputeW2, arr, 1);
}

void AssignTwiddle(Array& arr, int start, int end, int Nj) {
  int inv = (arr.n / arr.Lprev);
  for (int i = start; i < end; i += Nj) {
    arr.twiddle[i] = arr.ws[i / inv];
  }
}

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

void DivN(Array& arr, int start, int end, int inc) {
  for (int i = start; i < end; i += inc) {
    arr.X[i] /= double(arr.n);
  }
}

void Transform(jarray& X_orig, bool invert, size_t num_threads) {
  Array arr(X_orig, num_threads, invert, true, false);

  const dims& N = arr.N;
  jarray& X = arr.X;

  int& n = arr.n;
  int& Lj = arr.Lj;
  int& Lprev = arr.Lprev;
  int& Nj = arr.Nj;

  Lj = 1;

  Parallelize(ReverseIndex, arr, 1);
  arr.update();

  for (int jn = N.size() - 1; jn >= 0; --jn) {
    Nj = N[jn];
    Lprev = Lj;
    Lj *= Nj;

    // twiddle
    arr.wln = GetW(arr.Lj, arr.invert);
    ComputeWL(arr);
    Parallelize(AssignTwiddle, arr, Nj);

    // fourier
    Parallelize(SeqTransform, arr, Nj);

    // transpose
    Parallelize2(Transpose, arr, Lj, n / Lj);
    arr.update();
  }

  if (invert) {
    Parallelize(DivN, arr, 1);
  }

  std::copy(std::execution::par_unseq, X.begin(), X.end(), X_orig.begin());
}

}  // namespace MFFT
