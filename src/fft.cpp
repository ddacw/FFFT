#include "fft.h"

namespace FFT {

// Peforms the DFT by defintion in O(N^2).
void DFT(jarray& X, bool invert) {
  int N = X.size();
  jarray next(N, 0.);
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < N; ++i) {
      next[k] += X[i] * GetW(N, invert, k * i);
    }
  }
  X = next;
  if (invert) {
    for (int i = 0; i < N; ++i) {
      X[i] /= N;
    }
  }
}

// Implements the Radix-2 FFT algorihm.
void FFTRecAux(jarray& X, bool invert) {
  int n = X.size();
  if (n == 1) {
    return;
  }
  jarray odd(n / 2), even(n / 2);
  for (int i = 0, j = 0; i < n; i += 2, ++j) {
    even[j] = X[i];
    odd[j] = X[i + 1];
  }
  FFTRecAux(even, invert);
  FFTRecAux(odd, invert);

  cmplx wn = GetW(n, invert, 1);
  cmplx w(1.);

  int half = n / 2;

  for (int i = 0; i < half; ++i) {
    X[i] = even[i] + w * odd[i];
    X[i + half] = even[i] - w * odd[i];
    w *= wn;
  }
}

// Wrapper for the Radix-2 algorithm.
void FFTRec(jarray& X, bool invert) {
  PadZero(X);
  FFTRecAux(X, invert);
  if (invert) {
    for (size_t i = 0; i < X.size(); ++i) {
      X[i] /= X.size();
    }
  }
}

// Performs FFTRec without recursion.
void FFTSeq(jarray& X, bool invert) {
  if (X.size() == 1) {
    return;
  }

  // padding
  int n_bit = PadZero(X);
  int n = (1 << n_bit);

  // reverse bit
  for (int i = 0; i < n; ++i) {
    int j = RevBit(i, n_bit);
    if (i < j) {
      std::swap(X[i], X[j]);
    }
  }

  jarray next(n);
  for (int step = 1; step < n; step <<= 1) {
    cmplx w(1.);
    cmplx wn = GetW(2 * step, invert, 1);
    jarray W;
    for (int i = 0; i < step; ++i) {
      W.push_back(w);
      w *= wn;
    }

    int start_even = 0;
    int start_odd = start_even + step;
    while (start_even < n) {
      for (int i = 0; i < step; ++i) {
        next[start_even + i] = X[start_even + i] + W[i] * X[start_odd + i];
        next[start_odd + i] = X[start_even + i] - W[i] * X[start_odd + i];
      }
      start_even += 2 * step;
      start_odd = start_even + step;
    }
    for (int i = 0; i < n; ++i) {
      X[i] = next[i];
    }
  }

  if (invert) {
    for (int i = 0; i < n; ++i) {
      X[i] /= n;
    }
  }
}

// FFTParallel aux function.
// Swap X[i] with X[j], given that i and j have 'reversed' binary
// representations.
void RevBitThread(Array& arr, size_t start, size_t end, int inc) {
  for (size_t i = start; i < end; ++i) {
    int j = RevBit(int(i), arr.n_bit);
    // thread-safe, since each pair is swapped only once
    if (int(i) < j) {
      std::swap(arr.X[i], arr.X[j]);
    }
  }
}

// FFTParallel aux function.
// Computes the current DFT based on the results of the last step.
void ComputeThread(Array& arr, int begin, int end, int inc) {
  int& step = arr.step;
  int start_even = begin;
  int start_odd = start_even + arr.step;
  int step2_per_thread = (end - begin) / step;
  jarray& X = arr.X;
  jarray& next = arr.next;
  jarray& W = arr.ws;
  for (int iter = 0; iter < step2_per_thread; ++iter) {
    for (int i = 0; i < arr.step; i += inc) {
      next[start_even + i] = X[start_even + i] + W[i] * X[start_odd + i];
      next[start_odd + i] = X[start_even + i] - W[i] * X[start_odd + i];
    }
    start_even += step * 2;
    start_odd = start_even + step;
  }
}

// Performs FFTSeq in parallel.
void FFTParallel(jarray& X_orig, bool invert, size_t num_threads) {
  if (X_orig.size() == 1) {
    return;
  }
  Array arr(X_orig, num_threads, invert, false, true);

  jarray& X = arr.X;
  int& n = arr.n;

  // Rearrange.
  Parallelize(RevBitThread, arr, 1);

  int& step = arr.step;
  // Log(N) steps.
  for (step = 1; step < arr.n; step <<= 1) {
    
    // Precomputed roots because it might be faster.
    arr.ws = jarray(step);
    for (int i = 0; i < step; ++i) {
      arr.ws[i] = GetW(2 * step, invert, i);
    }

    // Reduce the number of theads used for the last steps.
    while (n / (2 * step * num_threads) == 0) {
      num_threads /= 2;
    }

    Parallelize(ComputeThread, arr, 1);
    arr.update();
  }

  // Epilogue.
  if (invert) {
    std::transform(std::execution::par_unseq, X.begin(), X.end(), X.begin(),
                   [=](cmplx x) -> cmplx { return x / double(n); });
  }
  std::copy(std::execution::par_unseq, X.begin(), X.end(), X_orig.begin());
}

}  // namespace FFT
