#pragma once

#include "main.h"

// 1-thread FFT

int CeilBit(int n) {
  int pow2 = 1, count = 0;
  while (pow2 < n) {
    pow2 <<= 1;
    count++;
  }
  return count;
}

int GetBit(int mask, int n) {
  return (mask >> n) & 1;
}

int RevBit(int mask, int count) {
  for (int i = 0, j = count - 1; i < j; ++i, --j) {
    if (GetBit(mask, i) != GetBit(mask, j)) {
      mask ^= (1 << i);
      mask ^= (1 << j);
    }
  }
  return mask;
}

void RefFFT(jarray& X, bool invert) {
  if (X.size() == 1) {
    return;
  }

  int n_bit = CeilBit(X.size());
  int n = (1 << n_bit);
  while (X.size() < size_t(n)) {
    X.emplace_back(0);
  }

  for (int i = 0; i < n; ++i) {
    int j = RevBit(i, n_bit);
    if (i < j) {
      std::swap(X[i], X[j]);
    }
  }

  jarray next(n);
  for (int step = 1; step < n; step <<= 1) {
    double angle = M_PI / step;
    if (invert) {
      angle = -angle;
    }

    cmplx w(1);
    cmplx wn = exp(cmplx(0, -1) * angle);
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

// Returns the shortened & transformed data.
jarray Compress(const jarray& original, size_t new_size) {
  jarray compressed(original);
  RefFFT(compressed, 0);
  compressed.resize(new_size);
  return compressed;
}

jarray Uncompress(const jarray& compressed, size_t original_size) {
  jarray recovered(original_size, cmplx(0, 0));
  for (size_t i = 0; i < compressed.size(); ++i) {
    recovered[i] = compressed[i];
  }
  RefFFT(recovered, 1);
  return recovered;
}

// Computes the mean square error.
double MSE(const jarray& original, const jarray& recovered) {
  size_t n = original.size();
  if (recovered.size() != n) {
    return -1;
  }
  cmplx mse = 0.;
  for (size_t i = 0; i < n; ++i) {
    cmplx diff = original[i] - recovered[i];
    mse += (diff * diff);
  }
  return mse.real() / double(n);
}
