#include "fft.h"

namespace FFT {

// ----------------- RECURSIVE ---------------------- //

void FFTRecAux(jarray& X) {
  int n = X.size();
  if (n == 1) {
    return;
  }
  jarray odd(n / 2), even(n / 2);
  for (int i = 0, j = 0; i < n; i += 2, ++j) {
    even[j] = X[i];
    odd[j] = X[i + 1];
  }
  FFTRecAux(even);
  FFTRecAux(odd);

  cmplx wn = GetW(n);
  cmplx w(1.);

  int half = n / 2;

  for (int i = 0; i < half; ++i) {
    X[i] = even[i] + w * odd[i];
    X[i + half] = even[i] - w * odd[i];
    w *= wn;
  }
}

void FFTRec(jarray& X) {
  PadZero(X);
  FFTRecAux(X);
}

// ----------------- ITERATIVE ------------------ //
void FFTSeq(jarray& X, bool invert) {
  if (X.size() == 1) {
    return;
  }

  int n_bit = PadZero(X);
  int n = (1 << n_bit);

  for (int i = 0; i < n; ++i) {
    int j = RevBit(i, n_bit);
    if (i < j) {
      std::swap(X[i], X[j]);
    }
  }

  jarray next(n);
  for (int step = 1; step < n; step <<= 1) {
    double angle = M_PI / step;  // 2*PI / (2*step)
    if (invert) {
      angle = -angle;
    }

    cmplx w(1.), wn(cos(angle), sin(angle));

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

// -------------------- PARALLEL -------------------- //
void RevBitThread(jarray& X, size_t start, size_t end, int n_bit) {
  for (size_t i = start; i < end; ++i) {
    int j = RevBit(int(i), n_bit);
    // thread-safe, since each pair is swapped only once
    if (int(i) < j) {
      std::swap(X[i], X[j]);
    }
  }
}

void RevBitParallel(jarray& X, size_t num_threads, int n_bit) {
  size_t block_size = X.size() / num_threads;
  std::vector<std::thread> workers(num_threads - 1);
  size_t start_block = 0;
  for (size_t i = 0; i < num_threads - 1; ++i) {
    size_t end_block = start_block + block_size;
    workers[i] =
        std::thread(RevBitThread, std::ref(X), start_block, end_block, n_bit);
    start_block = end_block;
  }
  RevBitThread(X, start_block, X.size(), n_bit);
  for (size_t i = 0; i < num_threads - 1; ++i) {
    workers[i].join();
  }
}

void TwiddleThread(const jarray& X, const jarray& W, jarray& next, int begin,
                   int step2_per_thread, int step) {
  std::cerr << begin << ' ';
  int start_even = begin;
  int start_odd = start_even + step;
  for (int iter = 0; iter < step2_per_thread; ++iter) {
    for (int i = 0; i < step; ++i) {
      next[start_even + i] = X[start_even + i] + W[i] * X[start_odd + i];
      next[start_odd + i] = X[start_even + i] - W[i] * X[start_odd + i];
    }
    start_even += step * 2;
    start_odd = start_even + step;
  }
}

void TwiddleParallel(jarray& X, const jarray& W, int step, size_t num_threads) {
  int n = X.size();
  // (step | step) | (step | step )
  // (even | odd ) | (even | odd  )
  int step2_per_thread;
  while ((step2_per_thread = n / (2 * step * num_threads)) == 0) {
    num_threads /= 2;
  }

  std::cerr << num_threads << std::endl;
  jarray next(n);
  int begin = 0;
  // int start_odd = start_even + step;

  std::vector<std::thread> workers(num_threads - 1);
  for (size_t i = 0; i < num_threads - 1; ++i) {
    workers[i] = std::thread(TwiddleThread, std::cref(X), std::cref(W),
                             std::ref(next), begin, step2_per_thread, step);
    begin += step2_per_thread * 2 * step;
  }
  TwiddleThread(X, W, next, begin, step2_per_thread, step);
  for (size_t i = 0; i < num_threads - 1; ++i) {
    workers[i].join();
  }
  std::copy(std::execution::par_unseq, next.begin(), next.end(), X.begin());
  // for (int i = 0; i < n; ++i) {
  //   X[i] = next[i];
  // }
}

void FFTParallel(jarray& X, bool invert, size_t num_threads) {
  if (X.size() == 1) {
    return;
  }

  // Padding
  int n_bit = PadZero(X);
  int n = (1 << n_bit);

  // Rearrange
  RevBitParallel(X, num_threads, n_bit);

  jarray next(n);
  for (int step = 1; step < n; step <<= 1) {
    double angle = M_PI / step;
    if (invert) {
      angle = -angle;
    }

    cmplx w(1);
    cmplx wn = exp(cmplx(0, 1) * angle);
    jarray W;
    for (int i = 0; i < step; ++i) {
      W.push_back(w);
      w *= wn;
    }

    TwiddleParallel(X, W, step, num_threads);
  }

  if (invert) {
    std::transform(std::execution::par_unseq, X.begin(), X.end(), X.begin(),
                   [=](cmplx x) -> cmplx { return x / double(n); });
    // for (int i = 0; i < n; ++i) {
    //   X[i] /= n;
    // }
  }
}

// -------------------- APPLICATION -------------------- //

// Returns the shortened & transformed data.
jarray Compress(const jarray& original, size_t new_size) {
  jarray compressed(original);
  FFTSeq(compressed, 0);
  compressed.resize(new_size);
  return compressed;
}

jarray Uncompress(const jarray& compressed, size_t original_size) {
  jarray recovered(original_size, cmplx(0, 0));
  for (size_t i = 0; i < compressed.size(); ++i) {
    recovered[i] = compressed[i];
  }
  FFTSeq(recovered, 1);
  return recovered;
}
}  // namespace FFT
