#include "util.h"

Array::Array(jarray X, int num_threads, bool invert, bool init_dim = false,
             bool pad_zero = true)
    : X(X), num_threads(num_threads), invert(invert) {
  if (pad_zero) {
    n_bit = PadZero(this->X);
  }
  n = this->X.size();
  twiddle = jarray(this->X.size());
  if (init_dim) {
    if (!InitN(this->N, this->X.size(), num_threads)) {
      throw std::invalid_argument("Unable to factorize n by num_threads.");
    }
  }
  next = jarray(n);
  ws = jarray(n);
  wp = jarray(n);
}

void Array::update() {
  std::copy(std::execution::par_unseq, next.begin(), next.end(), X.begin());
}

void Parallelize(std::function<void(Array&, int, int, int)> f, Array& arr,
                 int inc) {
  size_t num_threads = arr.num_threads;
  size_t block_size = (arr.n / inc) / num_threads;
  std::vector<std::thread> workers(num_threads - 1);
  size_t start_block = 0;

  if (block_size) {
    for (size_t i = 0; i < num_threads - 1; ++i) {
      size_t end_block = start_block + block_size * inc;
      workers[i] = std::thread(f, std::ref(arr), start_block, end_block, inc);
      start_block = end_block;
    }
  }
  f(arr, start_block, arr.n, inc);
  if (block_size) {
    for (size_t i = 0; i < num_threads - 1; ++i) {
      workers[i].join();
    }
  }
}

void Parallelize2(std::function<void(Array&, int, int, int, int)> f, Array& arr,
                  int N, int M) {
  size_t num_threads = arr.num_threads;
  size_t start_block = 0;
  std::vector<std::thread> workers(num_threads - 1);

  if (N > M) {
    size_t block_size = N / num_threads;
    for (size_t i = 0; i < num_threads - 1; ++i) {
      size_t end_block = start_block + block_size;
      workers[i] = std::thread(f, std::ref(arr), start_block, end_block, 0, M);
      start_block = end_block;
    }
    f(arr, start_block, N, 0, M);
  } else {
    size_t block_size = M / num_threads;
    for (size_t i = 0; i < num_threads - 1; ++i) {
      size_t end_block = start_block + block_size;
      workers[i] = std::thread(f, std::ref(arr), 0, N, start_block, end_block);
      start_block = end_block;
    }
    f(arr, 0, N, start_block, M);
  }
  for (size_t i = 0; i < num_threads - 1; ++i) {
    workers[i].join();
  }
}

bool InitN(dims& N, size_t n, size_t num_threads) {
  size_t tmp = n;
  while (tmp / num_threads) {
    N.push_back(num_threads);
    tmp /= num_threads;
  }
  if (tmp == 2) {
    N.push_back(2);
  }

  size_t n_check = 1;
  for (size_t x : N) {
    n_check *= x;
  }
  return n == n_check;
}

cmplx GetW(int N, bool invert) {
  double angle = 2. * M_PI / double(N);
  if (invert) {
    angle = -angle;
  }
  return exp(cmplx(0, -1) * angle);
  // return cmplx(cos(angle), sin(angle));
}

int CeilBit(int n) {
  int pow2 = 1, count = 0;
  while (pow2 < n) {
    pow2 <<= 1;
    count++;
  }
  return count;
}

int GetBit(int mask, int n) { return (mask >> n) & 1; }

int RevBit(int mask, int count) {
  for (int i = 0, j = count - 1; i < j; ++i, --j) {
    if (GetBit(mask, i) != GetBit(mask, j)) {
      mask ^= (1 << i);
      mask ^= (1 << j);
    }
  }
  return mask;
}

int PadZero(jarray& X) {
  int n_bit = CeilBit(X.size());
  int n = (1 << n_bit);
  while (X.size() < size_t(n)) {
    X.emplace_back(0);
  }
  return n_bit;
}
