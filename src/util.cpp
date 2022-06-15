#include "util.h"

Array::Array(jarray X, int num_threads, bool invert, bool init_dim = false)
    : X(X), num_threads(num_threads), invert(invert) {
  n_bit = PadZero(this->X);
  n = (1 << n_bit);
  twiddle = jarray(X.size());
  if (init_dim) {
    InitN(this->N, X.size(), num_threads);
  }
  next = jarray(n);
}

void Array::update() {
  std::copy(std::execution::par_unseq, next.begin(), next.end(), X.begin());
}

void Parallelize(std::function<void(Array&, int, int)> f, Array& arr) {
  size_t num_threads = arr.num_threads;
  size_t block_size = arr.X.size() / num_threads;
  std::vector<std::thread> workers(num_threads - 1);
  size_t start_block = 0;
  for (size_t i = 0; i < num_threads - 1; ++i) {
    size_t end_block = start_block + block_size;
    workers[i] = std::thread(f, std::ref(arr), start_block, end_block);
    start_block = end_block;
  }
  f(arr, start_block, arr.X.size());
  for (size_t i = 0; i < num_threads - 1; ++i) {
    workers[i].join();
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
