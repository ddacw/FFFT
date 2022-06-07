#include "main.h"
#include "refft.h"
#include "pfft.h"

int main(int argc, char* argv[]) {
  jarray X;
  size_t n_threads = 2;
  for (int i = 0; i < 64; ++i) {
    X.emplace_back(double(i % 4));
  }
  PFFT pfft(X, n_threads);
  pfft.Print();
}
