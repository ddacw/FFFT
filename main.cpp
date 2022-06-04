#include "main.h"
#include "refft.h"

int main(int argc, char* argv[]) {
  jarray X;

  for (int i = 0; i < 128; ++i) {
    X.emplace_back(double(i % 4));
  }
  for (int csize = 1; csize <= 128; csize <<= 1) {
    std::cout << "Compressed size = "
              << csize
              << "; Error = "
              << MSE(X, Uncompress(Compress(X, csize), 128))
              << std::endl;
  }

  std::cout << std::endl;
}
