#pragma once

#include <complex>
#include <iostream>
#include <utility>
#include <vector>

typedef std::complex<double> cmplx;
typedef std::vector<cmplx> jarray;
const size_t DATA_SIZE = 1 << 18;

cmplx GetW(int N, bool invert);

int CeilBit(int n);

int GetBit(int mask, int n);

int RevBit(int mask, int count);

int PadZero(jarray& X);
