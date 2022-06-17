#pragma once

#include <algorithm>
#include <execution>
#include <thread>

#include "util.h"

namespace FFT {

void DFT(jarray& X, bool invert);

void FFTRec(jarray& X, bool invert);

void FFTSeq(jarray& X, bool invert);

void FFTParallel(jarray& X, bool invert, size_t num_threads);

}  // namespace FFT
