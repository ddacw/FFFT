#pragma once

#include <algorithm>
#include <execution>
#include <thread>

#include "util.h"

namespace FFT {

void FFTSlow(jarray& X);

void FFTRec(jarray& X);

void FFTSeq(jarray& X, bool invert);

void FFTParallel(jarray& X, bool invert, size_t num_threads);

}  // namespace FFT
