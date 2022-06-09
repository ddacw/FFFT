#pragma once

#include <algorithm>
#include <execution>
#include <thread>

#include "utils.h"

namespace FFT {

void FFTSeq(jarray& X, bool invert);

void FFTParallel(jarray& X, bool invert, size_t num_threads);

}  // namespace FFT
