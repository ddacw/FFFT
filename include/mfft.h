#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "util.h"

namespace MFFT {

void Transform(jarray& X, bool invert, size_t num_threads);

}
