#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "util.h"

namespace MFFT {

void Transform(jarray& X, size_t num_threads);

}
