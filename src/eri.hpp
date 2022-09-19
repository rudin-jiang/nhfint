#pragma once

#include "basis.hpp"
#include <vector>

namespace nhfInt {

// calculate eri using HGP method
std::vector<double> hgp_eri(const BasisSet &bs);


}   // namespace (nhfInt)