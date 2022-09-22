#pragma once

#include <vector>

namespace nhfInt {

using VecReal = std::vector<double>;

VecReal hgp_hrr_0_0_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_1_0_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_1_1_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_2_0_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_2_1_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_2_2_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_3_0_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_3_1_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_3_2_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_3_3_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_4_0_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_4_1_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_4_2_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_4_3_(const VecReal &hrrInp, double x, double y, double z);
VecReal hgp_hrr_4_4_(const VecReal &hrrInp, double x, double y, double z);

}