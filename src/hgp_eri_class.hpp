#pragma once

#include "basis.hpp"
#include "eri_class.hpp"
#include <cstddef>

namespace nhfInt {

using nhfMath::Vec3d;

// use HGP method to calculate a eri class
// doi: 10.1063/1.455553
EriClass hgp_eri_class(
    const Basis &a, const Basis &b,
    const Basis &c, const Basis &d,
    const PairDataSet &pds
);

// generate (e0|f0)
EriClass hgp_vrr_contracted(
    std::size_t e, std::size_t f,
    const Basis &a, const Basis &b,
    const Basis &c, const Basis &d,
    const PairDataSet &pds
);

// generate [e0|f0]
void hgp_vrr_primitive(
    std::vector<double> &result, std::size_t e, std::size_t f,
    double zetaAB, double kAB, double Px, double Py, double Pz, double PAx, double PAy, double PAz,
    double zetaCD, double kCD, double Qx, double Qy, double Qz, double QCx, double QCy, double QCz
);

// (e0|xx) --> (ab|xx)
EriClass hgp_hrr_bra(
    const std::vector<EriClass> &eriVec, Vec3d v,
    std::size_t a, std::size_t b
);

// (xx|f0) --> (xx|cd)
EriClass hgp_hrr_ket(
    const std::vector<EriClass> &eriVec, Vec3d v,
    std::size_t c, std::size_t d
);

// single hgp hrr
std::vector<double> hgp_hrr(
    const std::vector<double> &hrrInp, 
    double x, double y, double z,
    std::size_t a, std::size_t b
);


} // namespace (nhfInt)