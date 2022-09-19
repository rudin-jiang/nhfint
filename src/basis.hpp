#pragma once

#include <vector>
#include <string>
#include "vec3d.hpp"


namespace nhfInt {

using nhfMath::Vec3d;
using VecReal = std::vector<double>;

class Basis {
public:
    std::size_t     bsBeg;      // basis func begin id
    std::size_t     gsBeg;      // gauss func begin id
    std::size_t     angMom;     // angular momentum
    VecReal         alpha;      // gaussian exponent
    VecReal         coeff;      // contracted coeff
    Vec3d           centre;     // gaussian centre

    Basis();
    Basis(std::size_t bsBeg, std::size_t gsBeg, std::size_t angMom,
          const VecReal &alpha, const VecReal &coeff, const Vec3d &centre);
    Basis(std::string atomType, std::string basisFile, const Vec3d &centre);

    std::size_t basis_size() const
    { return (angMom+1) * (angMom+2) / 2; }

    std::size_t gauss_size() const
    { return alpha.size(); }
};


class BasisSet {
public:
    std::vector<Basis>  bsData;

    BasisSet();
    BasisSet(std::size_t nBasis);

    // access
    Basis   operator[](std::size_t i) const;
    Basis&  operator[](std::size_t i);

    std::size_t basis_size() const;
    std::size_t gauss_size() const;
    std::size_t size() const;
};


// PairData: record shell-pair data
// see HGP-1988 eq 8, 9 and 15
class PairData {
public:
    // Each pair of pimitive gaussians will have a PairData.
    // Order-independent pair data will be recorded.
    // That says (a,b) and (b,a) is the same

    double  zeta;   // alpha + beta (HGP eq 8)
    double  Kab;    // (HGP eq 15)
    double  Cab;    // coeffa * coeffb
    Vec3d   P;      // (HGP eq 9)

    PairData();
    PairData(double zeta, double Kab, double Cab, const Vec3d &P);
};


class PairDataSet {
public:
    std::vector<PairData> pdData;

    PairDataSet();
    PairDataSet(std::size_t nPairData);
    PairDataSet(const BasisSet &bss);

    PairData    operator[](std::size_t i) const;
    PairData&   operator[](std::size_t i);

    std::size_t size() const;
};


// When we use a one-dimensional array to store a symmetric matrix, 
// idx2 calculates the position of the matrix element in the array.
std::size_t idx2(std::size_t i, std::size_t j);
std::size_t idx4(std::size_t i, std::size_t j, 
                 std::size_t k, std::size_t l);

}   // namespace (nhfInt)
