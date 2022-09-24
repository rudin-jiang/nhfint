#include "basis.hpp"
#include "vec3d.hpp"
#include "stropt.hpp"
#include <vector>
#include <string>
#include <cstddef>
#include <cmath>
#include <fstream>

// pow(2.0, 0.5) * pow(pi, 1.25), used in KAB and KCD
// in function hgp_vrr_contracted
static const double preCoeff = 5.9149671727956132;

namespace nhfInt {

Basis::Basis()
: bsBeg(0), gsBeg(0), angMom(0) {}

Basis::Basis(
    std::size_t bsBeg, std::size_t gsBeg, std::size_t angMom,
    const VecReal &alpha, const VecReal &coeff, const Vec3d &centre
) : bsBeg(bsBeg), gsBeg(gsBeg), angMom(angMom),
    alpha(alpha), coeff(coeff), centre(centre) {}


/*          BasisSet            */
BasisSet::BasisSet() {}

BasisSet::BasisSet(std::size_t nBasis)
: bsData(std::vector<Basis>(nBasis)) {}

// BasisSet::BasisSet(std::string basisFile, std::string atomType, const Vec3d &centre) {}


Basis   BasisSet::operator[](std::size_t i) const
{ return bsData[i]; }

Basis&  BasisSet::operator[](std::size_t i)
{ return bsData[i]; }

std::size_t BasisSet::basis_size() const{
    std::size_t nBs = 0;
    for (const Basis &b : bsData)
        nBs += b.basis_size();
    return nBs;
}

std::size_t BasisSet::gauss_size() const {
    std::size_t nGs = 0;
    for (const Basis &b : bsData)
        nGs += b.gauss_size();
    return nGs;
}

std::size_t BasisSet::size() const 
{ return bsData.size(); }


/*          PairData            */
PairData::PairData()
: zeta(0.0), Kab(0.0), Cab(0.0) {}

PairData::PairData(double zeta, double Kab, double Cab, const Vec3d &P)
: zeta(zeta), Kab(Kab), Cab(Cab), P(P) {}


/*          PairDataSet         */
PairDataSet::PairDataSet() {}

PairDataSet::PairDataSet(std::size_t nPairData)
: pdData(std::vector<PairData>(nPairData)) {}

PairDataSet::PairDataSet(const BasisSet &bss) {
    std::size_t nGs = bss.gauss_size();
    std::size_t nPr = idx2(nGs-1, nGs-1) + 1;
    pdData = std::vector<PairData>(nPr);

    for (std::size_t i = 0; i < bss.size(); ++i) {
        const Basis &a = bss[i];

    for (std::size_t j = 0; j <= i; ++j) {
        const Basis &b = bss[j];

        for (std::size_t ia = 0; ia < a.gauss_size(); ++ia) {
        for (std::size_t jb = 0; jb < b.gauss_size(); ++jb) {

            double zeta = a.alpha[ia] + b.alpha[jb];
            double invz = 1.0 / zeta;
            double len2 = (a.centre - b.centre).len2();
            double inp  = -a.alpha[ia] * b.alpha[jb] * invz * len2;
            double Kab  = preCoeff * invz * std::exp(inp);
            double Cab  = a.coeff[ia] * b.coeff[jb];
            Vec3d  P    = (a.alpha[ia] * a.centre + b.alpha[jb] * b.centre) * invz;

            // position of this pair in PairDataSet
            std::size_t iPos = idx2(a.gsBeg + ia, b.gsBeg + jb);
            pdData[iPos] = PairData(zeta, Kab, Cab, P);
        }}

    }}
}

PairData PairDataSet::operator[](std::size_t i) const
{ return pdData[i]; }

PairData& PairDataSet::operator[](std::size_t i)
{ return pdData[i]; }

std::size_t PairDataSet::size() const
{ return pdData.size(); }


/*          idx2 and idx4           */
std::size_t idx2(std::size_t i, std::size_t j)
{ return  i>j ? i * (i+1) / 2 + j : j * (j+1) / 2 + i; }

std::size_t idx4(std::size_t i, std::size_t j, 
                 std::size_t k, std::size_t l)
{ return  idx2(idx2(i,j), idx2(k,l)); }


}   // namespace (nhfInt)