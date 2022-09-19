#include "hgp_eri_class.hpp"
#include "basis.hpp"
#include "eri_class.hpp"
#include "vec3d.hpp"
#include <vector>
#include <cstddef>
#include <cassert>

namespace nhfInt {


EriClass hgp_eri_class(
    const Basis &a, const Basis &b,
    const Basis &c, const Basis &d,
    const PairDataSet &pds
) {
    assert(a.angMom >= b.angMom);
    assert(c.angMom >= d.angMom);

    Vec3d AB = a.centre - b.centre;
    Vec3d CD = c.centre - d.centre;

    // (e0|f0)  f: c.angMom, ..., c.angMom + d.angMom
    std::vector<EriClass> e0f0(d.angMom + 1);

    // (e0|cd)  e: a.angMom, ..., a.angMom + b.angMom
    std::vector<EriClass> e0cd(b.angMom + 1);

    for (std::size_t i = 0; i <= b.angMom; ++i) {
        std::size_t e = i + a.angMom;
        for (std::size_t j = 0; j <= d.angMom; ++j) {
            std::size_t f = j + c.angMom;
            e0f0[j] = hgp_vrr_contracted(e, f, a, b, c, d, pds);
        }

        e0cd[i] = hgp_hrr_ket(e0f0, CD, c.angMom, d.angMom);
    }

    EriClass abcd = hgp_hrr_bra(e0cd, AB, a.angMom, b.angMom);

    return abcd;
}


EriClass hgp_vrr_contracted(
    std::size_t e, std::size_t f,
    const Basis &a, const Basis &b,
    const Basis &c, const Basis &d,
    const PairDataSet &pds
) {
    EriClass ret(e, 0, f, 0);

    // store each primitive integral
    std::size_t nInt = ret.size();
    std::vector<double> primitiveIntegral(nInt);

    for (std::size_t i = 0; i < a.gauss_size(); ++i) {
    for (std::size_t j = 0; j < b.gauss_size(); ++j) {
        std::size_t ij = idx2(i + a.gsBeg, j + b.gsBeg);
        const PairData &pdij = pds[ij];
        Vec3d PA = pdij.P - a.centre;

        for (std::size_t k = 0; k < c.gauss_size(); ++k) {
        for (std::size_t l = 0; l < d.gauss_size(); ++l) {
            std::size_t kl = idx2(k + c.gsBeg, l + d.gsBeg);
            const PairData &pdkl = pds[kl];
            Vec3d QC = pdkl.P - c.centre;

            hgp_vrr_primitive(
                primitiveIntegral, e, f,
                pdij.zeta, pdij.Kab, pdij.P.x, pdij.P.y, pdij.P.z, PA.x, PA.y, PA.z,
                pdkl.zeta, pdkl.Kab, pdkl.P.x, pdkl.P.y, pdkl.P.z, QC.x, QC.y, QC.z
            );

            for (std::size_t ic = 0; ic < nInt; ++ic) {
                ret[ic] += pdij.Cab * pdkl.Cab * primitiveIntegral[ic];
            }
        }}
    }}

    return ret;
}

// result has already been allocateed when hgp_vrr_primitive is called 
void hgp_vrr_primitive(
    std::vector<double> result, std::size_t e, std::size_t f,
    double zetaAB, double kAB, double Px, double Py, double Pz, double PAx, double PAy, double PAz,
    double zetaCD, double kCD, double Qx, double Qy, double Qz, double QCx, double QCy, double QCz
) {
    assert(e >= f);
    if (e = 0 && f = 0) hgp_vrr_0_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 1 && f = 0) hgp_vrr_1_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 1 && f = 1) hgp_vrr_1_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 2 && f = 0) hgp_vrr_2_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 2 && f = 1) hgp_vrr_2_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 2 && f = 2) hgp_vrr_2_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 3 && f = 0) hgp_vrr_3_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 3 && f = 1) hgp_vrr_3_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 3 && f = 2) hgp_vrr_3_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 3 && f = 3) hgp_vrr_3_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 4 && f = 0) hgp_vrr_4_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 4 && f = 1) hgp_vrr_4_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 4 && f = 2) hgp_vrr_4_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 4 && f = 3) hgp_vrr_4_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 4 && f = 4) hgp_vrr_4_4_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 0) hgp_vrr_5_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 1) hgp_vrr_5_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 2) hgp_vrr_5_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 3) hgp_vrr_5_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 4) hgp_vrr_5_4_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 5 && f = 5) hgp_vrr_5_5_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 0) hgp_vrr_6_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 1) hgp_vrr_6_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 2) hgp_vrr_6_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 3) hgp_vrr_6_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 4) hgp_vrr_6_4_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 5) hgp_vrr_6_5_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 6 && f = 6) hgp_vrr_6_6_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 0) hgp_vrr_7_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 1) hgp_vrr_7_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 2) hgp_vrr_7_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 3) hgp_vrr_7_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 4) hgp_vrr_7_4_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 5) hgp_vrr_7_5_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 6) hgp_vrr_7_6_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 7 && f = 7) hgp_vrr_7_7_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 0) hgp_vrr_8_0_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 1) hgp_vrr_8_1_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 2) hgp_vrr_8_2_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 3) hgp_vrr_8_3_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 4) hgp_vrr_8_4_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 5) hgp_vrr_8_5_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 6) hgp_vrr_8_6_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 7) hgp_vrr_8_7_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
    if (e = 8 && f = 8) hgp_vrr_8_8_(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);
}



// (e0|xx) --> (ab|xx)
EriClass hgp_hrr_bra(
    const std::vector<EriClass> &eriVec, Vec3d v,
    std::size_t a, std::size_t b
) {
    assert(eriVec.size() == b + 1);
    assert(eriVec.front().angA == a);
    assert(eriVec.back().angA == a + b);

    // determine the angular momentum of the returned EriClass
    std::size_t angC = eriVec.front().angC;
    std::size_t angD = eriVec.front().angD;

    // result EriClass
    EriClass ret(a, b, angC, angD);

    // size of hrr input vector
    // mma: Sum[(n+1)*(n+2)/2, {n, a, a+b}]
    std::size_t hrrLen = (1+b) * (6+9*a+3*a*a+5*b+3*a*b+b*b) / 6;
    std::vector<double> hrrInp(hrrLen);

    for (std::size_t j = 0; j < ret.nbsC * ret.nbsD; ++j) {
        
        // create hrr input vector
        std::size_t hrrPos = 0;
        for (std::size_t ie = 0; ie < eriVec.size(); ++ie) {
            std::size_t numEri = eriVec[ie].nbsA * eriVec[ie].nbsB;
            for (std::size_t i = 0; i < numEri; ++i) {
                hrrInp[hrrPos++] = eriVec[ie](i, j);
            }
        }

        // calculate (ab|xx)
        std::vector<double> abxx = hgp_hrr(hrrInp, v.x, v.y, v.z, a, b);

        // copy data to returned EriClass
        for (std::size_t i = 0; i < abxx.size(); ++i) {
            ret(i, j) = abxx[i];
        }
    }

    return ret;
}


// (xx|f0) --> (xx|cd)
EriClass hgp_hrr_ket(
    const std::vector<EriClass> &eriVec, Vec3d v,
    std::size_t c, std::size_t d
) {
    assert(eriVec.size() == d + 1);
    assert(eriVec.front().angC == c);
    assert(eriVec.back().angC == c + d);

    // determine the angular momentum of the returned EriClass
    std::size_t angA = eriVec.front().angA;
    std::size_t angB = eriVec.front().angB;

    // result EriClass
    EriClass ret(angA, angB, c, d);

    // size of hrr input vector
    // mma: Sum[(n+1)*(n+2)/2, {n, c, c+d}]
    std::size_t hrrLen = (1+d) * (6+9*c+3*c*c+5*d+3*c*d+d*d) / 6;
    std::vector<double> hrrInp(hrrLen);

    for (std::size_t i = 0; i < ret.nbsA * ret.nbsB; ++i) {

        // create hrr input vector
        std::size_t hrrPos = 0;
        for (std::size_t ie = 0; ie < eriVec.size(); ++ie) {
            std::size_t numEri = eriVec[ie].nbsC * eriVec[ie].nbsD;
            for (std::size_t j = 0; j < numEri; ++j) {
                hrrInp[hrrPos++] = eriVec[ie](i, j);
            }
        }

        // calculate (xx|cd)
        std::vector<double> xxcd = hgp_hrr(hrrInp, v.x, v.y, v.z, c, d);

        // copy data to returned EriClass
        for (std::size_t j = 0; j < xxcd.size(); ++j) {
            ret(i, j) = abxx[j];
        }
    }

    return ret;
}

// single hgp hrr
std::vector<double> hgp_hrr(
    const std::vector<double> &hrrInp, 
    double x, double y, double z,
    std::size_t a, std::size_t b
) {
    assert(a >= b);
    if (a = 0 && b = 0) return hgp_hrr_0_0_(hrrInp, x, y, z);
    if (a = 1 && b = 0) return hgp_hrr_1_0_(hrrInp, x, y, z);
    if (a = 1 && b = 1) return hgp_hrr_1_1_(hrrInp, x, y, z);
    if (a = 2 && b = 0) return hgp_hrr_2_0_(hrrInp, x, y, z);
    if (a = 2 && b = 1) return hgp_hrr_2_1_(hrrInp, x, y, z);
    if (a = 2 && b = 2) return hgp_hrr_2_2_(hrrInp, x, y, z);
    if (a = 3 && b = 0) return hgp_hrr_3_0_(hrrInp, x, y, z);
    if (a = 3 && b = 1) return hgp_hrr_3_1_(hrrInp, x, y, z);
    if (a = 3 && b = 2) return hgp_hrr_3_2_(hrrInp, x, y, z);
    if (a = 3 && b = 3) return hgp_hrr_3_3_(hrrInp, x, y, z);
    if (a = 4 && b = 0) return hgp_hrr_4_0_(hrrInp, x, y, z);
    if (a = 4 && b = 1) return hgp_hrr_4_1_(hrrInp, x, y, z);
    if (a = 4 && b = 2) return hgp_hrr_4_2_(hrrInp, x, y, z);
    if (a = 4 && b = 3) return hgp_hrr_4_3_(hrrInp, x, y, z);
    if (a = 4 && b = 4) return hgp_hrr_4_4_(hrrInp, x, y, z);
    return std::vector<double>();
}


}   // namespace (nhfInt)

