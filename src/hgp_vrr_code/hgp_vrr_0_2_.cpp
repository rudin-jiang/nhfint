#include "hgp_vrr.hpp"
#include "boysfun.hpp"
#include <vector>
#include <cmath>

namespace nhfInt {

using VecReal = std::vector<double>;

void hgp_vrr_0_2_(
    VecReal &a0b2, 
    double zetaAB, double kAB, double Px, double Py, double Pz,
    double PAx, double PAy, double PAz, 
    double zetaCD, double kCD, double Qx, double Qy, double Qz, 
    double QCx, double QCy, double QCz
) { 
    double zetaSum   = zetaAB + zetaCD;
    double invZeta   = 1.0 / zetaSum;
    double invZetaAB = 1.0 / zetaAB;
    double invZetaCD = 1.0 / zetaCD;
    double PQx = Px - Qx;
    double PQy = Py - Qy;
    double PQz = Pz - Qz;
    double Wx = (zetaAB * Px + zetaCD * Qx) * invZeta;
    double Wy = (zetaAB * Py + zetaCD * Qy) * invZeta;
    double Wz = (zetaAB * Pz + zetaCD * Qz) * invZeta;
    double WPx = Wx - Px;
    double WPy = Wy - Py;
    double WPz = Wz - Pz;
    double WQx = Wx - Qx;
    double WQy = Wy - Qy;
    double WQz = Wz - Qz;

    double len2PQ = PQx*PQx + PQy*PQy + PQz*PQz;
    double T = zetaAB * zetaCD * invZeta * len2PQ;
    double preBoysFunCoeff = std::sqrt(invZeta) * kAB * kCD;

    double a0b0[3];
    double a0b1[6];

    a0b0[0] = preBoysFunCoeff * boysfun(0, T);
    a0b0[1] = preBoysFunCoeff * boysfun(1, T);
    a0b0[2] = preBoysFunCoeff * boysfun(2, T);

    a0b1[0] = QCz * a0b0[0] + WQz * a0b0[1];
    a0b1[1] = QCy * a0b0[0] + WQy * a0b0[1];
    a0b1[2] = QCx * a0b0[0] + WQx * a0b0[1];
    a0b1[3] = QCz * a0b0[1] + WQz * a0b0[2];
    a0b1[4] = QCy * a0b0[1] + WQy * a0b0[2];
    a0b1[5] = QCx * a0b0[1] + WQx * a0b0[2];
    a0b2[0] = QCz * a0b1[0] + WQz * a0b1[3] + 0.5 * 1 * invZetaCD * (a0b0[0] - zetaAB * invZeta * a0b0[1]);
    a0b2[1] = QCy * a0b1[0] + WQy * a0b1[3];
    a0b2[2] = QCy * a0b1[1] + WQy * a0b1[4] + 0.5 * 1 * invZetaCD * (a0b0[0] - zetaAB * invZeta * a0b0[1]);
    a0b2[3] = QCx * a0b1[0] + WQx * a0b1[3];
    a0b2[4] = QCx * a0b1[1] + WQx * a0b1[4];
    a0b2[5] = QCx * a0b1[2] + WQx * a0b1[5] + 0.5 * 1 * invZetaCD * (a0b0[0] - zetaAB * invZeta * a0b0[1]);
} // function (hgp_vrr_0_2_)

} // namespace (nhfInt)

