#include "hgp_vrr.hpp"
#include "boysfun.hpp"
#include <vector>
#include <cmath>

namespace nhfInt {

using VecReal = std::vector<double>;

void hgp_vrr_1_2_(
    VecReal &a1b2, 
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

    double a0b0[ 4];
    double a0b1[ 9];
    double a1b0[ 9];
    double a0b2[12];
    double a1b1[18];

    a0b0[0] = preBoysFunCoeff * boysfun(0, T);
    a0b0[1] = preBoysFunCoeff * boysfun(1, T);
    a0b0[2] = preBoysFunCoeff * boysfun(2, T);
    a0b0[3] = preBoysFunCoeff * boysfun(3, T);

    a0b1[ 0] = QCz * a0b0[ 0] + WQz * a0b0[ 1];
    a0b1[ 1] = QCy * a0b0[ 0] + WQy * a0b0[ 1];
    a0b1[ 2] = QCx * a0b0[ 0] + WQx * a0b0[ 1];
    a0b1[ 3] = QCz * a0b0[ 1] + WQz * a0b0[ 2];
    a0b1[ 4] = QCy * a0b0[ 1] + WQy * a0b0[ 2];
    a0b1[ 5] = QCx * a0b0[ 1] + WQx * a0b0[ 2];
    a0b1[ 6] = QCz * a0b0[ 2] + WQz * a0b0[ 3];
    a0b1[ 7] = QCy * a0b0[ 2] + WQy * a0b0[ 3];
    a0b1[ 8] = QCx * a0b0[ 2] + WQx * a0b0[ 3];
    a1b0[ 1] = PAy * a0b0[ 0] + WPy * a0b0[ 1];
    a1b0[ 2] = PAx * a0b0[ 0] + WPx * a0b0[ 1];
    a1b0[ 4] = PAy * a0b0[ 1] + WPy * a0b0[ 2];
    a1b0[ 5] = PAx * a0b0[ 1] + WPx * a0b0[ 2];
    a1b0[ 7] = PAy * a0b0[ 2] + WPy * a0b0[ 3];
    a1b0[ 8] = PAx * a0b0[ 2] + WPx * a0b0[ 3];
    a0b2[ 0] = QCz * a0b1[ 0] + WQz * a0b1[ 3] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 1] = QCy * a0b1[ 0] + WQy * a0b1[ 3];
    a0b2[ 2] = QCy * a0b1[ 1] + WQy * a0b1[ 4] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 4] = QCx * a0b1[ 1] + WQx * a0b1[ 4];
    a0b2[ 5] = QCx * a0b1[ 2] + WQx * a0b1[ 5] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 6] = QCz * a0b1[ 3] + WQz * a0b1[ 6] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 7] = QCy * a0b1[ 3] + WQy * a0b1[ 6];
    a0b2[ 8] = QCy * a0b1[ 4] + WQy * a0b1[ 7] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[10] = QCx * a0b1[ 4] + WQx * a0b1[ 7];
    a0b2[11] = QCx * a0b1[ 5] + WQx * a0b1[ 8] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a1b1[ 0] = PAz * a0b1[ 0] + WPz * a0b1[ 3] + 0.5 * 1 * invZeta * a0b0[ 1];
    a1b1[ 3] = PAy * a0b1[ 0] + WPy * a0b1[ 3];
    a1b1[ 6] = PAx * a0b1[ 0] + WPx * a0b1[ 3];
    a1b1[ 9] = PAz * a0b1[ 3] + WPz * a0b1[ 6] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[12] = PAy * a0b1[ 3] + WPy * a0b1[ 6];
    a1b1[15] = PAx * a0b1[ 3] + WPx * a0b1[ 6];
    a1b2[ 0] = PAz * a0b2[ 0] + WPz * a0b2[ 6] + 0.5 * 2 * invZeta * a0b1[ 3];
    a1b2[ 1] = PAz * a0b2[ 1] + WPz * a0b2[ 7] + 0.5 * 1 * invZeta * a0b1[ 4];
    a1b2[ 2] = PAz * a0b2[ 2] + WPz * a0b2[ 8];
    a1b2[ 3] = QCx * a1b1[ 0] + WQx * a1b1[ 9];
    a1b2[ 4] = PAz * a0b2[ 4] + WPz * a0b2[10];
    a1b2[ 5] = PAz * a0b2[ 5] + WPz * a0b2[11];
    a1b2[ 6] = PAy * a0b2[ 0] + WPy * a0b2[ 6];
    a1b2[ 7] = PAy * a0b2[ 1] + WPy * a0b2[ 7] + 0.5 * 1 * invZeta * a0b1[ 3];
    a1b2[ 8] = PAy * a0b2[ 2] + WPy * a0b2[ 8] + 0.5 * 2 * invZeta * a0b1[ 4];
    a1b2[ 9] = QCx * a1b1[ 3] + WQx * a1b1[12];
    a1b2[10] = PAy * a0b2[ 4] + WPy * a0b2[10] + 0.5 * 1 * invZeta * a0b1[ 5];
    a1b2[11] = PAy * a0b2[ 5] + WPy * a0b2[11];
    a1b2[12] = PAx * a0b2[ 0] + WPx * a0b2[ 6];
    a1b2[13] = PAx * a0b2[ 1] + WPx * a0b2[ 7];
    a1b2[14] = PAx * a0b2[ 2] + WPx * a0b2[ 8];
    a1b2[15] = QCx * a1b1[ 6] + WQx * a1b1[15] + 0.5 * 1 * invZeta * a0b1[ 3];
    a1b2[16] = PAx * a0b2[ 4] + WPx * a0b2[10] + 0.5 * 1 * invZeta * a0b1[ 4];
    a1b2[17] = PAx * a0b2[ 5] + WPx * a0b2[11] + 0.5 * 2 * invZeta * a0b1[ 5];
} // function (hgp_vrr_1_2_)

} // namespace (nhfInt)

