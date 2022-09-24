#include "hgp_vrr.hpp"
#include "boysfun.hpp"
#include <vector>
#include <cmath>

namespace nhfInt {

using VecReal = std::vector<double>;

void hgp_vrr_0_4_(
    VecReal &a0b4, 
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

    double a0b0[ 5];
    double a0b1[12];
    double a0b2[18];
    double a0b3[20];

    a0b0[0] = preBoysFunCoeff * boysfun(0, T);
    a0b0[1] = preBoysFunCoeff * boysfun(1, T);
    a0b0[2] = preBoysFunCoeff * boysfun(2, T);
    a0b0[3] = preBoysFunCoeff * boysfun(3, T);
    a0b0[4] = preBoysFunCoeff * boysfun(4, T);

    a0b1[ 0] = QCz * a0b0[ 0] + WQz * a0b0[ 1];
    a0b1[ 1] = QCy * a0b0[ 0] + WQy * a0b0[ 1];
    a0b1[ 2] = QCx * a0b0[ 0] + WQx * a0b0[ 1];
    a0b1[ 3] = QCz * a0b0[ 1] + WQz * a0b0[ 2];
    a0b1[ 4] = QCy * a0b0[ 1] + WQy * a0b0[ 2];
    a0b1[ 5] = QCx * a0b0[ 1] + WQx * a0b0[ 2];
    a0b1[ 6] = QCz * a0b0[ 2] + WQz * a0b0[ 3];
    a0b1[ 7] = QCy * a0b0[ 2] + WQy * a0b0[ 3];
    a0b1[ 8] = QCx * a0b0[ 2] + WQx * a0b0[ 3];
    a0b1[ 9] = QCz * a0b0[ 3] + WQz * a0b0[ 4];
    a0b1[10] = QCy * a0b0[ 3] + WQy * a0b0[ 4];
    a0b1[11] = QCx * a0b0[ 3] + WQx * a0b0[ 4];
    a0b2[ 0] = QCz * a0b1[ 0] + WQz * a0b1[ 3] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 1] = QCy * a0b1[ 0] + WQy * a0b1[ 3];
    a0b2[ 2] = QCy * a0b1[ 1] + WQy * a0b1[ 4] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 3] = QCx * a0b1[ 0] + WQx * a0b1[ 3];
    a0b2[ 4] = QCx * a0b1[ 1] + WQx * a0b1[ 4];
    a0b2[ 5] = QCx * a0b1[ 2] + WQx * a0b1[ 5] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 6] = QCz * a0b1[ 3] + WQz * a0b1[ 6] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 7] = QCy * a0b1[ 3] + WQy * a0b1[ 6];
    a0b2[ 8] = QCy * a0b1[ 4] + WQy * a0b1[ 7] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 9] = QCx * a0b1[ 3] + WQx * a0b1[ 6];
    a0b2[10] = QCx * a0b1[ 4] + WQx * a0b1[ 7];
    a0b2[11] = QCx * a0b1[ 5] + WQx * a0b1[ 8] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[12] = QCz * a0b1[ 6] + WQz * a0b1[ 9] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b2[14] = QCy * a0b1[ 7] + WQy * a0b1[10] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b2[16] = QCx * a0b1[ 7] + WQx * a0b1[10];
    a0b2[17] = QCx * a0b1[ 8] + WQx * a0b1[11] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b3[ 0] = QCz * a0b2[ 0] + WQz * a0b2[ 6] + 0.5 * 2 * invZetaCD * (a0b1[ 0] - zetaAB * invZeta * a0b1[ 3]);
    a0b3[ 1] = QCy * a0b2[ 0] + WQy * a0b2[ 6];
    a0b3[ 3] = QCy * a0b2[ 2] + WQy * a0b2[ 8] + 0.5 * 2 * invZetaCD * (a0b1[ 1] - zetaAB * invZeta * a0b1[ 4]);
    a0b3[ 6] = QCx * a0b2[ 2] + WQx * a0b2[ 8];
    a0b3[ 7] = QCx * a0b2[ 3] + WQx * a0b2[ 9] + 0.5 * 1 * invZetaCD * (a0b1[ 0] - zetaAB * invZeta * a0b1[ 3]);
    a0b3[ 9] = QCx * a0b2[ 5] + WQx * a0b2[11] + 0.5 * 2 * invZetaCD * (a0b1[ 2] - zetaAB * invZeta * a0b1[ 5]);
    a0b3[10] = QCz * a0b2[ 6] + WQz * a0b2[12] + 0.5 * 2 * invZetaCD * (a0b1[ 3] - zetaAB * invZeta * a0b1[ 6]);
    a0b3[11] = QCy * a0b2[ 6] + WQy * a0b2[12];
    a0b3[13] = QCy * a0b2[ 8] + WQy * a0b2[14] + 0.5 * 2 * invZetaCD * (a0b1[ 4] - zetaAB * invZeta * a0b1[ 7]);
    a0b3[16] = QCx * a0b2[ 8] + WQx * a0b2[14];
    a0b3[17] = QCz * a0b2[11] + WQz * a0b2[17];
    a0b3[19] = QCx * a0b2[11] + WQx * a0b2[17] + 0.5 * 2 * invZetaCD * (a0b1[ 5] - zetaAB * invZeta * a0b1[ 8]);
    a0b4[ 0] = QCz * a0b3[ 0] + WQz * a0b3[10] + 0.5 * 3 * invZetaCD * (a0b2[ 0] - zetaAB * invZeta * a0b2[ 6]);
    a0b4[ 1] = QCy * a0b3[ 0] + WQy * a0b3[10];
    a0b4[ 2] = QCy * a0b3[ 1] + WQy * a0b3[11] + 0.5 * 1 * invZetaCD * (a0b2[ 0] - zetaAB * invZeta * a0b2[ 6]);
    a0b4[ 3] = QCz * a0b3[ 3] + WQz * a0b3[13];
    a0b4[ 4] = QCy * a0b3[ 3] + WQy * a0b3[13] + 0.5 * 3 * invZetaCD * (a0b2[ 2] - zetaAB * invZeta * a0b2[ 8]);
    a0b4[ 5] = QCx * a0b3[ 0] + WQx * a0b3[10];
    a0b4[ 6] = QCx * a0b3[ 1] + WQx * a0b3[11];
    a0b4[ 7] = QCz * a0b3[ 6] + WQz * a0b3[16];
    a0b4[ 8] = QCx * a0b3[ 3] + WQx * a0b3[13];
    a0b4[ 9] = QCz * a0b3[ 7] + WQz * a0b3[17] + 0.5 * 1 * invZetaCD * (a0b2[ 5] - zetaAB * invZeta * a0b2[11]);
    a0b4[10] = QCy * a0b3[ 7] + WQy * a0b3[17];
    a0b4[11] = QCx * a0b3[ 6] + WQx * a0b3[16] + 0.5 * 1 * invZetaCD * (a0b2[ 2] - zetaAB * invZeta * a0b2[ 8]);
    a0b4[12] = QCx * a0b3[ 7] + WQx * a0b3[17] + 0.5 * 2 * invZetaCD * (a0b2[ 3] - zetaAB * invZeta * a0b2[ 9]);
    a0b4[13] = QCy * a0b3[ 9] + WQy * a0b3[19];
    a0b4[14] = QCx * a0b3[ 9] + WQx * a0b3[19] + 0.5 * 3 * invZetaCD * (a0b2[ 5] - zetaAB * invZeta * a0b2[11]);
} // function (hgp_vrr_0_4_)

} // namespace (nhfInt)

