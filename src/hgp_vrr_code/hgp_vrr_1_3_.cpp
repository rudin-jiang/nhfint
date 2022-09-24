#include "hgp_vrr.hpp"
#include "boysfun.hpp"
#include <vector>
#include <cmath>

namespace nhfInt {

using VecReal = std::vector<double>;

void hgp_vrr_1_3_(
    VecReal &a1b3, 
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
    double a1b0[12];
    double a0b2[18];
    double a1b1[27];
    double a0b3[20];
    double a1b2[36];

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
    a1b0[ 0] = PAz * a0b0[ 0] + WPz * a0b0[ 1];
    a1b0[ 1] = PAy * a0b0[ 0] + WPy * a0b0[ 1];
    a1b0[ 3] = PAz * a0b0[ 1] + WPz * a0b0[ 2];
    a1b0[ 4] = PAy * a0b0[ 1] + WPy * a0b0[ 2];
    a1b0[ 6] = PAz * a0b0[ 2] + WPz * a0b0[ 3];
    a1b0[ 7] = PAy * a0b0[ 2] + WPy * a0b0[ 3];
    a1b0[10] = PAy * a0b0[ 3] + WPy * a0b0[ 4];
    a0b2[ 0] = QCz * a0b1[ 0] + WQz * a0b1[ 3] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 1] = QCy * a0b1[ 0] + WQy * a0b1[ 3];
    a0b2[ 2] = QCy * a0b1[ 1] + WQy * a0b1[ 4] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 4] = QCx * a0b1[ 1] + WQx * a0b1[ 4];
    a0b2[ 5] = QCx * a0b1[ 2] + WQx * a0b1[ 5] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 6] = QCz * a0b1[ 3] + WQz * a0b1[ 6] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 7] = QCy * a0b1[ 3] + WQy * a0b1[ 6];
    a0b2[ 8] = QCy * a0b1[ 4] + WQy * a0b1[ 7] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 9] = QCx * a0b1[ 3] + WQx * a0b1[ 6];
    a0b2[10] = QCx * a0b1[ 4] + WQx * a0b1[ 7];
    a0b2[11] = QCx * a0b1[ 5] + WQx * a0b1[ 8] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[12] = QCz * a0b1[ 6] + WQz * a0b1[ 9] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b2[13] = QCy * a0b1[ 6] + WQy * a0b1[ 9];
    a0b2[14] = QCy * a0b1[ 7] + WQy * a0b1[10] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b2[15] = QCx * a0b1[ 6] + WQx * a0b1[ 9];
    a0b2[16] = QCx * a0b1[ 7] + WQx * a0b1[10];
    a0b2[17] = QCx * a0b1[ 8] + WQx * a0b1[11] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a1b1[ 0] = PAz * a0b1[ 0] + WPz * a0b1[ 3] + 0.5 * 1 * invZeta * a0b0[ 1];
    a1b1[ 1] = PAz * a0b1[ 1] + WPz * a0b1[ 4];
    a1b1[ 3] = PAy * a0b1[ 0] + WPy * a0b1[ 3];
    a1b1[ 4] = PAy * a0b1[ 1] + WPy * a0b1[ 4] + 0.5 * 1 * invZeta * a0b0[ 1];
    a1b1[ 5] = PAy * a0b1[ 2] + WPy * a0b1[ 5];
    a1b1[ 9] = PAz * a0b1[ 3] + WPz * a0b1[ 6] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[10] = PAz * a0b1[ 4] + WPz * a0b1[ 7];
    a1b1[12] = PAy * a0b1[ 3] + WPy * a0b1[ 6];
    a1b1[13] = PAy * a0b1[ 4] + WPy * a0b1[ 7] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[14] = PAy * a0b1[ 5] + WPy * a0b1[ 8];
    a1b1[18] = PAz * a0b1[ 6] + WPz * a0b1[ 9] + 0.5 * 1 * invZeta * a0b0[ 3];
    a1b1[22] = PAy * a0b1[ 7] + WPy * a0b1[10] + 0.5 * 1 * invZeta * a0b0[ 3];
    a1b1[23] = PAy * a0b1[ 8] + WPy * a0b1[11];
    a0b3[ 0] = QCz * a0b2[ 0] + WQz * a0b2[ 6] + 0.5 * 2 * invZetaCD * (a0b1[ 0] - zetaAB * invZeta * a0b1[ 3]);
    a0b3[ 1] = QCy * a0b2[ 0] + WQy * a0b2[ 6];
    a0b3[ 2] = QCy * a0b2[ 1] + WQy * a0b2[ 7] + 0.5 * 1 * invZetaCD * (a0b1[ 0] - zetaAB * invZeta * a0b1[ 3]);
    a0b3[ 3] = QCy * a0b2[ 2] + WQy * a0b2[ 8] + 0.5 * 2 * invZetaCD * (a0b1[ 1] - zetaAB * invZeta * a0b1[ 4]);
    a0b3[ 4] = QCx * a0b2[ 0] + WQx * a0b2[ 6];
    a0b3[ 5] = QCx * a0b2[ 1] + WQx * a0b2[ 7];
    a0b3[ 6] = QCx * a0b2[ 2] + WQx * a0b2[ 8];
    a0b3[ 7] = QCz * a0b2[ 5] + WQz * a0b2[11];
    a0b3[ 8] = QCx * a0b2[ 4] + WQx * a0b2[10] + 0.5 * 1 * invZetaCD * (a0b1[ 1] - zetaAB * invZeta * a0b1[ 4]);
    a0b3[ 9] = QCx * a0b2[ 5] + WQx * a0b2[11] + 0.5 * 2 * invZetaCD * (a0b1[ 2] - zetaAB * invZeta * a0b1[ 5]);
    a0b3[10] = QCz * a0b2[ 6] + WQz * a0b2[12] + 0.5 * 2 * invZetaCD * (a0b1[ 3] - zetaAB * invZeta * a0b1[ 6]);
    a0b3[11] = QCy * a0b2[ 6] + WQy * a0b2[12];
    a0b3[12] = QCy * a0b2[ 7] + WQy * a0b2[13] + 0.5 * 1 * invZetaCD * (a0b1[ 3] - zetaAB * invZeta * a0b1[ 6]);
    a0b3[13] = QCy * a0b2[ 8] + WQy * a0b2[14] + 0.5 * 2 * invZetaCD * (a0b1[ 4] - zetaAB * invZeta * a0b1[ 7]);
    a0b3[14] = QCx * a0b2[ 6] + WQx * a0b2[12];
    a0b3[15] = QCx * a0b2[ 7] + WQx * a0b2[13];
    a0b3[16] = QCx * a0b2[ 8] + WQx * a0b2[14];
    a0b3[17] = QCx * a0b2[ 9] + WQx * a0b2[15] + 0.5 * 1 * invZetaCD * (a0b1[ 3] - zetaAB * invZeta * a0b1[ 6]);
    a0b3[18] = QCx * a0b2[10] + WQx * a0b2[16] + 0.5 * 1 * invZetaCD * (a0b1[ 4] - zetaAB * invZeta * a0b1[ 7]);
    a0b3[19] = QCx * a0b2[11] + WQx * a0b2[17] + 0.5 * 2 * invZetaCD * (a0b1[ 5] - zetaAB * invZeta * a0b1[ 8]);
    a1b2[ 0] = PAz * a0b2[ 0] + WPz * a0b2[ 6] + 0.5 * 2 * invZeta * a0b1[ 3];
    a1b2[ 2] = PAz * a0b2[ 2] + WPz * a0b2[ 8];
    a1b2[ 3] = QCx * a1b1[ 0] + WQx * a1b1[ 9];
    a1b2[ 7] = PAy * a0b2[ 1] + WPy * a0b2[ 7] + 0.5 * 1 * invZeta * a0b1[ 3];
    a1b2[ 8] = PAy * a0b2[ 2] + WPy * a0b2[ 8] + 0.5 * 2 * invZeta * a0b1[ 4];
    a1b2[ 9] = QCx * a1b1[ 3] + WQx * a1b1[12];
    a1b2[10] = PAy * a0b2[ 4] + WPy * a0b2[10] + 0.5 * 1 * invZeta * a0b1[ 5];
    a1b2[12] = PAx * a0b2[ 0] + WPx * a0b2[ 6];
    a1b2[17] = PAx * a0b2[ 5] + WPx * a0b2[11] + 0.5 * 2 * invZeta * a0b1[ 5];
    a1b2[18] = PAz * a0b2[ 6] + WPz * a0b2[12] + 0.5 * 2 * invZeta * a0b1[ 6];
    a1b2[20] = PAz * a0b2[ 8] + WPz * a0b2[14];
    a1b2[21] = PAz * a0b2[ 9] + WPz * a0b2[15] + 0.5 * 1 * invZeta * a0b1[ 8];
    a1b2[25] = PAy * a0b2[ 7] + WPy * a0b2[13] + 0.5 * 1 * invZeta * a0b1[ 6];
    a1b2[26] = PAy * a0b2[ 8] + WPy * a0b2[14] + 0.5 * 2 * invZeta * a0b1[ 7];
    a1b2[27] = PAy * a0b2[ 9] + WPy * a0b2[15];
    a1b2[28] = PAy * a0b2[10] + WPy * a0b2[16] + 0.5 * 1 * invZeta * a0b1[ 8];
    a1b2[30] = PAx * a0b2[ 6] + WPx * a0b2[12];
    a1b2[35] = PAx * a0b2[11] + WPx * a0b2[17] + 0.5 * 2 * invZeta * a0b1[ 8];
    a1b3[ 0] = PAz * a0b3[ 0] + WPz * a0b3[10] + 0.5 * 3 * invZeta * a0b2[ 6];
    a1b3[ 1] = PAz * a0b3[ 1] + WPz * a0b3[11] + 0.5 * 2 * invZeta * a0b2[ 7];
    a1b3[ 2] = PAz * a0b3[ 2] + WPz * a0b3[12] + 0.5 * 1 * invZeta * a0b2[ 8];
    a1b3[ 3] = PAz * a0b3[ 3] + WPz * a0b3[13];
    a1b3[ 4] = PAz * a0b3[ 4] + WPz * a0b3[14] + 0.5 * 2 * invZeta * a0b2[ 9];
    a1b3[ 5] = PAz * a0b3[ 5] + WPz * a0b3[15] + 0.5 * 1 * invZeta * a0b2[10];
    a1b3[ 6] = PAz * a0b3[ 6] + WPz * a0b3[16];
    a1b3[ 7] = PAz * a0b3[ 7] + WPz * a0b3[17] + 0.5 * 1 * invZeta * a0b2[11];
    a1b3[ 8] = PAz * a0b3[ 8] + WPz * a0b3[18];
    a1b3[ 9] = PAz * a0b3[ 9] + WPz * a0b3[19];
    a1b3[10] = PAy * a0b3[ 0] + WPy * a0b3[10];
    a1b3[11] = PAy * a0b3[ 1] + WPy * a0b3[11] + 0.5 * 1 * invZeta * a0b2[ 6];
    a1b3[12] = PAy * a0b3[ 2] + WPy * a0b3[12] + 0.5 * 2 * invZeta * a0b2[ 7];
    a1b3[13] = PAy * a0b3[ 3] + WPy * a0b3[13] + 0.5 * 3 * invZeta * a0b2[ 8];
    a1b3[14] = PAy * a0b3[ 4] + WPy * a0b3[14];
    a1b3[15] = PAy * a0b3[ 5] + WPy * a0b3[15] + 0.5 * 1 * invZeta * a0b2[ 9];
    a1b3[16] = PAy * a0b3[ 6] + WPy * a0b3[16] + 0.5 * 2 * invZeta * a0b2[10];
    a1b3[17] = PAy * a0b3[ 7] + WPy * a0b3[17];
    a1b3[18] = PAy * a0b3[ 8] + WPy * a0b3[18] + 0.5 * 1 * invZeta * a0b2[11];
    a1b3[19] = PAy * a0b3[ 9] + WPy * a0b3[19];
    a1b3[20] = PAx * a0b3[ 0] + WPx * a0b3[10];
    a1b3[21] = PAx * a0b3[ 1] + WPx * a0b3[11];
    a1b3[22] = PAx * a0b3[ 2] + WPx * a0b3[12];
    a1b3[23] = PAx * a0b3[ 3] + WPx * a0b3[13];
    a1b3[24] = PAx * a0b3[ 4] + WPx * a0b3[14] + 0.5 * 1 * invZeta * a0b2[ 6];
    a1b3[25] = PAx * a0b3[ 5] + WPx * a0b3[15] + 0.5 * 1 * invZeta * a0b2[ 7];
    a1b3[26] = PAx * a0b3[ 6] + WPx * a0b3[16] + 0.5 * 1 * invZeta * a0b2[ 8];
    a1b3[27] = PAx * a0b3[ 7] + WPx * a0b3[17] + 0.5 * 2 * invZeta * a0b2[ 9];
    a1b3[28] = PAx * a0b3[ 8] + WPx * a0b3[18] + 0.5 * 2 * invZeta * a0b2[10];
    a1b3[29] = PAx * a0b3[ 9] + WPx * a0b3[19] + 0.5 * 3 * invZeta * a0b2[11];
} // function (hgp_vrr_1_3_)

} // namespace (nhfInt)

