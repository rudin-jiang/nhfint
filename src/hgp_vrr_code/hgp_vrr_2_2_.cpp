#include "hgp_vrr.hpp"
#include "boysfun.hpp"
#include <vector>
#include <cmath>

namespace nhfInt {

using VecReal = std::vector<double>;

void hgp_vrr_2_2_(
    VecReal &a2b2, 
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
    double a2b0[18];
    double a1b2[36];
    double a2b1[36];

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
    a1b0[ 2] = PAx * a0b0[ 0] + WPx * a0b0[ 1];
    a1b0[ 3] = PAz * a0b0[ 1] + WPz * a0b0[ 2];
    a1b0[ 4] = PAy * a0b0[ 1] + WPy * a0b0[ 2];
    a1b0[ 5] = PAx * a0b0[ 1] + WPx * a0b0[ 2];
    a1b0[ 6] = PAz * a0b0[ 2] + WPz * a0b0[ 3];
    a1b0[ 7] = PAy * a0b0[ 2] + WPy * a0b0[ 3];
    a1b0[ 8] = PAx * a0b0[ 2] + WPx * a0b0[ 3];
    a1b0[ 9] = PAz * a0b0[ 3] + WPz * a0b0[ 4];
    a1b0[10] = PAy * a0b0[ 3] + WPy * a0b0[ 4];
    a1b0[11] = PAx * a0b0[ 3] + WPx * a0b0[ 4];
    a0b2[ 0] = QCz * a0b1[ 0] + WQz * a0b1[ 3] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 2] = QCy * a0b1[ 1] + WQy * a0b1[ 4] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 5] = QCx * a0b1[ 2] + WQx * a0b1[ 5] + 0.5 * 1 * invZetaCD * (a0b0[ 0] - zetaAB * invZeta * a0b0[ 1]);
    a0b2[ 6] = QCz * a0b1[ 3] + WQz * a0b1[ 6] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[ 7] = QCy * a0b1[ 3] + WQy * a0b1[ 6];
    a0b2[ 8] = QCy * a0b1[ 4] + WQy * a0b1[ 7] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[11] = QCx * a0b1[ 5] + WQx * a0b1[ 8] + 0.5 * 1 * invZetaCD * (a0b0[ 1] - zetaAB * invZeta * a0b0[ 2]);
    a0b2[12] = QCz * a0b1[ 6] + WQz * a0b1[ 9] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a0b2[13] = QCy * a0b1[ 6] + WQy * a0b1[ 9];
    a0b2[17] = QCx * a0b1[ 8] + WQx * a0b1[11] + 0.5 * 1 * invZetaCD * (a0b0[ 2] - zetaAB * invZeta * a0b0[ 3]);
    a1b1[ 0] = PAz * a0b1[ 0] + WPz * a0b1[ 3] + 0.5 * 1 * invZeta * a0b0[ 1];
    a1b1[ 3] = PAy * a0b1[ 0] + WPy * a0b1[ 3];
    a1b1[ 5] = PAy * a0b1[ 2] + WPy * a0b1[ 5];
    a1b1[ 6] = PAx * a0b1[ 0] + WPx * a0b1[ 3];
    a1b1[ 7] = PAx * a0b1[ 1] + WPx * a0b1[ 4];
    a1b1[ 8] = PAx * a0b1[ 2] + WPx * a0b1[ 5] + 0.5 * 1 * invZeta * a0b0[ 1];
    a1b1[ 9] = PAz * a0b1[ 3] + WPz * a0b1[ 6] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[10] = PAz * a0b1[ 4] + WPz * a0b1[ 7];
    a1b1[11] = PAz * a0b1[ 5] + WPz * a0b1[ 8];
    a1b1[12] = PAy * a0b1[ 3] + WPy * a0b1[ 6];
    a1b1[13] = PAy * a0b1[ 4] + WPy * a0b1[ 7] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[14] = PAy * a0b1[ 5] + WPy * a0b1[ 8];
    a1b1[15] = PAx * a0b1[ 3] + WPx * a0b1[ 6];
    a1b1[16] = PAx * a0b1[ 4] + WPx * a0b1[ 7];
    a1b1[17] = PAx * a0b1[ 5] + WPx * a0b1[ 8] + 0.5 * 1 * invZeta * a0b0[ 2];
    a1b1[18] = PAz * a0b1[ 6] + WPz * a0b1[ 9] + 0.5 * 1 * invZeta * a0b0[ 3];
    a1b1[19] = PAz * a0b1[ 7] + WPz * a0b1[10];
    a1b1[20] = PAz * a0b1[ 8] + WPz * a0b1[11];
    a1b1[21] = PAy * a0b1[ 6] + WPy * a0b1[ 9];
    a1b1[22] = PAy * a0b1[ 7] + WPy * a0b1[10] + 0.5 * 1 * invZeta * a0b0[ 3];
    a1b1[23] = PAy * a0b1[ 8] + WPy * a0b1[11];
    a1b1[24] = PAx * a0b1[ 6] + WPx * a0b1[ 9];
    a1b1[25] = PAx * a0b1[ 7] + WPx * a0b1[10];
    a1b1[26] = PAx * a0b1[ 8] + WPx * a0b1[11] + 0.5 * 1 * invZeta * a0b0[ 3];
    a2b0[ 0] = PAz * a1b0[ 0] + WPz * a1b0[ 3] + 0.5 * 1 * invZetaAB * (a0b0[ 0] - zetaCD * invZeta * a0b0[ 1]);
    a2b0[ 1] = PAy * a1b0[ 0] + WPy * a1b0[ 3];
    a2b0[ 2] = PAy * a1b0[ 1] + WPy * a1b0[ 4] + 0.5 * 1 * invZetaAB * (a0b0[ 0] - zetaCD * invZeta * a0b0[ 1]);
    a2b0[ 3] = PAx * a1b0[ 0] + WPx * a1b0[ 3];
    a2b0[ 5] = PAx * a1b0[ 2] + WPx * a1b0[ 5] + 0.5 * 1 * invZetaAB * (a0b0[ 0] - zetaCD * invZeta * a0b0[ 1]);
    a2b0[ 6] = PAz * a1b0[ 3] + WPz * a1b0[ 6] + 0.5 * 1 * invZetaAB * (a0b0[ 1] - zetaCD * invZeta * a0b0[ 2]);
    a2b0[ 7] = PAy * a1b0[ 3] + WPy * a1b0[ 6];
    a2b0[ 8] = PAy * a1b0[ 4] + WPy * a1b0[ 7] + 0.5 * 1 * invZetaAB * (a0b0[ 1] - zetaCD * invZeta * a0b0[ 2]);
    a2b0[ 9] = PAx * a1b0[ 3] + WPx * a1b0[ 6];
    a2b0[11] = PAx * a1b0[ 5] + WPx * a1b0[ 8] + 0.5 * 1 * invZetaAB * (a0b0[ 1] - zetaCD * invZeta * a0b0[ 2]);
    a2b0[12] = PAz * a1b0[ 6] + WPz * a1b0[ 9] + 0.5 * 1 * invZetaAB * (a0b0[ 2] - zetaCD * invZeta * a0b0[ 3]);
    a2b0[13] = PAy * a1b0[ 6] + WPy * a1b0[ 9];
    a2b0[14] = PAy * a1b0[ 7] + WPy * a1b0[10] + 0.5 * 1 * invZetaAB * (a0b0[ 2] - zetaCD * invZeta * a0b0[ 3]);
    a2b0[15] = PAx * a1b0[ 6] + WPx * a1b0[ 9];
    a2b0[17] = PAx * a1b0[ 8] + WPx * a1b0[11] + 0.5 * 1 * invZetaAB * (a0b0[ 2] - zetaCD * invZeta * a0b0[ 3]);
    a1b2[ 5] = PAz * a0b2[ 5] + WPz * a0b2[11];
    a1b2[ 6] = PAy * a0b2[ 0] + WPy * a0b2[ 6];
    a1b2[ 7] = QCy * a1b1[ 3] + WQy * a1b1[12] + 0.5 * 1 * invZeta * a0b1[ 3];
    a1b2[ 8] = PAy * a0b2[ 2] + WPy * a0b2[ 8] + 0.5 * 2 * invZeta * a0b1[ 4];
    a1b2[11] = PAy * a0b2[ 5] + WPy * a0b2[11];
    a1b2[12] = PAx * a0b2[ 0] + WPx * a0b2[ 6];
    a1b2[14] = PAx * a0b2[ 2] + WPx * a0b2[ 8];
    a1b2[17] = PAx * a0b2[ 5] + WPx * a0b2[11] + 0.5 * 2 * invZeta * a0b1[ 5];
    a1b2[23] = PAz * a0b2[11] + WPz * a0b2[17];
    a1b2[24] = PAy * a0b2[ 6] + WPy * a0b2[12];
    a1b2[25] = PAy * a0b2[ 7] + WPy * a0b2[13] + 0.5 * 1 * invZeta * a0b1[ 6];
    a1b2[26] = QCy * a1b1[13] + WQy * a1b1[22] + 0.5 * 1 * invZetaCD * (a1b0[ 4] - zetaAB * invZeta * a1b0[ 7]) + 0.5 * 1 * invZeta * a0b1[ 7];
    a1b2[29] = PAy * a0b2[11] + WPy * a0b2[17];
    a1b2[30] = PAx * a0b2[ 6] + WPx * a0b2[12];
    a1b2[32] = QCy * a1b1[16] + WQy * a1b1[25] + 0.5 * 1 * invZetaCD * (a1b0[ 5] - zetaAB * invZeta * a1b0[ 8]);
    a1b2[35] = PAx * a0b2[11] + WPx * a0b2[17] + 0.5 * 2 * invZeta * a0b1[ 8];
    a2b1[ 0] = PAz * a1b1[ 0] + WPz * a1b1[ 9] + 0.5 * 1 * invZetaAB * (a0b1[ 0] - zetaCD * invZeta * a0b1[ 3]) + 0.5 * 1 * invZeta * a1b0[ 3];
    a2b1[ 1] = QCy * a2b0[ 0] + WQy * a2b0[ 6];
    a2b1[ 2] = QCx * a2b0[ 0] + WQx * a2b0[ 6];
    a2b1[ 5] = PAz * a1b1[ 5] + WPz * a1b1[14];
    a2b1[ 6] = PAy * a1b1[ 3] + WPy * a1b1[12] + 0.5 * 1 * invZetaAB * (a0b1[ 0] - zetaCD * invZeta * a0b1[ 3]);
    a2b1[ 8] = PAy * a1b1[ 5] + WPy * a1b1[14] + 0.5 * 1 * invZetaAB * (a0b1[ 2] - zetaCD * invZeta * a0b1[ 5]);
    a2b1[ 9] = PAx * a1b1[ 0] + WPx * a1b1[ 9];
    a2b1[10] = PAz * a1b1[ 7] + WPz * a1b1[16];
    a2b1[11] = PAz * a1b1[ 8] + WPz * a1b1[17];
    a2b1[12] = PAx * a1b1[ 3] + WPx * a1b1[12];
    a2b1[13] = PAy * a1b1[ 7] + WPy * a1b1[16] + 0.5 * 1 * invZeta * a1b0[ 5];
    a2b1[14] = PAx * a1b1[ 5] + WPx * a1b1[14] + 0.5 * 1 * invZeta * a1b0[ 4];
    a2b1[16] = PAx * a1b1[ 7] + WPx * a1b1[16] + 0.5 * 1 * invZetaAB * (a0b1[ 1] - zetaCD * invZeta * a0b1[ 4]);
    a2b1[17] = PAx * a1b1[ 8] + WPx * a1b1[17] + 0.5 * 1 * invZetaAB * (a0b1[ 2] - zetaCD * invZeta * a0b1[ 5]) + 0.5 * 1 * invZeta * a1b0[ 5];
    a2b1[18] = PAz * a1b1[ 9] + WPz * a1b1[18] + 0.5 * 1 * invZetaAB * (a0b1[ 3] - zetaCD * invZeta * a0b1[ 6]) + 0.5 * 1 * invZeta * a1b0[ 6];
    a2b1[19] = PAz * a1b1[10] + WPz * a1b1[19] + 0.5 * 1 * invZetaAB * (a0b1[ 4] - zetaCD * invZeta * a0b1[ 7]);
    a2b1[20] = PAz * a1b1[11] + WPz * a1b1[20] + 0.5 * 1 * invZetaAB * (a0b1[ 5] - zetaCD * invZeta * a0b1[ 8]);
    a2b1[23] = PAy * a1b1[11] + WPy * a1b1[20];
    a2b1[24] = PAy * a1b1[12] + WPy * a1b1[21] + 0.5 * 1 * invZetaAB * (a0b1[ 3] - zetaCD * invZeta * a0b1[ 6]);
    a2b1[26] = PAy * a1b1[14] + WPy * a1b1[23] + 0.5 * 1 * invZetaAB * (a0b1[ 5] - zetaCD * invZeta * a0b1[ 8]);
    a2b1[27] = PAx * a1b1[ 9] + WPx * a1b1[18];
    a2b1[28] = PAx * a1b1[10] + WPx * a1b1[19];
    a2b1[29] = PAx * a1b1[11] + WPx * a1b1[20] + 0.5 * 1 * invZeta * a1b0[ 6];
    a2b1[30] = PAx * a1b1[12] + WPx * a1b1[21];
    a2b1[31] = PAx * a1b1[13] + WPx * a1b1[22];
    a2b1[32] = PAx * a1b1[14] + WPx * a1b1[23] + 0.5 * 1 * invZeta * a1b0[ 7];
    a2b1[34] = PAx * a1b1[16] + WPx * a1b1[25] + 0.5 * 1 * invZetaAB * (a0b1[ 4] - zetaCD * invZeta * a0b1[ 7]);
    a2b1[35] = PAx * a1b1[17] + WPx * a1b1[26] + 0.5 * 1 * invZetaAB * (a0b1[ 5] - zetaCD * invZeta * a0b1[ 8]) + 0.5 * 1 * invZeta * a1b0[ 8];
    a2b2[ 0] = QCz * a2b1[ 0] + WQz * a2b1[18] + 0.5 * 1 * invZetaCD * (a2b0[ 0] - zetaAB * invZeta * a2b0[ 6]) + 0.5 * 2 * invZeta * a1b1[ 9];
    a2b2[ 1] = QCy * a2b1[ 0] + WQy * a2b1[18];
    a2b2[ 2] = QCy * a2b1[ 1] + WQy * a2b1[19] + 0.5 * 1 * invZetaCD * (a2b0[ 0] - zetaAB * invZeta * a2b0[ 6]);
    a2b2[ 3] = QCx * a2b1[ 0] + WQx * a2b1[18];
    a2b2[ 4] = QCx * a2b1[ 1] + WQx * a2b1[19];
    a2b2[ 5] = PAz * a1b2[ 5] + WPz * a1b2[23] + 0.5 * 1 * invZetaAB * (a0b2[ 5] - zetaCD * invZeta * a0b2[11]);
    a2b2[ 6] = PAz * a1b2[ 6] + WPz * a1b2[24] + 0.5 * 2 * invZeta * a1b1[12];
    a2b2[ 7] = PAz * a1b2[ 7] + WPz * a1b2[25] + 0.5 * 1 * invZeta * a1b1[13];
    a2b2[ 8] = PAz * a1b2[ 8] + WPz * a1b2[26];
    a2b2[ 9] = QCz * a2b1[ 5] + WQz * a2b1[23] + 0.5 * 1 * invZeta * a1b1[14];
    a2b2[10] = QCy * a2b1[ 5] + WQy * a2b1[23] + 0.5 * 1 * invZeta * a1b1[11];
    a2b2[11] = PAy * a1b2[ 5] + WPy * a1b2[23];
    a2b2[12] = PAy * a1b2[ 6] + WPy * a1b2[24] + 0.5 * 1 * invZetaAB * (a0b2[ 0] - zetaCD * invZeta * a0b2[ 6]);
    a2b2[13] = QCy * a2b1[ 6] + WQy * a2b1[24] + 0.5 * 2 * invZeta * a1b1[12];
    a2b2[14] = PAy * a1b2[ 8] + WPy * a1b2[26] + 0.5 * 1 * invZetaAB * (a0b2[ 2] - zetaCD * invZeta * a0b2[ 8]) + 0.5 * 2 * invZeta * a1b1[13];
    a2b2[15] = QCx * a2b1[ 6] + WQx * a2b1[24];
    a2b2[16] = QCy * a2b1[ 8] + WQy * a2b1[26] + 0.5 * 2 * invZeta * a1b1[14];
    a2b2[17] = PAy * a1b2[11] + WPy * a1b2[29] + 0.5 * 1 * invZetaAB * (a0b2[ 5] - zetaCD * invZeta * a0b2[11]);
    a2b2[18] = PAz * a1b2[12] + WPz * a1b2[30] + 0.5 * 2 * invZeta * a1b1[15];
    a2b2[19] = QCy * a2b1[ 9] + WQy * a2b1[27];
    a2b2[20] = PAz * a1b2[14] + WPz * a1b2[32];
    a2b2[21] = QCx * a2b1[ 9] + WQx * a2b1[27] + 0.5 * 1 * invZeta * a1b1[ 9];
    a2b2[22] = QCx * a2b1[10] + WQx * a2b1[28] + 0.5 * 1 * invZeta * a1b1[10];
    a2b2[23] = PAx * a1b2[ 5] + WPx * a1b2[23] + 0.5 * 2 * invZeta * a1b1[11];
    a2b2[24] = PAx * a1b2[ 6] + WPx * a1b2[24];
    a2b2[25] = PAx * a1b2[ 7] + WPx * a1b2[25];
    a2b2[26] = PAx * a1b2[ 8] + WPx * a1b2[26];
    a2b2[27] = QCx * a2b1[12] + WQx * a2b1[30] + 0.5 * 1 * invZeta * a1b1[12];
    a2b2[28] = QCx * a2b1[13] + WQx * a2b1[31] + 0.5 * 1 * invZeta * a1b1[13];
    a2b2[29] = PAx * a1b2[11] + WPx * a1b2[29] + 0.5 * 2 * invZeta * a1b1[14];
    a2b2[30] = PAx * a1b2[12] + WPx * a1b2[30] + 0.5 * 1 * invZetaAB * (a0b2[ 0] - zetaCD * invZeta * a0b2[ 6]);
    a2b2[31] = QCz * a2b1[16] + WQz * a2b1[34];
    a2b2[32] = PAx * a1b2[14] + WPx * a1b2[32] + 0.5 * 1 * invZetaAB * (a0b2[ 2] - zetaCD * invZeta * a0b2[ 8]);
    a2b2[33] = QCz * a2b1[17] + WQz * a2b1[35];
    a2b2[34] = QCx * a2b1[16] + WQx * a2b1[34] + 0.5 * 2 * invZeta * a1b1[16];
    a2b2[35] = PAx * a1b2[17] + WPx * a1b2[35] + 0.5 * 1 * invZetaAB * (a0b2[ 5] - zetaCD * invZeta * a0b2[11]) + 0.5 * 2 * invZeta * a1b1[17];
} // function (hgp_vrr_2_2_)

} // namespace (nhfInt)

