#include "hgp_hrr.hpp"
#include <vector>

namespace nhfInt {

using VecReal = std::vector<double>;

VecReal hgp_hrr_3_3_(const VecReal &hrrInp, double x, double y, double z) { 
    const double *a3b0 = hrrInp.data() +  0;
    const double *a4b0 = hrrInp.data() + 10;
    const double *a5b0 = hrrInp.data() + 25;
    const double *a6b0 = hrrInp.data() + 46;

    // intermediate arrays
    double a3b1[30];
    double a4b1[45];
    double a5b1[63];
    double a3b2[60];
    double a4b2[90];

    // resulted VecReal
    VecReal a3b3(100);

    a3b1[ 0] = a4b0[ 0] + z * a3b0[ 0];    // [(003)(001)] = [(004)(000)] + ABz * [(003)(000)]
    a3b1[ 1] = a4b0[ 1] + y * a3b0[ 0];    // [(003)(010)] = [(013)(000)] + ABy * [(003)(000)]
    a3b1[ 2] = a4b0[ 5] + x * a3b0[ 0];    // [(003)(100)] = [(103)(000)] + ABx * [(003)(000)]
    a3b1[ 3] = a4b0[ 1] + z * a3b0[ 1];    // [(012)(001)] = [(013)(000)] + ABz * [(012)(000)]
    a3b1[ 4] = a4b0[ 2] + y * a3b0[ 1];    // [(012)(010)] = [(022)(000)] + ABy * [(012)(000)]
    a3b1[ 5] = a4b0[ 6] + x * a3b0[ 1];    // [(012)(100)] = [(112)(000)] + ABx * [(012)(000)]
    a3b1[ 6] = a4b0[ 2] + z * a3b0[ 2];    // [(021)(001)] = [(022)(000)] + ABz * [(021)(000)]
    a3b1[ 7] = a4b0[ 3] + y * a3b0[ 2];    // [(021)(010)] = [(031)(000)] + ABy * [(021)(000)]
    a3b1[ 8] = a4b0[ 7] + x * a3b0[ 2];    // [(021)(100)] = [(121)(000)] + ABx * [(021)(000)]
    a3b1[ 9] = a4b0[ 3] + z * a3b0[ 3];    // [(030)(001)] = [(031)(000)] + ABz * [(030)(000)]
    a3b1[10] = a4b0[ 4] + y * a3b0[ 3];    // [(030)(010)] = [(040)(000)] + ABy * [(030)(000)]
    a3b1[11] = a4b0[ 8] + x * a3b0[ 3];    // [(030)(100)] = [(130)(000)] + ABx * [(030)(000)]
    a3b1[12] = a4b0[ 5] + z * a3b0[ 4];    // [(102)(001)] = [(103)(000)] + ABz * [(102)(000)]
    a3b1[13] = a4b0[ 6] + y * a3b0[ 4];    // [(102)(010)] = [(112)(000)] + ABy * [(102)(000)]
    a3b1[14] = a4b0[ 9] + x * a3b0[ 4];    // [(102)(100)] = [(202)(000)] + ABx * [(102)(000)]
    a3b1[15] = a4b0[ 6] + z * a3b0[ 5];    // [(111)(001)] = [(112)(000)] + ABz * [(111)(000)]
    a3b1[16] = a4b0[ 7] + y * a3b0[ 5];    // [(111)(010)] = [(121)(000)] + ABy * [(111)(000)]
    a3b1[17] = a4b0[10] + x * a3b0[ 5];    // [(111)(100)] = [(211)(000)] + ABx * [(111)(000)]
    a3b1[18] = a4b0[ 7] + z * a3b0[ 6];    // [(120)(001)] = [(121)(000)] + ABz * [(120)(000)]
    a3b1[19] = a4b0[ 8] + y * a3b0[ 6];    // [(120)(010)] = [(130)(000)] + ABy * [(120)(000)]
    a3b1[20] = a4b0[11] + x * a3b0[ 6];    // [(120)(100)] = [(220)(000)] + ABx * [(120)(000)]
    a3b1[21] = a4b0[ 9] + z * a3b0[ 7];    // [(201)(001)] = [(202)(000)] + ABz * [(201)(000)]
    a3b1[22] = a4b0[10] + y * a3b0[ 7];    // [(201)(010)] = [(211)(000)] + ABy * [(201)(000)]
    a3b1[23] = a4b0[12] + x * a3b0[ 7];    // [(201)(100)] = [(301)(000)] + ABx * [(201)(000)]
    a3b1[24] = a4b0[10] + z * a3b0[ 8];    // [(210)(001)] = [(211)(000)] + ABz * [(210)(000)]
    a3b1[25] = a4b0[11] + y * a3b0[ 8];    // [(210)(010)] = [(220)(000)] + ABy * [(210)(000)]
    a3b1[26] = a4b0[13] + x * a3b0[ 8];    // [(210)(100)] = [(310)(000)] + ABx * [(210)(000)]
    a3b1[27] = a4b0[12] + z * a3b0[ 9];    // [(300)(001)] = [(301)(000)] + ABz * [(300)(000)]
    a3b1[28] = a4b0[13] + y * a3b0[ 9];    // [(300)(010)] = [(310)(000)] + ABy * [(300)(000)]
    a3b1[29] = a4b0[14] + x * a3b0[ 9];    // [(300)(100)] = [(400)(000)] + ABx * [(300)(000)]
    a4b1[ 0] = a5b0[ 0] + z * a4b0[ 0];    // [(004)(001)] = [(005)(000)] + ABz * [(004)(000)]
    a4b1[ 2] = a5b0[ 6] + x * a4b0[ 0];    // [(004)(100)] = [(104)(000)] + ABx * [(004)(000)]
    a4b1[ 3] = a5b0[ 1] + z * a4b0[ 1];    // [(013)(001)] = [(014)(000)] + ABz * [(013)(000)]
    a4b1[ 4] = a5b0[ 2] + y * a4b0[ 1];    // [(013)(010)] = [(023)(000)] + ABy * [(013)(000)]
    a4b1[ 5] = a5b0[ 7] + x * a4b0[ 1];    // [(013)(100)] = [(113)(000)] + ABx * [(013)(000)]
    a4b1[ 6] = a5b0[ 2] + z * a4b0[ 2];    // [(022)(001)] = [(023)(000)] + ABz * [(022)(000)]
    a4b1[ 7] = a5b0[ 3] + y * a4b0[ 2];    // [(022)(010)] = [(032)(000)] + ABy * [(022)(000)]
    a4b1[ 8] = a5b0[ 8] + x * a4b0[ 2];    // [(022)(100)] = [(122)(000)] + ABx * [(022)(000)]
    a4b1[ 9] = a5b0[ 3] + z * a4b0[ 3];    // [(031)(001)] = [(032)(000)] + ABz * [(031)(000)]
    a4b1[10] = a5b0[ 4] + y * a4b0[ 3];    // [(031)(010)] = [(041)(000)] + ABy * [(031)(000)]
    a4b1[11] = a5b0[ 9] + x * a4b0[ 3];    // [(031)(100)] = [(131)(000)] + ABx * [(031)(000)]
    a4b1[12] = a5b0[ 4] + z * a4b0[ 4];    // [(040)(001)] = [(041)(000)] + ABz * [(040)(000)]
    a4b1[13] = a5b0[ 5] + y * a4b0[ 4];    // [(040)(010)] = [(050)(000)] + ABy * [(040)(000)]
    a4b1[14] = a5b0[10] + x * a4b0[ 4];    // [(040)(100)] = [(140)(000)] + ABx * [(040)(000)]
    a4b1[15] = a5b0[ 6] + z * a4b0[ 5];    // [(103)(001)] = [(104)(000)] + ABz * [(103)(000)]
    a4b1[16] = a5b0[ 7] + y * a4b0[ 5];    // [(103)(010)] = [(113)(000)] + ABy * [(103)(000)]
    a4b1[17] = a5b0[11] + x * a4b0[ 5];    // [(103)(100)] = [(203)(000)] + ABx * [(103)(000)]
    a4b1[18] = a5b0[ 7] + z * a4b0[ 6];    // [(112)(001)] = [(113)(000)] + ABz * [(112)(000)]
    a4b1[19] = a5b0[ 8] + y * a4b0[ 6];    // [(112)(010)] = [(122)(000)] + ABy * [(112)(000)]
    a4b1[20] = a5b0[12] + x * a4b0[ 6];    // [(112)(100)] = [(212)(000)] + ABx * [(112)(000)]
    a4b1[21] = a5b0[ 8] + z * a4b0[ 7];    // [(121)(001)] = [(122)(000)] + ABz * [(121)(000)]
    a4b1[22] = a5b0[ 9] + y * a4b0[ 7];    // [(121)(010)] = [(131)(000)] + ABy * [(121)(000)]
    a4b1[23] = a5b0[13] + x * a4b0[ 7];    // [(121)(100)] = [(221)(000)] + ABx * [(121)(000)]
    a4b1[24] = a5b0[ 9] + z * a4b0[ 8];    // [(130)(001)] = [(131)(000)] + ABz * [(130)(000)]
    a4b1[25] = a5b0[10] + y * a4b0[ 8];    // [(130)(010)] = [(140)(000)] + ABy * [(130)(000)]
    a4b1[26] = a5b0[14] + x * a4b0[ 8];    // [(130)(100)] = [(230)(000)] + ABx * [(130)(000)]
    a4b1[27] = a5b0[11] + z * a4b0[ 9];    // [(202)(001)] = [(203)(000)] + ABz * [(202)(000)]
    a4b1[28] = a5b0[12] + y * a4b0[ 9];    // [(202)(010)] = [(212)(000)] + ABy * [(202)(000)]
    a4b1[29] = a5b0[15] + x * a4b0[ 9];    // [(202)(100)] = [(302)(000)] + ABx * [(202)(000)]
    a4b1[30] = a5b0[12] + z * a4b0[10];    // [(211)(001)] = [(212)(000)] + ABz * [(211)(000)]
    a4b1[31] = a5b0[13] + y * a4b0[10];    // [(211)(010)] = [(221)(000)] + ABy * [(211)(000)]
    a4b1[32] = a5b0[16] + x * a4b0[10];    // [(211)(100)] = [(311)(000)] + ABx * [(211)(000)]
    a4b1[33] = a5b0[13] + z * a4b0[11];    // [(220)(001)] = [(221)(000)] + ABz * [(220)(000)]
    a4b1[34] = a5b0[14] + y * a4b0[11];    // [(220)(010)] = [(230)(000)] + ABy * [(220)(000)]
    a4b1[35] = a5b0[17] + x * a4b0[11];    // [(220)(100)] = [(320)(000)] + ABx * [(220)(000)]
    a4b1[36] = a5b0[15] + z * a4b0[12];    // [(301)(001)] = [(302)(000)] + ABz * [(301)(000)]
    a4b1[37] = a5b0[16] + y * a4b0[12];    // [(301)(010)] = [(311)(000)] + ABy * [(301)(000)]
    a4b1[38] = a5b0[18] + x * a4b0[12];    // [(301)(100)] = [(401)(000)] + ABx * [(301)(000)]
    a4b1[39] = a5b0[16] + z * a4b0[13];    // [(310)(001)] = [(311)(000)] + ABz * [(310)(000)]
    a4b1[40] = a5b0[17] + y * a4b0[13];    // [(310)(010)] = [(320)(000)] + ABy * [(310)(000)]
    a4b1[41] = a5b0[19] + x * a4b0[13];    // [(310)(100)] = [(410)(000)] + ABx * [(310)(000)]
    a4b1[42] = a5b0[18] + z * a4b0[14];    // [(400)(001)] = [(401)(000)] + ABz * [(400)(000)]
    a4b1[43] = a5b0[19] + y * a4b0[14];    // [(400)(010)] = [(410)(000)] + ABy * [(400)(000)]
    a4b1[44] = a5b0[20] + x * a4b0[14];    // [(400)(100)] = [(500)(000)] + ABx * [(400)(000)]
    a5b1[ 0] = a6b0[ 0] + z * a5b0[ 0];    // [(005)(001)] = [(006)(000)] + ABz * [(005)(000)]
    a5b1[ 2] = a6b0[ 7] + x * a5b0[ 0];    // [(005)(100)] = [(105)(000)] + ABx * [(005)(000)]
    a5b1[ 3] = a6b0[ 1] + z * a5b0[ 1];    // [(014)(001)] = [(015)(000)] + ABz * [(014)(000)]
    a5b1[ 6] = a6b0[ 2] + z * a5b0[ 2];    // [(023)(001)] = [(024)(000)] + ABz * [(023)(000)]
    a5b1[ 7] = a6b0[ 3] + y * a5b0[ 2];    // [(023)(010)] = [(033)(000)] + ABy * [(023)(000)]
    a5b1[ 9] = a6b0[ 3] + z * a5b0[ 3];    // [(032)(001)] = [(033)(000)] + ABz * [(032)(000)]
    a5b1[10] = a6b0[ 4] + y * a5b0[ 3];    // [(032)(010)] = [(042)(000)] + ABy * [(032)(000)]
    a5b1[12] = a6b0[ 4] + z * a5b0[ 4];    // [(041)(001)] = [(042)(000)] + ABz * [(041)(000)]
    a5b1[13] = a6b0[ 5] + y * a5b0[ 4];    // [(041)(010)] = [(051)(000)] + ABy * [(041)(000)]
    a5b1[14] = a6b0[11] + x * a5b0[ 4];    // [(041)(100)] = [(141)(000)] + ABx * [(041)(000)]
    a5b1[16] = a6b0[ 6] + y * a5b0[ 5];    // [(050)(010)] = [(060)(000)] + ABy * [(050)(000)]
    a5b1[18] = a6b0[ 7] + z * a5b0[ 6];    // [(104)(001)] = [(105)(000)] + ABz * [(104)(000)]
    a5b1[19] = a6b0[ 8] + y * a5b0[ 6];    // [(104)(010)] = [(114)(000)] + ABy * [(104)(000)]
    a5b1[20] = a6b0[13] + x * a5b0[ 6];    // [(104)(100)] = [(204)(000)] + ABx * [(104)(000)]
    a5b1[21] = a6b0[ 8] + z * a5b0[ 7];    // [(113)(001)] = [(114)(000)] + ABz * [(113)(000)]
    a5b1[22] = a6b0[ 9] + y * a5b0[ 7];    // [(113)(010)] = [(123)(000)] + ABy * [(113)(000)]
    a5b1[23] = a6b0[14] + x * a5b0[ 7];    // [(113)(100)] = [(213)(000)] + ABx * [(113)(000)]
    a5b1[24] = a6b0[ 9] + z * a5b0[ 8];    // [(122)(001)] = [(123)(000)] + ABz * [(122)(000)]
    a5b1[25] = a6b0[10] + y * a5b0[ 8];    // [(122)(010)] = [(132)(000)] + ABy * [(122)(000)]
    a5b1[26] = a6b0[15] + x * a5b0[ 8];    // [(122)(100)] = [(222)(000)] + ABx * [(122)(000)]
    a5b1[27] = a6b0[10] + z * a5b0[ 9];    // [(131)(001)] = [(132)(000)] + ABz * [(131)(000)]
    a5b1[28] = a6b0[11] + y * a5b0[ 9];    // [(131)(010)] = [(141)(000)] + ABy * [(131)(000)]
    a5b1[29] = a6b0[16] + x * a5b0[ 9];    // [(131)(100)] = [(231)(000)] + ABx * [(131)(000)]
    a5b1[31] = a6b0[12] + y * a5b0[10];    // [(140)(010)] = [(150)(000)] + ABy * [(140)(000)]
    a5b1[32] = a6b0[17] + x * a5b0[10];    // [(140)(100)] = [(240)(000)] + ABx * [(140)(000)]
    a5b1[33] = a6b0[13] + z * a5b0[11];    // [(203)(001)] = [(204)(000)] + ABz * [(203)(000)]
    a5b1[35] = a6b0[18] + x * a5b0[11];    // [(203)(100)] = [(303)(000)] + ABx * [(203)(000)]
    a5b1[36] = a6b0[14] + z * a5b0[12];    // [(212)(001)] = [(213)(000)] + ABz * [(212)(000)]
    a5b1[37] = a6b0[15] + y * a5b0[12];    // [(212)(010)] = [(222)(000)] + ABy * [(212)(000)]
    a5b1[38] = a6b0[19] + x * a5b0[12];    // [(212)(100)] = [(312)(000)] + ABx * [(212)(000)]
    a5b1[39] = a6b0[15] + z * a5b0[13];    // [(221)(001)] = [(222)(000)] + ABz * [(221)(000)]
    a5b1[40] = a6b0[16] + y * a5b0[13];    // [(221)(010)] = [(231)(000)] + ABy * [(221)(000)]
    a5b1[41] = a6b0[20] + x * a5b0[13];    // [(221)(100)] = [(321)(000)] + ABx * [(221)(000)]
    a5b1[43] = a6b0[17] + y * a5b0[14];    // [(230)(010)] = [(240)(000)] + ABy * [(230)(000)]
    a5b1[44] = a6b0[21] + x * a5b0[14];    // [(230)(100)] = [(330)(000)] + ABx * [(230)(000)]
    a5b1[45] = a6b0[18] + z * a5b0[15];    // [(302)(001)] = [(303)(000)] + ABz * [(302)(000)]
    a5b1[46] = a6b0[19] + y * a5b0[15];    // [(302)(010)] = [(312)(000)] + ABy * [(302)(000)]
    a5b1[47] = a6b0[22] + x * a5b0[15];    // [(302)(100)] = [(402)(000)] + ABx * [(302)(000)]
    a5b1[48] = a6b0[19] + z * a5b0[16];    // [(311)(001)] = [(312)(000)] + ABz * [(311)(000)]
    a5b1[49] = a6b0[20] + y * a5b0[16];    // [(311)(010)] = [(321)(000)] + ABy * [(311)(000)]
    a5b1[50] = a6b0[23] + x * a5b0[16];    // [(311)(100)] = [(411)(000)] + ABx * [(311)(000)]
    a5b1[51] = a6b0[20] + z * a5b0[17];    // [(320)(001)] = [(321)(000)] + ABz * [(320)(000)]
    a5b1[52] = a6b0[21] + y * a5b0[17];    // [(320)(010)] = [(330)(000)] + ABy * [(320)(000)]
    a5b1[53] = a6b0[24] + x * a5b0[17];    // [(320)(100)] = [(420)(000)] + ABx * [(320)(000)]
    a5b1[56] = a6b0[25] + x * a5b0[18];    // [(401)(100)] = [(501)(000)] + ABx * [(401)(000)]
    a5b1[57] = a6b0[23] + z * a5b0[19];    // [(410)(001)] = [(411)(000)] + ABz * [(410)(000)]
    a5b1[58] = a6b0[24] + y * a5b0[19];    // [(410)(010)] = [(420)(000)] + ABy * [(410)(000)]
    a5b1[59] = a6b0[26] + x * a5b0[19];    // [(410)(100)] = [(510)(000)] + ABx * [(410)(000)]
    a5b1[62] = a6b0[27] + x * a5b0[20];    // [(500)(100)] = [(600)(000)] + ABx * [(500)(000)]
    a3b2[ 0] = a4b1[ 0] + z * a3b1[ 0];    // [(003)(002)] = [(004)(001)] + ABz * [(003)(001)]
    a3b2[ 1] = a4b1[ 3] + y * a3b1[ 0];    // [(003)(011)] = [(013)(001)] + ABy * [(003)(001)]
    a3b2[ 2] = a4b1[ 4] + y * a3b1[ 1];    // [(003)(020)] = [(013)(010)] + ABy * [(003)(010)]
    a3b2[ 3] = a4b1[15] + x * a3b1[ 0];    // [(003)(101)] = [(103)(001)] + ABx * [(003)(001)]
    a3b2[ 5] = a4b1[17] + x * a3b1[ 2];    // [(003)(200)] = [(103)(100)] + ABx * [(003)(100)]
    a3b2[ 6] = a4b1[ 3] + z * a3b1[ 3];    // [(012)(002)] = [(013)(001)] + ABz * [(012)(001)]
    a3b2[ 7] = a4b1[ 6] + y * a3b1[ 3];    // [(012)(011)] = [(022)(001)] + ABy * [(012)(001)]
    a3b2[ 8] = a4b1[ 7] + y * a3b1[ 4];    // [(012)(020)] = [(022)(010)] + ABy * [(012)(010)]
    a3b2[ 9] = a4b1[18] + x * a3b1[ 3];    // [(012)(101)] = [(112)(001)] + ABx * [(012)(001)]
    a3b2[11] = a4b1[20] + x * a3b1[ 5];    // [(012)(200)] = [(112)(100)] + ABx * [(012)(100)]
    a3b2[12] = a4b1[ 6] + z * a3b1[ 6];    // [(021)(002)] = [(022)(001)] + ABz * [(021)(001)]
    a3b2[13] = a4b1[ 9] + y * a3b1[ 6];    // [(021)(011)] = [(031)(001)] + ABy * [(021)(001)]
    a3b2[14] = a4b1[10] + y * a3b1[ 7];    // [(021)(020)] = [(031)(010)] + ABy * [(021)(010)]
    a3b2[15] = a4b1[21] + x * a3b1[ 6];    // [(021)(101)] = [(121)(001)] + ABx * [(021)(001)]
    a3b2[17] = a4b1[23] + x * a3b1[ 8];    // [(021)(200)] = [(121)(100)] + ABx * [(021)(100)]
    a3b2[18] = a4b1[ 9] + z * a3b1[ 9];    // [(030)(002)] = [(031)(001)] + ABz * [(030)(001)]
    a3b2[20] = a4b1[13] + y * a3b1[10];    // [(030)(020)] = [(040)(010)] + ABy * [(030)(010)]
    a3b2[21] = a4b1[24] + x * a3b1[ 9];    // [(030)(101)] = [(130)(001)] + ABx * [(030)(001)]
    a3b2[23] = a4b1[26] + x * a3b1[11];    // [(030)(200)] = [(130)(100)] + ABx * [(030)(100)]
    a3b2[24] = a4b1[15] + z * a3b1[12];    // [(102)(002)] = [(103)(001)] + ABz * [(102)(001)]
    a3b2[25] = a4b1[18] + y * a3b1[12];    // [(102)(011)] = [(112)(001)] + ABy * [(102)(001)]
    a3b2[26] = a4b1[19] + y * a3b1[13];    // [(102)(020)] = [(112)(010)] + ABy * [(102)(010)]
    a3b2[27] = a4b1[27] + x * a3b1[12];    // [(102)(101)] = [(202)(001)] + ABx * [(102)(001)]
    a3b2[29] = a4b1[29] + x * a3b1[14];    // [(102)(200)] = [(202)(100)] + ABx * [(102)(100)]
    a3b2[30] = a4b1[18] + z * a3b1[15];    // [(111)(002)] = [(112)(001)] + ABz * [(111)(001)]
    a3b2[31] = a4b1[21] + y * a3b1[15];    // [(111)(011)] = [(121)(001)] + ABy * [(111)(001)]
    a3b2[32] = a4b1[22] + y * a3b1[16];    // [(111)(020)] = [(121)(010)] + ABy * [(111)(010)]
    a3b2[35] = a4b1[32] + x * a3b1[17];    // [(111)(200)] = [(211)(100)] + ABx * [(111)(100)]
    a3b2[36] = a4b1[21] + z * a3b1[18];    // [(120)(002)] = [(121)(001)] + ABz * [(120)(001)]
    a3b2[37] = a4b1[24] + y * a3b1[18];    // [(120)(011)] = [(130)(001)] + ABy * [(120)(001)]
    a3b2[38] = a4b1[25] + y * a3b1[19];    // [(120)(020)] = [(130)(010)] + ABy * [(120)(010)]
    a3b2[40] = a4b1[34] + x * a3b1[19];    // [(120)(110)] = [(220)(010)] + ABx * [(120)(010)]
    a3b2[41] = a4b1[35] + x * a3b1[20];    // [(120)(200)] = [(220)(100)] + ABx * [(120)(100)]
    a3b2[42] = a4b1[27] + z * a3b1[21];    // [(201)(002)] = [(202)(001)] + ABz * [(201)(001)]
    a3b2[43] = a4b1[30] + y * a3b1[21];    // [(201)(011)] = [(211)(001)] + ABy * [(201)(001)]
    a3b2[44] = a4b1[31] + y * a3b1[22];    // [(201)(020)] = [(211)(010)] + ABy * [(201)(010)]
    a3b2[46] = a4b1[37] + x * a3b1[22];    // [(201)(110)] = [(301)(010)] + ABx * [(201)(010)]
    a3b2[47] = a4b1[38] + x * a3b1[23];    // [(201)(200)] = [(301)(100)] + ABx * [(201)(100)]
    a3b2[48] = a4b1[30] + z * a3b1[24];    // [(210)(002)] = [(211)(001)] + ABz * [(210)(001)]
    a3b2[50] = a4b1[34] + y * a3b1[25];    // [(210)(020)] = [(220)(010)] + ABy * [(210)(010)]
    a3b2[52] = a4b1[40] + x * a3b1[25];    // [(210)(110)] = [(310)(010)] + ABx * [(210)(010)]
    a3b2[53] = a4b1[41] + x * a3b1[26];    // [(210)(200)] = [(310)(100)] + ABx * [(210)(100)]
    a3b2[54] = a4b1[36] + z * a3b1[27];    // [(300)(002)] = [(301)(001)] + ABz * [(300)(001)]
    a3b2[55] = a4b1[39] + y * a3b1[27];    // [(300)(011)] = [(310)(001)] + ABy * [(300)(001)]
    a3b2[56] = a4b1[40] + y * a3b1[28];    // [(300)(020)] = [(310)(010)] + ABy * [(300)(010)]
    a3b2[57] = a4b1[42] + x * a3b1[27];    // [(300)(101)] = [(400)(001)] + ABx * [(300)(001)]
    a3b2[59] = a4b1[44] + x * a3b1[29];    // [(300)(200)] = [(400)(100)] + ABx * [(300)(100)]
    a4b2[ 0] = a5b1[ 0] + z * a4b1[ 0];    // [(004)(002)] = [(005)(001)] + ABz * [(004)(001)]
    a4b2[ 3] = a5b1[18] + x * a4b1[ 0];    // [(004)(101)] = [(104)(001)] + ABx * [(004)(001)]
    a4b2[ 6] = a5b1[ 3] + z * a4b1[ 3];    // [(013)(002)] = [(014)(001)] + ABz * [(013)(001)]
    a4b2[ 7] = a5b1[ 6] + y * a4b1[ 3];    // [(013)(011)] = [(023)(001)] + ABy * [(013)(001)]
    a4b2[ 8] = a5b1[ 7] + y * a4b1[ 4];    // [(013)(020)] = [(023)(010)] + ABy * [(013)(010)]
    a4b2[11] = a5b1[23] + x * a4b1[ 5];    // [(013)(200)] = [(113)(100)] + ABx * [(013)(100)]
    a4b2[12] = a5b1[ 6] + z * a4b1[ 6];    // [(022)(002)] = [(023)(001)] + ABz * [(022)(001)]
    a4b2[13] = a5b1[ 9] + y * a4b1[ 6];    // [(022)(011)] = [(032)(001)] + ABy * [(022)(001)]
    a4b2[14] = a5b1[10] + y * a4b1[ 7];    // [(022)(020)] = [(032)(010)] + ABy * [(022)(010)]
    a4b2[17] = a5b1[26] + x * a4b1[ 8];    // [(022)(200)] = [(122)(100)] + ABx * [(022)(100)]
    a4b2[18] = a5b1[ 9] + z * a4b1[ 9];    // [(031)(002)] = [(032)(001)] + ABz * [(031)(001)]
    a4b2[20] = a5b1[13] + y * a4b1[10];    // [(031)(020)] = [(041)(010)] + ABy * [(031)(010)]
    a4b2[21] = a5b1[27] + x * a4b1[ 9];    // [(031)(101)] = [(131)(001)] + ABx * [(031)(001)]
    a4b2[23] = a5b1[29] + x * a4b1[11];    // [(031)(200)] = [(131)(100)] + ABx * [(031)(100)]
    a4b2[24] = a5b1[12] + z * a4b1[12];    // [(040)(002)] = [(041)(001)] + ABz * [(040)(001)]
    a4b2[26] = a5b1[16] + y * a4b1[13];    // [(040)(020)] = [(050)(010)] + ABy * [(040)(010)]
    a4b2[27] = a5b1[14] + z * a4b1[14];    // [(040)(101)] = [(041)(100)] + ABz * [(040)(100)]
    a4b2[29] = a5b1[32] + x * a4b1[14];    // [(040)(200)] = [(140)(100)] + ABx * [(040)(100)]
    a4b2[30] = a5b1[18] + z * a4b1[15];    // [(103)(002)] = [(104)(001)] + ABz * [(103)(001)]
    a4b2[31] = a5b1[21] + y * a4b1[15];    // [(103)(011)] = [(113)(001)] + ABy * [(103)(001)]
    a4b2[32] = a5b1[22] + y * a4b1[16];    // [(103)(020)] = [(113)(010)] + ABy * [(103)(010)]
    a4b2[33] = a5b1[33] + x * a4b1[15];    // [(103)(101)] = [(203)(001)] + ABx * [(103)(001)]
    a4b2[35] = a5b1[35] + x * a4b1[17];    // [(103)(200)] = [(203)(100)] + ABx * [(103)(100)]
    a4b2[36] = a5b1[21] + z * a4b1[18];    // [(112)(002)] = [(113)(001)] + ABz * [(112)(001)]
    a4b2[37] = a5b1[24] + y * a4b1[18];    // [(112)(011)] = [(122)(001)] + ABy * [(112)(001)]
    a4b2[38] = a5b1[25] + y * a4b1[19];    // [(112)(020)] = [(122)(010)] + ABy * [(112)(010)]
    a4b2[39] = a5b1[36] + x * a4b1[18];    // [(112)(101)] = [(212)(001)] + ABx * [(112)(001)]
    a4b2[41] = a5b1[38] + x * a4b1[20];    // [(112)(200)] = [(212)(100)] + ABx * [(112)(100)]
    a4b2[42] = a5b1[24] + z * a4b1[21];    // [(121)(002)] = [(122)(001)] + ABz * [(121)(001)]
    a4b2[43] = a5b1[27] + y * a4b1[21];    // [(121)(011)] = [(131)(001)] + ABy * [(121)(001)]
    a4b2[44] = a5b1[28] + y * a4b1[22];    // [(121)(020)] = [(131)(010)] + ABy * [(121)(010)]
    a4b2[46] = a5b1[40] + x * a4b1[22];    // [(121)(110)] = [(221)(010)] + ABx * [(121)(010)]
    a4b2[47] = a5b1[41] + x * a4b1[23];    // [(121)(200)] = [(221)(100)] + ABx * [(121)(100)]
    a4b2[48] = a5b1[27] + z * a4b1[24];    // [(130)(002)] = [(131)(001)] + ABz * [(130)(001)]
    a4b2[49] = a5b1[28] + z * a4b1[25];    // [(130)(011)] = [(131)(010)] + ABz * [(130)(010)]
    a4b2[50] = a5b1[31] + y * a4b1[25];    // [(130)(020)] = [(140)(010)] + ABy * [(130)(010)]
    a4b2[51] = a5b1[29] + z * a4b1[26];    // [(130)(101)] = [(131)(100)] + ABz * [(130)(100)]
    a4b2[53] = a5b1[44] + x * a4b1[26];    // [(130)(200)] = [(230)(100)] + ABx * [(130)(100)]
    a4b2[54] = a5b1[33] + z * a4b1[27];    // [(202)(002)] = [(203)(001)] + ABz * [(202)(001)]
    a4b2[56] = a5b1[37] + y * a4b1[28];    // [(202)(020)] = [(212)(010)] + ABy * [(202)(010)]
    a4b2[58] = a5b1[46] + x * a4b1[28];    // [(202)(110)] = [(302)(010)] + ABx * [(202)(010)]
    a4b2[59] = a5b1[47] + x * a4b1[29];    // [(202)(200)] = [(302)(100)] + ABx * [(202)(100)]
    a4b2[60] = a5b1[36] + z * a4b1[30];    // [(211)(002)] = [(212)(001)] + ABz * [(211)(001)]
    a4b2[61] = a5b1[39] + y * a4b1[30];    // [(211)(011)] = [(221)(001)] + ABy * [(211)(001)]
    a4b2[62] = a5b1[40] + y * a4b1[31];    // [(211)(020)] = [(221)(010)] + ABy * [(211)(010)]
    a4b2[64] = a5b1[49] + x * a4b1[31];    // [(211)(110)] = [(311)(010)] + ABx * [(211)(010)]
    a4b2[65] = a5b1[50] + x * a4b1[32];    // [(211)(200)] = [(311)(100)] + ABx * [(211)(100)]
    a4b2[66] = a5b1[39] + z * a4b1[33];    // [(220)(002)] = [(221)(001)] + ABz * [(220)(001)]
    a4b2[68] = a5b1[43] + y * a4b1[34];    // [(220)(020)] = [(230)(010)] + ABy * [(220)(010)]
    a4b2[71] = a5b1[53] + x * a4b1[35];    // [(220)(200)] = [(320)(100)] + ABx * [(220)(100)]
    a4b2[72] = a5b1[45] + z * a4b1[36];    // [(301)(002)] = [(302)(001)] + ABz * [(301)(001)]
    a4b2[75] = a5b1[47] + z * a4b1[38];    // [(301)(101)] = [(302)(100)] + ABz * [(301)(100)]
    a4b2[77] = a5b1[56] + x * a4b1[38];    // [(301)(200)] = [(401)(100)] + ABx * [(301)(100)]
    a4b2[78] = a5b1[48] + z * a4b1[39];    // [(310)(002)] = [(311)(001)] + ABz * [(310)(001)]
    a4b2[79] = a5b1[51] + y * a4b1[39];    // [(310)(011)] = [(320)(001)] + ABy * [(310)(001)]
    a4b2[80] = a5b1[52] + y * a4b1[40];    // [(310)(020)] = [(320)(010)] + ABy * [(310)(010)]
    a4b2[81] = a5b1[57] + x * a4b1[39];    // [(310)(101)] = [(410)(001)] + ABx * [(310)(001)]
    a4b2[83] = a5b1[59] + x * a4b1[41];    // [(310)(200)] = [(410)(100)] + ABx * [(310)(100)]
    a4b2[86] = a5b1[58] + y * a4b1[43];    // [(400)(020)] = [(410)(010)] + ABy * [(400)(010)]
    a4b2[89] = a5b1[62] + x * a4b1[44];    // [(400)(200)] = [(500)(100)] + ABx * [(400)(100)]
    a3b3[ 0] = a4b2[ 0] + z * a3b2[ 0];    // [(003)(003)] = [(004)(002)] + ABz * [(003)(002)]
    a3b3[ 1] = a4b2[ 6] + y * a3b2[ 0];    // [(003)(012)] = [(013)(002)] + ABy * [(003)(002)]
    a3b3[ 2] = a4b2[ 7] + y * a3b2[ 1];    // [(003)(021)] = [(013)(011)] + ABy * [(003)(011)]
    a3b3[ 3] = a4b2[ 8] + y * a3b2[ 2];    // [(003)(030)] = [(013)(020)] + ABy * [(003)(020)]
    a3b3[ 4] = a4b2[30] + x * a3b2[ 0];    // [(003)(102)] = [(103)(002)] + ABx * [(003)(002)]
    a3b3[ 5] = a4b2[31] + x * a3b2[ 1];    // [(003)(111)] = [(103)(011)] + ABx * [(003)(011)]
    a3b3[ 6] = a4b2[32] + x * a3b2[ 2];    // [(003)(120)] = [(103)(020)] + ABx * [(003)(020)]
    a3b3[ 7] = a4b2[33] + x * a3b2[ 3];    // [(003)(201)] = [(103)(101)] + ABx * [(003)(101)]
    a3b3[ 8] = a4b2[11] + y * a3b2[ 5];    // [(003)(210)] = [(013)(200)] + ABy * [(003)(200)]
    a3b3[ 9] = a4b2[35] + x * a3b2[ 5];    // [(003)(300)] = [(103)(200)] + ABx * [(003)(200)]
    a3b3[10] = a4b2[ 6] + z * a3b2[ 6];    // [(012)(003)] = [(013)(002)] + ABz * [(012)(002)]
    a3b3[11] = a4b2[12] + y * a3b2[ 6];    // [(012)(012)] = [(022)(002)] + ABy * [(012)(002)]
    a3b3[12] = a4b2[13] + y * a3b2[ 7];    // [(012)(021)] = [(022)(011)] + ABy * [(012)(011)]
    a3b3[13] = a4b2[14] + y * a3b2[ 8];    // [(012)(030)] = [(022)(020)] + ABy * [(012)(020)]
    a3b3[14] = a4b2[36] + x * a3b2[ 6];    // [(012)(102)] = [(112)(002)] + ABx * [(012)(002)]
    a3b3[15] = a4b2[37] + x * a3b2[ 7];    // [(012)(111)] = [(112)(011)] + ABx * [(012)(011)]
    a3b3[16] = a4b2[38] + x * a3b2[ 8];    // [(012)(120)] = [(112)(020)] + ABx * [(012)(020)]
    a3b3[17] = a4b2[39] + x * a3b2[ 9];    // [(012)(201)] = [(112)(101)] + ABx * [(012)(101)]
    a3b3[18] = a4b2[17] + y * a3b2[11];    // [(012)(210)] = [(022)(200)] + ABy * [(012)(200)]
    a3b3[19] = a4b2[41] + x * a3b2[11];    // [(012)(300)] = [(112)(200)] + ABx * [(012)(200)]
    a3b3[20] = a4b2[12] + z * a3b2[12];    // [(021)(003)] = [(022)(002)] + ABz * [(021)(002)]
    a3b3[21] = a4b2[18] + y * a3b2[12];    // [(021)(012)] = [(031)(002)] + ABy * [(021)(002)]
    a3b3[22] = a4b2[14] + z * a3b2[14];    // [(021)(021)] = [(022)(020)] + ABz * [(021)(020)]
    a3b3[23] = a4b2[20] + y * a3b2[14];    // [(021)(030)] = [(031)(020)] + ABy * [(021)(020)]
    a3b3[24] = a4b2[42] + x * a3b2[12];    // [(021)(102)] = [(121)(002)] + ABx * [(021)(002)]
    a3b3[25] = a4b2[43] + x * a3b2[13];    // [(021)(111)] = [(121)(011)] + ABx * [(021)(011)]
    a3b3[26] = a4b2[44] + x * a3b2[14];    // [(021)(120)] = [(121)(020)] + ABx * [(021)(020)]
    a3b3[27] = a4b2[17] + z * a3b2[17];    // [(021)(201)] = [(022)(200)] + ABz * [(021)(200)]
    a3b3[28] = a4b2[46] + x * a3b2[16];    // [(021)(210)] = [(121)(110)] + ABx * [(021)(110)]
    a3b3[29] = a4b2[47] + x * a3b2[17];    // [(021)(300)] = [(121)(200)] + ABx * [(021)(200)]
    a3b3[30] = a4b2[18] + z * a3b2[18];    // [(030)(003)] = [(031)(002)] + ABz * [(030)(002)]
    a3b3[31] = a4b2[24] + y * a3b2[18];    // [(030)(012)] = [(040)(002)] + ABy * [(030)(002)]
    a3b3[32] = a4b2[20] + z * a3b2[20];    // [(030)(021)] = [(031)(020)] + ABz * [(030)(020)]
    a3b3[33] = a4b2[26] + y * a3b2[20];    // [(030)(030)] = [(040)(020)] + ABy * [(030)(020)]
    a3b3[34] = a4b2[48] + x * a3b2[18];    // [(030)(102)] = [(130)(002)] + ABx * [(030)(002)]
    a3b3[35] = a4b2[49] + x * a3b2[19];    // [(030)(111)] = [(130)(011)] + ABx * [(030)(011)]
    a3b3[36] = a4b2[50] + x * a3b2[20];    // [(030)(120)] = [(130)(020)] + ABx * [(030)(020)]
    a3b3[37] = a4b2[51] + x * a3b2[21];    // [(030)(201)] = [(130)(101)] + ABx * [(030)(101)]
    a3b3[38] = a4b2[29] + y * a3b2[23];    // [(030)(210)] = [(040)(200)] + ABy * [(030)(200)]
    a3b3[39] = a4b2[53] + x * a3b2[23];    // [(030)(300)] = [(130)(200)] + ABx * [(030)(200)]
    a3b3[40] = a4b2[30] + z * a3b2[24];    // [(102)(003)] = [(103)(002)] + ABz * [(102)(002)]
    a3b3[41] = a4b2[36] + y * a3b2[24];    // [(102)(012)] = [(112)(002)] + ABy * [(102)(002)]
    a3b3[42] = a4b2[37] + y * a3b2[25];    // [(102)(021)] = [(112)(011)] + ABy * [(102)(011)]
    a3b3[43] = a4b2[38] + y * a3b2[26];    // [(102)(030)] = [(112)(020)] + ABy * [(102)(020)]
    a3b3[44] = a4b2[54] + x * a3b2[24];    // [(102)(102)] = [(202)(002)] + ABx * [(102)(002)]
    a3b3[45] = a4b2[39] + y * a3b2[27];    // [(102)(111)] = [(112)(101)] + ABy * [(102)(101)]
    a3b3[46] = a4b2[56] + x * a3b2[26];    // [(102)(120)] = [(202)(020)] + ABx * [(102)(020)]
    a3b3[47] = a4b2[35] + z * a3b2[29];    // [(102)(201)] = [(103)(200)] + ABz * [(102)(200)]
    a3b3[48] = a4b2[58] + x * a3b2[28];    // [(102)(210)] = [(202)(110)] + ABx * [(102)(110)]
    a3b3[49] = a4b2[59] + x * a3b2[29];    // [(102)(300)] = [(202)(200)] + ABx * [(102)(200)]
    a3b3[50] = a4b2[36] + z * a3b2[30];    // [(111)(003)] = [(112)(002)] + ABz * [(111)(002)]
    a3b3[51] = a4b2[42] + y * a3b2[30];    // [(111)(012)] = [(121)(002)] + ABy * [(111)(002)]
    a3b3[52] = a4b2[43] + y * a3b2[31];    // [(111)(021)] = [(121)(011)] + ABy * [(111)(011)]
    a3b3[53] = a4b2[44] + y * a3b2[32];    // [(111)(030)] = [(121)(020)] + ABy * [(111)(020)]
    a3b3[54] = a4b2[60] + x * a3b2[30];    // [(111)(102)] = [(211)(002)] + ABx * [(111)(002)]
    a3b3[55] = a4b2[61] + x * a3b2[31];    // [(111)(111)] = [(211)(011)] + ABx * [(111)(011)]
    a3b3[56] = a4b2[62] + x * a3b2[32];    // [(111)(120)] = [(211)(020)] + ABx * [(111)(020)]
    a3b3[57] = a4b2[41] + z * a3b2[35];    // [(111)(201)] = [(112)(200)] + ABz * [(111)(200)]
    a3b3[58] = a4b2[64] + x * a3b2[34];    // [(111)(210)] = [(211)(110)] + ABx * [(111)(110)]
    a3b3[59] = a4b2[65] + x * a3b2[35];    // [(111)(300)] = [(211)(200)] + ABx * [(111)(200)]
    a3b3[60] = a4b2[42] + z * a3b2[36];    // [(120)(003)] = [(121)(002)] + ABz * [(120)(002)]
    a3b3[61] = a4b2[48] + y * a3b2[36];    // [(120)(012)] = [(130)(002)] + ABy * [(120)(002)]
    a3b3[62] = a4b2[49] + y * a3b2[37];    // [(120)(021)] = [(130)(011)] + ABy * [(120)(011)]
    a3b3[63] = a4b2[50] + y * a3b2[38];    // [(120)(030)] = [(130)(020)] + ABy * [(120)(020)]
    a3b3[64] = a4b2[66] + x * a3b2[36];    // [(120)(102)] = [(220)(002)] + ABx * [(120)(002)]
    a3b3[65] = a4b2[51] + y * a3b2[39];    // [(120)(111)] = [(130)(101)] + ABy * [(120)(101)]
    a3b3[66] = a4b2[68] + x * a3b2[38];    // [(120)(120)] = [(220)(020)] + ABx * [(120)(020)]
    a3b3[67] = a4b2[47] + z * a3b2[41];    // [(120)(201)] = [(121)(200)] + ABz * [(120)(200)]
    a3b3[68] = a4b2[53] + y * a3b2[41];    // [(120)(210)] = [(130)(200)] + ABy * [(120)(200)]
    a3b3[69] = a4b2[71] + x * a3b2[41];    // [(120)(300)] = [(220)(200)] + ABx * [(120)(200)]
    a3b3[70] = a4b2[54] + z * a3b2[42];    // [(201)(003)] = [(202)(002)] + ABz * [(201)(002)]
    a3b3[71] = a4b2[60] + y * a3b2[42];    // [(201)(012)] = [(211)(002)] + ABy * [(201)(002)]
    a3b3[72] = a4b2[61] + y * a3b2[43];    // [(201)(021)] = [(211)(011)] + ABy * [(201)(011)]
    a3b3[73] = a4b2[62] + y * a3b2[44];    // [(201)(030)] = [(211)(020)] + ABy * [(201)(020)]
    a3b3[74] = a4b2[72] + x * a3b2[42];    // [(201)(102)] = [(301)(002)] + ABx * [(201)(002)]
    a3b3[75] = a4b2[58] + z * a3b2[46];    // [(201)(111)] = [(202)(110)] + ABz * [(201)(110)]
    a3b3[76] = a4b2[64] + y * a3b2[46];    // [(201)(120)] = [(211)(110)] + ABy * [(201)(110)]
    a3b3[77] = a4b2[75] + x * a3b2[45];    // [(201)(201)] = [(301)(101)] + ABx * [(201)(101)]
    a3b3[78] = a4b2[65] + y * a3b2[47];    // [(201)(210)] = [(211)(200)] + ABy * [(201)(200)]
    a3b3[79] = a4b2[77] + x * a3b2[47];    // [(201)(300)] = [(301)(200)] + ABx * [(201)(200)]
    a3b3[80] = a4b2[60] + z * a3b2[48];    // [(210)(003)] = [(211)(002)] + ABz * [(210)(002)]
    a3b3[81] = a4b2[66] + y * a3b2[48];    // [(210)(012)] = [(220)(002)] + ABy * [(210)(002)]
    a3b3[82] = a4b2[62] + z * a3b2[50];    // [(210)(021)] = [(211)(020)] + ABz * [(210)(020)]
    a3b3[83] = a4b2[68] + y * a3b2[50];    // [(210)(030)] = [(220)(020)] + ABy * [(210)(020)]
    a3b3[84] = a4b2[78] + x * a3b2[48];    // [(210)(102)] = [(310)(002)] + ABx * [(210)(002)]
    a3b3[85] = a4b2[79] + x * a3b2[49];    // [(210)(111)] = [(310)(011)] + ABx * [(210)(011)]
    a3b3[86] = a4b2[80] + x * a3b2[50];    // [(210)(120)] = [(310)(020)] + ABx * [(210)(020)]
    a3b3[87] = a4b2[81] + x * a3b2[51];    // [(210)(201)] = [(310)(101)] + ABx * [(210)(101)]
    a3b3[88] = a4b2[71] + y * a3b2[53];    // [(210)(210)] = [(220)(200)] + ABy * [(210)(200)]
    a3b3[89] = a4b2[83] + x * a3b2[53];    // [(210)(300)] = [(310)(200)] + ABx * [(210)(200)]
    a3b3[90] = a4b2[72] + z * a3b2[54];    // [(300)(003)] = [(301)(002)] + ABz * [(300)(002)]
    a3b3[91] = a4b2[78] + y * a3b2[54];    // [(300)(012)] = [(310)(002)] + ABy * [(300)(002)]
    a3b3[92] = a4b2[79] + y * a3b2[55];    // [(300)(021)] = [(310)(011)] + ABy * [(300)(011)]
    a3b3[93] = a4b2[80] + y * a3b2[56];    // [(300)(030)] = [(310)(020)] + ABy * [(300)(020)]
    a3b3[94] = a4b2[75] + z * a3b2[57];    // [(300)(102)] = [(301)(101)] + ABz * [(300)(101)]
    a3b3[95] = a4b2[81] + y * a3b2[57];    // [(300)(111)] = [(310)(101)] + ABy * [(300)(101)]
    a3b3[96] = a4b2[86] + x * a3b2[56];    // [(300)(120)] = [(400)(020)] + ABx * [(300)(020)]
    a3b3[97] = a4b2[77] + z * a3b2[59];    // [(300)(201)] = [(301)(200)] + ABz * [(300)(200)]
    a3b3[98] = a4b2[83] + y * a3b2[59];    // [(300)(210)] = [(310)(200)] + ABy * [(300)(200)]
    a3b3[99] = a4b2[89] + x * a3b2[59];    // [(300)(300)] = [(400)(200)] + ABx * [(300)(200)]

    return a3b3;

}  // function (hgp_hrr_3_3_)

}  // namespace (nhfInt)

