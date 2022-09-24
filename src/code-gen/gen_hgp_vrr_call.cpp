#include <iostream>
#include <cstddef>

int main() {
    std::size_t eMax = 0;
    std::cin >> eMax;

    for (std::size_t e = 0; e <= eMax; ++e) {
        for (std::size_t f = 0; f <= eMax; ++f) {
            std::cout << "if (e = " << e << " && f = " << f << ") hgp_vrr_" << e << "_" << f << "_"
                << "(result, zetaAB, kAB, Px, Py, Pz, PAx, PAy, PAz, zetaCD, kCD, Qx, Qy, Qz, QCx, QCy, QCz);" << std::endl;
        }
    }

    return 0;
}
