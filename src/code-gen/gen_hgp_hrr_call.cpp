#include <iostream>
#include <cstddef>

int main() {
    std::size_t aMax = 0;
    std::cin >> aMax;

    for (std::size_t a = 0; a <= aMax; ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            std::cout << "if (a = " << a << " && b = " << b << ") return hgp_hrr_" << a << "_" << b << "_"
                      << "(hrrInp, x, y, z);" << std::endl;
        }
    }

    return 0;
}