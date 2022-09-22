#include "eri_class.hpp"
#include <cstddef>
#include <vector>
#include <cassert>

namespace nhfInt {

EriClass::EriClass()
: angA(0), angB(0), angC(0), angD(0),
  nbsA(0), nbsB(0), nbsC(0), nbsD(0)
{
    num4 = nbsD;
    num3 = nbsC * num4;
    num2 = nbsB * num3;

    eri = std::vector<double>(size());
}

EriClass::EriClass(
    std::size_t a, std::size_t b, 
    std::size_t c, std::size_t d )
: angA(a), angB(b), angC(c), angD(d) {
    nbsA = (a + 1) * (a + 2) / 2;
    nbsB = (b + 1) * (b + 2) / 2;
    nbsC = (c + 1) * (c + 2) / 2;
    nbsD = (d + 1) * (d + 2) / 2;

    num4 = nbsD;
    num3 = nbsC * num4;
    num2 = nbsB * num3;

    eri = std::vector<double>(size());
}

EriClass::EriClass(
    std::size_t a, std::size_t b,
    std::size_t c, std::size_t d, 
    const std::vector<double> &eriVal)
: angA(a), angB(b), angC(c), angD(d) {
    nbsA = (a + 1) * (a + 2) / 2;
    nbsB = (b + 1) * (b + 2) / 2;
    nbsC = (c + 1) * (c + 2) / 2;
    nbsD = (d + 1) * (d + 2) / 2;

    num4 = nbsD;
    num3 = nbsC * num4;
    num2 = nbsB * num3;

    assert(eriVal.size() == size());
    eri = eriVal;
}

double  EriClass::operator()(
    std::size_t i, std::size_t j,
    std::size_t k, std::size_t l) const
{ return eri[i * num2 + j * num3 + k * num4 + l]; }

double& EriClass::operator()(
    std::size_t i, std::size_t j,
    std::size_t k, std::size_t l)
{ return eri[i * num2 + j * num3 + k * num4 + l]; }

double  EriClass::operator()(std::size_t i, std::size_t j) const
{ return eri[i * num3 + j]; }

double& EriClass::operator()(std::size_t i, std::size_t j)
{ return eri[i * num3 + j]; }

double  EriClass::operator[](std::size_t i) const
{ return eri[i]; }

double& EriClass::operator[](std::size_t i)
{ return eri[i]; }


}  // namespace (nhfInt)