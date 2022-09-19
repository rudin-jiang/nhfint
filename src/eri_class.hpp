#pragma once

#include <cstddef>
#include <vector>

namespace nhfInt {

// use a vector<double> to record eri class
class EriClass {
public:
    public:
    std::size_t angA, angB, angC, angD;
    std::size_t nbsA, nbsB, nbsC, nbsD;
    std::size_t num2, num3, num4;
    std::vector<double> eri;

public:
    EriClass(std::size_t a, std::size_t b, 
             std::size_t c, std::size_t d);

    EriClass(std::size_t a, std::size_t b, 
             std::size_t c, std::size_t d, 
             const std::vector<double> &eriVal);

    // view EriClass as four-tuple
    // i range from 0 to nbsA
    // j range from 0 to nbsB
    // k range from 0 to nbsC
    // l range from 0 to nbsD
    double  operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l) const;
    double& operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l);

    // view EriClass as two-tuple
    // i range from 0 to nbsA * nbsB
    // j range from 0 to nbsC * nbsD
    double  operator()(std::size_t i, std::size_t j) const;
    double& operator()(std::size_t i, std::size_t j);

    // get eri data directly
    double  operator[](std::size_t i) const;
    double& operator[](std::size_t i);

    std::size_t size() const
    { return nbsA * nbsB * nbsC * nbsD; }
};


} // namespace (nhfInt)