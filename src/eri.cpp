#include "eri.hpp"
#include "basis.hpp"
#include "eri_class.hpp"
#include "hgp_eri_class.hpp"
#include <vector>
#include <cstddef>

namespace nhfInt {

// calculate eri using HGP method
std::vector<double> hgp_eri(const BasisSet &bs) {
    std::size_t nBs = bs.basis_size();
    std::size_t nInt = idx4(nBs-1, nBs-1, nBs-1, nBs-1) + 1;
    std::vector<double> eri(nInt);

    // calculate PairDataSet
    PairDataSet pds(bs);

    for (std::size_t i = 0; i < bs.size(); ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
        std::size_t ij = idx2(i,j);

        for (std::size_t k = 0; k < bs.size(); ++k) {
        for (std::size_t l = 0; l <= k; ++l) {
            std::size_t kl = idx2(k,l);

            if (ij > kl) {
                continue;
            }

            // adjust the order of basis functions
            // make sure that
            // a.angMom + b.angMom >= c.angMom + d.angMom
            // a.angMom >= b.angMom
            // c.angMom >= d.angMom
            std::size_t sum12 = bs[i].angMom + bs[j].angMom;
            std::size_t sum34 = bs[k].angMom + bs[l].angMom;
            bool norm00 = sum12 >= sum34;
            bool norm12 = bs[i].angMom >= bs[j].angMom;
            bool norm34 = bs[k].angMom >= bs[l].angMom;

            const Basis &a = (norm00 ? (norm12 ? bs[i] : bs[j]) : (norm34 ? bs[k] : bs[l]));
            const Basis &b = (norm00 ? (norm12 ? bs[j] : bs[i]) : (norm34 ? bs[l] : bs[k]));
            const Basis &c = (norm00 ? (norm34 ? bs[k] : bs[l]) : (norm12 ? bs[i] : bs[j]));
            const Basis &d = (norm00 ? (norm34 ? bs[l] : bs[k]) : (norm12 ? bs[j] : bs[i]));

            EriClass ec = hgp_eri_class(a, b, c, d, pds);

            for (std::size_t ai = 0; ai < ec.nbsA; ++ai) {
            for (std::size_t bj = 0; bj < ec.nbsA; ++bj) {
            for (std::size_t ck = 0; ck < ec.nbsA; ++ck) {
            for (std::size_t dl = 0; dl < ec.nbsA; ++dl) {
                // position of this eri in the global
                std::size_t iPos = idx4(a.bsBeg + ai, b.bsBeg + bj,
                                        c.bsBeg + ck, d.bsBeg + dl);
                eri[iPos] = ec(ai, bj, ck, dl);
            }}}}
        }}
    }}

    return eri;
}


}   // namespace (nhfInt)