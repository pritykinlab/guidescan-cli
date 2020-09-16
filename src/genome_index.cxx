
#include "genome_index.hpp"

namespace genomics {
    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    template <class t_data>
    void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::iterator begin,
                                                                typename std::string::iterator end,
                                                                size_t sp, size_t ep, size_t mismatches,
                                                                const std::function<void(size_t, size_t, t_data&)> callback,
                                                                t_data& data) {
        if (begin == end) {
            callback(sp, ep, data);
        }

        char c = *(end - 1);

        size_t occ_before = csa.rank_bwt(sp, c);
        size_t occ_within = csa.rank_bwt(ep + 1, c) - occ_before;

        if (occ_within > 0) {
            size_t sp_prime = csa.C[csa.char2comp[c]] + occ_before;
            size_t ep_prime = sp_prime + occ_within - 1;
            inexact_search(begin, end - 1, sp_prime, ep_prime, mismatches,
                           callback, data);
        }

        if (mismatches < 1) return;

        for (size_t i = 0; i < search_alphabet_size; i++) {
            if (search_alphabet[i] == c) continue;

            char a = search_alphabet[i];

            occ_before = csa.rank_bwt(sp, a);
            occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

            if (occ_within > 0) {
                size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
                size_t ep_prime = sp_prime + occ_within - 1;
                inexact_search(begin, end - 1, sp_prime, ep_prime, mismatches - 1,
                               callback, data);
            }
        }
    }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    template <class t_data>
    void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::iterator begin,
                                                                typename std::string::iterator end,
                                                                size_t mismatches,
                                                                const std::function<void(size_t, size_t, t_data&)> callback,
                                                                t_data& data) {

        inexact_search(begin, end, 0, csa.size() - 1, mismatches, callback, data);
    }
};
