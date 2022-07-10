#ifndef GENOME_INDEX_H
#define GENOME_INDEX_H

#include "structures.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <vector>

namespace genomics {
  namespace {
      enum class bulge_state {none, dna, rna};

      struct affinity {
          uint64_t mismatches;
          uint64_t dna_bulges;
          uint64_t rna_bulges;
          bulge_state state;
          uint64_t curr_bulge_size;
      };
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  class genome_index {
  public:
    typedef sdsl::csa_wt<t_wt, t_dens, t_inv_dens> t_csa;
      
    t_csa csa;
    genome_structure gs;

    const char* search_alphabet = "ATCG";
    size_t search_alphabet_size = strlen(search_alphabet);

    genome_index() {}
    genome_index(t_csa csa, genome_structure gs) : csa(csa), gs(gs)
    {}
    genome_index(const genome_index& other) : csa(other.csa), gs(other.gs)
    {}

    friend void swap(genome_index& first, genome_index& second)
    {
      using std::swap;
      swap(first.csa, second.csa);
      swap(first.gs, second.gs);
    }

    genome_index& operator=(genome_index other)
    {
      swap(*this, other);
      return *this;
    }

    size_t resolve(size_t bwt_position) const {
      return csa[bwt_position];
    }

    /* 
       Searches for all strings in the genome matching the given
       query, up to a certain number of mismatches and allowing
       the wildcard character 'N' in the query.

       When a set of matches are found, the callback is called
       with the start and end position in the BWT of the genome,
       along with the number of mismatches, and a reference to
       the passed in data.

       To get the position of the matches in the original string,
       use the resolve(...) method of the class.
    */
    template <class t_data>
    void inexact_search(typename std::string::const_iterator begin,
                        typename std::string::const_iterator end,
                        size_t mismatches, 
                        const std::function<void(size_t, size_t, size_t, std::string, t_data&)> &callback,
                        t_data& data) const;

    template <class t_data>
    void inexact_search(typename std::string::const_iterator begin,
                        typename std::string::const_iterator end,
                        size_t sp, size_t ep, std::string match, size_t mismatches, size_t k,
                        const std::function<void(size_t, size_t, size_t, std::string, t_data&)> &callback,
                        t_data& data) const;

    /*
      PAM aware variant of inexact search where the PAM must match
      exactly on the right end. 
    */
    template <class t_data>
    void inexact_search(const std::string& query,
                        ssize_t position,
                        size_t sp, size_t ep,
                        std::string match,
                        const std::vector<std::string> &pams,
                        size_t mismatches, affinity aff,
                        const std::function<void(struct match, t_data&)> &callback,
                        t_data& data) const;

    template <class t_data>
    void inexact_search(const std::string& query,
                        const std::vector<std::string> &pams,
                        size_t mismatches, 
                        const std::function<void(struct match, t_data&)> &callback,
                        t_data& data) const;
  };

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::const_iterator begin,
                                                              typename std::string::const_iterator end,
                                                              size_t sp, size_t ep, std::string match,
                                                              size_t mismatches,
                                                              size_t k,
                                                              const std::function<void(size_t, size_t, size_t, std::string, t_data&)> &callback,
                                                              t_data& data) const {
    if (begin == end) {
      callback(sp, ep, k, match, data);
      return;
    }

    char c = *(end - 1);

    size_t occ_before = csa.rank_bwt(sp, c);
    size_t occ_within = csa.rank_bwt(ep + 1, c) - occ_before;

    if (occ_within > 0) {
      size_t sp_prime = csa.C[csa.char2comp[c]] + occ_before;
      size_t ep_prime = sp_prime + occ_within - 1;
      inexact_search(begin, end - 1, sp_prime, ep_prime, match + c, mismatches,
                     k, callback, data);
    }

    size_t cost = 1;
    if (k >= mismatches && c != 'N') return;
    if (c == 'N') cost = 0;

    for (size_t i = 0; i < search_alphabet_size; i++) {
      if (search_alphabet[i] == c) continue;

      char a = search_alphabet[i];

      occ_before = csa.rank_bwt(sp, a);
      occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

      if (occ_within > 0) {
        size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
        size_t ep_prime = sp_prime + occ_within - 1;
        inexact_search(begin, end - 1, sp_prime, ep_prime, match + a, mismatches,
                       k + cost, callback, data);
      }
    }
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::const_iterator begin,
                                                              typename std::string::const_iterator end,
                                                              size_t mismatches,
                                                              const std::function<void(size_t, size_t, size_t, std::string, t_data&)> &callback,
                                                              t_data& data) const {

    inexact_search(begin, end, 0, csa.size() - 1, "", mismatches, 0, callback, data);
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(const std::string& query,
                                                              ssize_t position,
                                                              size_t sp, size_t ep,
                                                              std::string sequence,
                                                              const std::vector<std::string> &pams,
                                                              size_t mismatches, 
                                                              affinity aff, 
                                                              const std::function<void(struct match, t_data&)> &callback,
                                                              t_data& data) const {
    if (position < 0) {
      std::function<void(size_t, size_t, size_t, std::string, t_data&)> matching_callback =
        [aff, callback](size_t sp, size_t ep, size_t mismatches, std::string sequence, t_data& data) {
          (void) mismatches;

          match m = {
              sequence, sp, ep, aff.mismatches, aff.dna_bulges, aff.rna_bulges
          };

          return callback(m, data);
      };

      for (const auto& pam : pams) {
        inexact_search(pam.begin(), pam.end(), sp, ep, sequence, 0, 0, matching_callback, data);
      }

      return;
    }

    char c = query[position];

    size_t occ_before = csa.rank_bwt(sp, c);
    size_t occ_within = csa.rank_bwt(ep + 1, c) - occ_before;

    if (occ_within > 0) {
      size_t sp_prime = csa.C[csa.char2comp[c]] + occ_before;
      size_t ep_prime = sp_prime + occ_within - 1;
      inexact_search(query, position - 1, sp_prime, ep_prime, sequence + c, pams,
                     mismatches, aff, callback, data);
    }

    size_t cost = 1;
    if (aff.mismatches >= mismatches) return;

    for (size_t i = 0; i < search_alphabet_size; i++) {
      if (search_alphabet[i] == c) continue;

      char a = search_alphabet[i];

      occ_before = csa.rank_bwt(sp, a);
      occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

      if (occ_within > 0) {
        size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
        size_t ep_prime = sp_prime + occ_within - 1;
        affinity aff_cost = aff;
        aff_cost.mismatches += cost;
        inexact_search(query, position - 1, sp_prime, ep_prime, sequence + a, pams,
                       mismatches, aff_cost, callback, data);
      }
    }
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(const std::string& query,
                                                              const std::vector<std::string> &pams,
                                                              size_t mismatches, 
                                                              const std::function<void(struct match, t_data&)> &callback,
                                                              t_data& data) const {
    affinity aff = {0, 0, 0, bulge_state::none, 0};
    inexact_search(query, query.length() - 1, 0, csa.size() - 1, "",
                   pams, mismatches, aff, callback, data);
  }

};

#endif /* GENOME_INDEX_H */
