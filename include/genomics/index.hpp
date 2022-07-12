#ifndef GENOME_INDEX_H
#define GENOME_INDEX_H

#include "structures.hpp"

#include <sdsl/suffix_arrays.hpp>
#include <vector>
#include <locale>

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
      PAM and bulge aware variant of inexact search where the PAM must match
      exactly on the right (or left) end. 
    */
    template <class t_data>
    void inexact_search(const std::string& query,
                        ssize_t position,
                        size_t sp, size_t ep,
                        std::string match,
                        const std::vector<std::string> &pams,
                        size_t mismatches, 
                        size_t max_rna_bulges,
                        size_t max_dna_bulges,
                        size_t max_bulge_size,
                        affinity aff,
                        const std::function<void(struct match, t_data&)> &callback,
                        t_data& data) const;

    template <class t_data>
    void inexact_search(const std::string& query,
                        const std::vector<std::string> &pams,
                        size_t mismatches, 
                        size_t max_rna_bulges,
                        size_t max_dna_bulges,
                        size_t max_bulge_size,
                        const std::function<void(struct match, t_data&)> &callback,
                        t_data& data) const;

    // bulge unaware variant (performance optimization)
    template <class t_data>
    void inexact_search(const std::string& query,
                        ssize_t position,
                        size_t sp, size_t ep,
                        std::string sequence,
                        const std::vector<std::string> &pams,
                        size_t mismatches, 
                        uint64_t k, 
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
                                                              uint64_t k, 
                                                              const std::function<void(struct match, t_data&)> &callback,
                                                              t_data& data) const {
    if (position < 0) {
      std::function<void(size_t, size_t, size_t, std::string, t_data&)> matching_callback =
        [k, callback](size_t sp, size_t ep, size_t mismatches, std::string sequence, t_data& data) {
          (void) mismatches;

          match m = {
            sequence, sp, ep, k, 0, 0
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
                     mismatches, k + 1, callback, data);
    }

    if (k >= mismatches) return;

    for (size_t i = 0; i < search_alphabet_size; i++) {
      if (search_alphabet[i] == c) continue;

      char a = search_alphabet[i];

      occ_before = csa.rank_bwt(sp, a);
      occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

      if (occ_within > 0) {
        size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
        size_t ep_prime = sp_prime + occ_within - 1;
        char a_lower = std::tolower(a); // to identify mismatches as lower case characters
        inexact_search(query, position - 1, sp_prime, ep_prime, sequence + a_lower, pams,
                       mismatches, k + 1, callback, data);
      }
    }
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(const std::string& query,
                                                              ssize_t position,
                                                              size_t sp, size_t ep,
                                                              std::string sequence,
                                                              const std::vector<std::string> &pams,
                                                              size_t mismatches,
                                                              size_t max_rna_bulges,
                                                              size_t max_dna_bulges,
                                                              size_t max_bulge_size,
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

      affinity aff_orig = aff; // make a copy as to note mutate original
      aff_orig.state = bulge_state::none;

      inexact_search(query, position - 1, sp_prime, ep_prime, sequence + c, pams,
                     mismatches, max_rna_bulges, max_dna_bulges, max_bulge_size, aff_orig,
                     callback, data);
    }

    if (mismatches > aff.mismatches) {
      for (size_t i = 0; i < search_alphabet_size; i++) {
        if (search_alphabet[i] == c) continue;

        char a = search_alphabet[i];

        occ_before = csa.rank_bwt(sp, a);
        occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

        if (occ_within > 0) {
          size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
          size_t ep_prime = sp_prime + occ_within - 1;

          affinity aff_mismatch = aff; // make a copy as to note mutate original
          aff_mismatch.state = bulge_state::none;
          aff_mismatch.mismatches += 1;

          char a_lower = std::tolower(a); // to identify mismatches as lower case characters
          inexact_search(query, position - 1, sp_prime, ep_prime, sequence + a_lower, pams,
                         mismatches, max_rna_bulges, max_dna_bulges, max_bulge_size,
                         aff_mismatch, callback, data);
        }
      }
    }

    affinity rna_bulge_aff = aff;
    if (max_rna_bulges > aff.rna_bulges) {
      if (aff.state != bulge_state::rna || rna_bulge_aff.curr_bulge_size == max_bulge_size) {
        rna_bulge_aff.state = bulge_state::rna;
        rna_bulge_aff.curr_bulge_size = 0;
        rna_bulge_aff.rna_bulges += 1;
      } 
    }

    if (rna_bulge_aff.state == bulge_state::rna && rna_bulge_aff.curr_bulge_size < max_bulge_size) {
      rna_bulge_aff.curr_bulge_size += 1;
      inexact_search(query, position - 1, sp, ep, sequence + '.', pams, mismatches, 
                     max_rna_bulges, max_dna_bulges, max_bulge_size,
                     rna_bulge_aff, callback, data);
    }

    affinity dna_bulge_aff = aff;
    if (max_dna_bulges > aff.dna_bulges) {
      if (aff.state != bulge_state::dna || dna_bulge_aff.curr_bulge_size == max_bulge_size) {
        dna_bulge_aff.state = bulge_state::dna;
        dna_bulge_aff.curr_bulge_size = 0;
        dna_bulge_aff.dna_bulges += 1;
      } 
    }

    if (dna_bulge_aff.state == bulge_state::dna && dna_bulge_aff.curr_bulge_size < max_bulge_size) {
      dna_bulge_aff.curr_bulge_size += 1;

      for (size_t i = 0; i < search_alphabet_size; i++) {
        char a = search_alphabet[i];

        occ_before = csa.rank_bwt(sp, a);
        occ_within = csa.rank_bwt(ep + 1, a) - occ_before;

        if (occ_within > 0) {
          size_t sp_prime = csa.C[csa.char2comp[a]] + occ_before;
          size_t ep_prime = sp_prime + occ_within - 1;

          char a_lower = std::tolower(a); // to identify mismatches as lower case characters
          inexact_search(query, position, sp_prime, ep_prime, sequence + a_lower, pams,
                         mismatches, max_rna_bulges, max_dna_bulges, max_bulge_size,
                         dna_bulge_aff, callback, data);
        }
      }
    }
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  template <class t_data>
  void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(const std::string& query,
                                                              const std::vector<std::string> &pams,
                                                              size_t mismatches, 
                                                              size_t max_rna_bulges,
                                                              size_t max_dna_bulges,
                                                              size_t max_bulge_size,
                                                              const std::function<void(struct match, t_data&)> &callback,
                                                              t_data& data) const {
    // performance optimization
    if (false && max_rna_bulges == 0 && max_dna_bulges == 0) {
      inexact_search(query, query.length() - 1, 0, csa.size() - 1, "",
                     pams, mismatches, 0, callback, data);
      return;
    }

    affinity aff = {0, 0, 0, bulge_state::none, 0};
    inexact_search(query, query.length() - 1, 0, csa.size() - 1, "",
                   pams, mismatches, max_rna_bulges, max_dna_bulges, max_bulge_size,
                   aff, callback, data);
  }

};

#endif /* GENOME_INDEX_H */
