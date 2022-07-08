#ifndef PROCESS_H
#define PROCESS_H

#include <set>
#include <list>
#include <tuple>

#include "json.hpp"
#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/printer.hpp"

namespace genomics {
  namespace {
    void off_target_enumerator(size_t sp, size_t ep, size_t k, std::string match,
                               std::vector<std::set<std::tuple<size_t, size_t, std::string>>> &off_targets_bwt) {
      off_targets_bwt[k].insert(std::make_tuple(sp, ep, match));
    }

    std::function<void(size_t, size_t, size_t, std::string, std::vector<std::set<std::tuple<size_t, size_t, std::string>>>&)> callback = off_target_enumerator;

    void off_target_counter(size_t sp, size_t ep, size_t k, std::string match, size_t &count) {
      (void) k;
      (void) match;
      count += ep - sp + 1;
    }

    std::function<void(size_t, size_t, size_t, std::string, size_t&)> counting_callback = off_target_counter;
  };

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi_forward,
                              const genome_index<t_wt, t_dens, t_inv_dens>& gi_reverse,
                              std::vector<std::string> pams, size_t mismatches,
                              int threshold, bool start, std::string out_format,
                              const kmer& k,
                              std::ostream& output,
                              std::mutex& output_mtx) {
    size_t count = 0;

    /* Because of the way inexact searching is implemented (from
     * back-to-front) I search for the reverse complement of the
     * kmer on the reverse complement strand (which is essentially
     * searching forward) and I search for the reverse complement
     * of the kmer on the forward strand.
     */
        
    if (k.pam == "") {
      pams = {std::string("")};
    } else {
      pams.push_back(k.pam);
    }

    std::vector<std::string> pams_c;
    for (const auto& pam : pams) {
      pams_c.push_back(genomics::reverse_complement(pam));
    }
        
    std::string kmer = !start ? genomics::reverse_complement(k.sequence) : k.sequence;

    if (threshold > 0 && !start) {
      gi_forward.inexact_search(kmer, pams_c, threshold, counting_callback, count);
      if (count > 1) return;
      gi_reverse.inexact_search(kmer, pams_c, threshold, counting_callback, count);
      if (count > 1) return;
    } else if (threshold > 0 && start) {
      gi_forward.inexact_search(kmer, pams, threshold, counting_callback, count);
      if (count > 1) return;
      gi_reverse.inexact_search(kmer, pams, threshold, counting_callback, count);
      if (count > 1) return;
    }

    std::vector<std::set<std::tuple<size_t, size_t, std::string>>> forward_off_targets_bwt(mismatches + 1);
    std::vector<std::set<std::tuple<size_t, size_t, std::string>>> reverse_off_targets_bwt(mismatches + 1);

    if (!start) {
      gi_forward.inexact_search(kmer, pams_c, mismatches, callback, forward_off_targets_bwt);
      gi_reverse.inexact_search(kmer, pams_c, mismatches, callback, reverse_off_targets_bwt);
    } else {
      gi_forward.inexact_search(kmer, pams, mismatches, callback, forward_off_targets_bwt);
      gi_reverse.inexact_search(kmer, pams, mismatches, callback, reverse_off_targets_bwt);
    }

    size_t genome_length = 0;
    for (size_t i = 0; i < gi_forward.gs.size(); i++) {
      genome_length += gi_forward.gs[i].length;
    }
        
    /*
     * This code resolves the position of the guide on the FORWARD
     * strand, making guides on the antisense strand negative so
     * that they can be distinguished.
     */

    std::vector<std::list<std::tuple<int64_t, std::string>>> off_targets(mismatches + 1);
    for (size_t i = 0; i < mismatches + 1; i++) {
      for (const auto& sp_ep : forward_off_targets_bwt[i]) {
        size_t sp = std::get<0>(sp_ep);
        size_t ep = std::get<1>(sp_ep);
        std::string match = std::get<2>(sp_ep);

        for (size_t j = sp; j <= ep; j++) {
          int64_t absolute_pos = -gi_forward.resolve(j);
          off_targets[i].push_back(std::make_tuple(absolute_pos, match));
        }
      }

      for (const auto& sp_ep : reverse_off_targets_bwt[i]) {
        size_t sp = std::get<0>(sp_ep);
        size_t ep = std::get<1>(sp_ep);
        std::string match = std::get<2>(sp_ep);

        for (size_t j = sp; j <= ep; j++) {
          int64_t absolute_pos = genome_length - (gi_reverse.resolve(j) + 1);
          off_targets[i].push_back(std::make_tuple(absolute_pos, match));
        }
      }
    }

    if (out_format == "csv") {
      std::string csv_lines = genomics::get_csv_lines(gi_forward, k, start, off_targets);
      output_mtx.lock();
      output << csv_lines;
      output_mtx.unlock();
    } else {
      std::string sam_line = genomics::get_sam_line(gi_forward, k, start, off_targets);
      output_mtx.lock();
      output << sam_line << std::endl;
      output_mtx.unlock();
    }
  }


  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  nlohmann::json search_kmer(const genome_index<t_wt, t_dens, t_inv_dens>& gi_forward,
                             const genome_index<t_wt, t_dens, t_inv_dens>& gi_reverse,
                             std::string kmer, size_t mismatches, int offtarget_lim) {
    using json = nlohmann::json;
    std::vector<std::set<std::tuple<size_t, size_t, std::string>>> forward_matches(mismatches + 1);
    std::vector<std::set<std::tuple<size_t, size_t, std::string>>> reverse_matches(mismatches + 1);

    gi_forward.inexact_search(kmer.begin(), kmer.end(), mismatches, callback, forward_matches);
    gi_reverse.inexact_search(kmer.begin(), kmer.end(), mismatches, callback, reverse_matches);

    size_t genome_length = 0;
    for (size_t i = 0; i < gi_forward.gs.size(); i++) {
      genome_length += gi_forward.gs[i].length;
    }

    json matches;
    int matches_enumerated = 0;
    for (size_t i = 0; i < mismatches + 1; i++) {
      for (const auto& sp_ep : forward_matches[i]) {
        size_t sp = std::get<0>(sp_ep);
        size_t ep = std::get<1>(sp_ep);
        for (size_t j = sp; j <= ep; j++) {
          if (offtarget_lim > 0 && matches_enumerated >= offtarget_lim) return matches;

          size_t absolute_pos = gi_forward.resolve(j);
          coordinates pos = resolve_absolute(gi_forward.gs, absolute_pos);
          json match = {
            {"chr", pos.chr.name},
            {"pos", pos.offset},
            {"absolute_pos", absolute_pos},
            {"strand", "+"},
            {"distance", i}
          };

          matches_enumerated++;
          matches.push_back(match);
        }
      }

      for (const auto& sp_ep : reverse_matches[i]) {
        size_t sp = std::get<0>(sp_ep);
        size_t ep = std::get<1>(sp_ep);
        for (size_t j = sp; j <= ep; j++) {
          if (offtarget_lim > 0 && matches_enumerated >= offtarget_lim) return matches;

          size_t absolute_pos = genome_length - (gi_reverse.resolve(j) + 1);
          coordinates pos = resolve_absolute(gi_forward.gs, absolute_pos);
          json match = {
            {"chr", pos.chr.name},
            {"absolute_pos", absolute_pos},
            {"pos", pos.offset},
            {"strand", "-"},
            {"distance", i}
          };

          matches_enumerated++;
          matches.push_back(match);
        }
      }
    }

    return matches;
  }

  /* Processes the kmers in the file, collecting all information
     about off targets and outputting it to a stream in SAM format. */
  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  void process_kmers_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi_forward,
                               const genome_index<t_wt, t_dens, t_inv_dens>& gi_reverse,
                               const std::vector<std::string> &pams,
                               size_t mismatches, int threshold, bool start,
                               std::string out_format,
                               const std::vector<kmer> &kmers,
                               std::ostream& output, std::mutex& output_mtx) {
    for (auto &kmer : kmers) {
      process_kmer_to_stream(gi_forward, gi_reverse, pams, mismatches, threshold,
                             start, out_format, kmer, output, output_mtx);
    }
  }
}

#endif /* PROCESS_H */
