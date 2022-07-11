#ifndef PROCESS_H
#define PROCESS_H

#include <set>
#include <list>
#include <tuple>

#include "json.hpp"
#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/printer.hpp"
#include "genomics/structures.hpp"
#include "guidescan.hpp"

namespace genomics {
  namespace {
    void off_target_enumerator(match m, std::vector<std::set<match>> &off_targets_bwt) {
      off_targets_bwt[m.mismatches].insert(m);
    }

    std::function<void(match, std::vector<std::set<match>>&)> callback = off_target_enumerator;

    void off_target_counter(match m, size_t &count) {
      count += m.ep - m.sp + 1;
    }

    std::function<void(match, size_t&)> counting_callback = off_target_counter;
  };

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi_forward,
                              const genome_index<t_wt, t_dens, t_inv_dens>& gi_reverse,
                              const enumerate_cmd_options& opts,
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
        
    std::vector<std::string> pams = opts.alt_pams;
    if (k.pam == "") {
      pams = {std::string("")};
    } else {
      pams.push_back(k.pam);
    }

    std::vector<std::string> pams_c;
    for (const auto& pam : pams) {
      pams_c.push_back(genomics::reverse_complement(pam));
    }
        
    std::string kmer = !opts.start ? genomics::reverse_complement(k.sequence) : k.sequence;

    // for this threshold search we need to not do bulges...
    if (opts.threshold > 0 && !opts.start) {
      gi_forward.inexact_search(kmer, pams_c, opts.threshold, 0, 0, 0, counting_callback, count);
      if (count > 1) return;
      gi_reverse.inexact_search(kmer, pams_c, opts.threshold, 0, 0, 0, counting_callback, count);
      if (count > 1) return;
    } else if (opts.threshold > 0 && opts.start) {
      gi_forward.inexact_search(kmer, pams, opts.threshold, 0, 0, 0, counting_callback, count);
      if (count > 1) return;
      gi_reverse.inexact_search(kmer, pams, opts.threshold, 0, 0, 0, counting_callback, count);
      if (count > 1) return;
    }

    std::vector<std::set<match>> forward_off_targets_bwt(opts.mismatches + 1);
    std::vector<std::set<match>> reverse_off_targets_bwt(opts.mismatches + 1);

    if (!opts.start) {
        gi_forward.inexact_search(kmer, pams_c, opts.mismatches, opts.rna_bulges, opts.dna_bulges, 1, callback, forward_off_targets_bwt);
        gi_reverse.inexact_search(kmer, pams_c, opts.mismatches, opts.rna_bulges, opts.dna_bulges, 1, callback, reverse_off_targets_bwt);
    } else {
        gi_forward.inexact_search(kmer, pams, opts.mismatches, opts.rna_bulges, opts.dna_bulges, 1, callback, forward_off_targets_bwt);
        gi_reverse.inexact_search(kmer, pams, opts.mismatches, opts.rna_bulges, opts.dna_bulges, 1, callback, reverse_off_targets_bwt);
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

    std::vector<std::list<std::tuple<int64_t, match>>> off_targets(opts.mismatches + 1);
    for (size_t i = 0; i < opts.mismatches + 1; i++) {
      for (const auto& m : forward_off_targets_bwt[i]) {
        for (size_t j = m.sp; j <= m.ep; j++) {
          int64_t absolute_pos = -gi_forward.resolve(j);
          off_targets[i].push_back(std::make_tuple(absolute_pos, m));
        }
      }

      for (const auto& m : reverse_off_targets_bwt[i]) {
        for (size_t j = m.sp; j <= m.ep; j++) {
          int64_t absolute_pos = genome_length - (gi_reverse.resolve(j) + 1);
          off_targets[i].push_back(std::make_tuple(absolute_pos, m));
        }
      }
    }

    if (opts.out_format == "csv") {
      std::string csv_lines = genomics::get_csv_lines(gi_forward, k, opts.start, off_targets);
      output_mtx.lock();
      output << csv_lines;
      output_mtx.unlock();
    } else {
      std::string sam_line = genomics::get_sam_line(gi_forward, k, opts.start, off_targets);
      output_mtx.lock();
      output << sam_line << std::endl;
      output_mtx.unlock();
    }
  }

  /* Processes the kmers in the file, collecting all information
     about off targets and outputting it to a stream in SAM format. */
  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  void process_kmers_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi_forward,
                               const genome_index<t_wt, t_dens, t_inv_dens>& gi_reverse,
                               const enumerate_cmd_options& opts,
                               const std::vector<kmer>& kmers,
                               std::ostream& output, std::mutex& output_mtx) {
    for (auto &kmer : kmers) {
        process_kmer_to_stream(gi_forward, gi_reverse, opts, kmer, output, output_mtx);
    }
  }
}

#endif /* PROCESS_H */
