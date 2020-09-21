#ifndef PROCESS_H
#define PROCESS_H

#include "kmer.hpp"

namespace genomics {
    namespace {
        void off_target_counter(size_t sp, size_t ep, size_t k,
                                std::vector<size_t> &off_targets) {
            off_targets[k] += ep - sp + 1;
        }

        std::function<void(size_t, size_t, size_t, std::vector<size_t>&)> callback = off_target_counter;
    };

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                const kmer& k,
                                std::ostream& output) {
        // if (sdsl::count(gi.csa, kmer.begin(), kmer.end()) > 1) return;

        // std::vector<size_t> off_targets(4);
        // gi.inexact_search(kmer.begin(), kmer.end(), 3, callback, off_targets);

        // output << kmer << " 1-" << off_targets[1] << ":2-" << off_targets[2]
        //        << ":3-" << off_targets[3] << std::endl; // race condition if multithreading
        output << k.sequence << k.pam << ":"
               << (k.dir == direction::positive ? "+" : "-") << std::endl;
    }

    /* Processes the kmers in the file, collecting all information
       about off targets and outputting it to a stream in BED format. */
    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmers_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                 const std::string& raw_sequence_file,
                                 size_t k,
                                 const std::string& pam,
                                 std::ostream& output,
                                 size_t start_position,
                                 size_t step_size) {
        std::ifstream sequence(raw_sequence_file);
        kmer_producer kmer_p(sequence, k, pam);


        kmer out_kmer;
        for (size_t idx = 0; kmer_p.get_next_kmer(out_kmer); idx++) {
            if ((idx - start_position) % step_size != 0) continue;
            process_kmer_to_stream(gi, out_kmer, output);
        }
    }
}

#endif /* PROCESS_H */
