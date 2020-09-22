#ifndef PROCESS_H
#define PROCESS_H

#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"

namespace genomics {
    namespace {
        void off_target_enumerator(size_t sp, size_t ep, size_t k,
                                   std::vector<std::tuple<size_t, size_t>> &off_targets) {
            if (k == 1) off_targets.push_back(std::make_tuple(sp, ep));
        }

        void off_target_counter(size_t sp, size_t ep, size_t k,
                                std::vector<size_t> &off_targets) {
            off_targets[k] += ep - sp + 1;
        }

        std::function<void(size_t, size_t, size_t, std::vector<std::tuple<size_t, size_t>>&)> callback_enum = off_target_enumerator;
        std::function<void(size_t, size_t, size_t, std::vector<size_t>&)> callback = off_target_counter;
    };

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                const kmer& k,
                                std::ostream& output,
				std::mutex& output_mtx) {

        coordinates coords = resolve_absolute(gi.gs, k.absolute_coords);

        if (coords.chr.length < 10000) {
            return;
        }

        std::string kmer_pos = k.sequence + k.pam;
        std::string kmer_neg = reverse_complement(kmer_pos);

        std::vector<size_t> off_targets(4, 0);
        gi.inexact_search(kmer_pos, k.pam.length(), false, 3, callback, off_targets);
        gi.inexact_search(kmer_neg, k.pam.length(), true, 3, callback, off_targets);

        if (off_targets[0] > 1) return;
        if (off_targets[1] > 0) return;

	output_mtx.lock();
        output << k.sequence << k.pam << ":"
               << (k.dir == direction::positive ? "+" : "-")
               << " 1-" << off_targets[1] << ":2-" << off_targets[2]
               << ":3-" << off_targets[3] << std::endl; // race condition if multithreading
	output_mtx.unlock();

    }

    /* Processes the kmers in the file, collecting all information
       about off targets and outputting it to a stream in BED format. */
    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmers_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                 const std::string& raw_sequence_file,
                                 size_t k,
                                 const std::string& pam,
                                 std::ostream& output, std::mutex& output_mtx,
                                 size_t start_position,
                                 size_t step_size) {
        std::ifstream sequence(raw_sequence_file);
        kmer_producer kmer_p(sequence, k, pam);
        kmer out_kmer;
        for (size_t idx = 0; kmer_p.get_next_kmer(out_kmer); idx++) {
            if ((idx - start_position) % step_size != 0) continue;
            process_kmer_to_stream(gi, out_kmer, output, output_mtx);
        }
    }
}

#endif /* PROCESS_H */
