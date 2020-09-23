#ifndef PROCESS_H
#define PROCESS_H

#include <tuple>

#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/sam.hpp"

namespace genomics {
    namespace {
        void off_target_enumerator(size_t sp, size_t ep, size_t k,
				   std::vector<std::vector<std::tuple<size_t, size_t>>> &off_targets) {
            off_targets[k].push_back(std::make_tuple(sp, ep));
        }

        std::function<void(size_t, size_t, size_t, std::vector<std::vector<std::tuple<size_t, size_t>>>&)> callback = off_target_enumerator;
    };

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                const kmer& k,
                                std::ostream& output,
				std::mutex& output_mtx) {

        coordinates coords = resolve_absolute(gi.gs, k.absolute_coords);

        std::string kmer_pos = k.sequence + k.pam;
        std::string kmer_neg = reverse_complement(kmer_pos);

        std::vector<std::vector<std::tuple<size_t, size_t>>> off_targets_bwt(4);
        gi.inexact_search(kmer_pos, k.pam.length(), false, 3, callback, off_targets_bwt);
        if (off_targets_bwt[0].size() > 1 || off_targets_bwt[1].size() > 0) return;
        gi.inexact_search(kmer_neg, k.pam.length(), true, 3, callback, off_targets_bwt);
        if (off_targets_bwt[0].size() > 1 || off_targets_bwt[1].size() > 0) return;

        std::vector<std::vector<size_t>> off_targets(3);
	for (int i = 0; i < 3; i++) {
	    for (const auto& sp_ep : off_targets_bwt[i + 1]) {
		size_t sp = std::get<0>(sp_ep);
		size_t ep = std::get<1>(sp_ep);
		for (int j = sp; j <= ep; j++) {
		    off_targets[i].push_back(gi.resolve(j));
		}
	    }
	}
	
	output_mtx.lock();
	genomics::write_sam_line(output, gi, k, coords, off_targets);
	output_mtx.unlock();
    }

    /* Processes the kmers in the file, collecting all information
       about off targets and outputting it to a stream in SAM format. */
    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmers_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
				 std::unique_ptr<genomics::kmer_producer>& kmer_p, std::mutex& kmer_mtx,
                                 std::ostream& output, std::mutex& output_mtx) {
        kmer out_kmer;
        while (true) {
	    kmer_mtx.lock();
	    bool kmers_left = kmer_p->get_next_kmer(out_kmer);
	    kmer_mtx.unlock();

	    if (!kmers_left) break;
            process_kmer_to_stream(gi, out_kmer, output, output_mtx);
        }
    }
}

#endif /* PROCESS_H */
