#ifndef KMER_H
#define KMER_H

#include <ostream>
#include <functional>

#include "genome_index.hpp"

namespace genomics {
    class kmer {
    public:
        std::string sequence;
        size_t absolute_position;
    };

    namespace {
        bool pam_matches_reverse(const std::string& kmer,
                                 const std::string& pam) {
            for (size_t i = 0; i < pam.length(); i++) {
                if (pam[i] == 'N') continue;
                if (kmer[i] != pam[pam.length() - 1 - i]) return false;
            }

            return true;
        }

        bool pam_matches_foward(const std::string& kmer,
                                const std::string& pam) {
            for (size_t i = 0; i < pam.length(); i++) {
                if (pam[i] == 'N') continue;
                if (kmer[kmer.length() - pam.length() + i] != pam[i]) return false;
            }

            return true;
        }

        
        void off_target_counter(size_t sp, size_t ep, size_t k,
                                std::vector<size_t> &off_targets) {
            off_targets[k] += ep - sp + 1;
        }

        std::function<void(size_t, size_t, size_t, std::vector<size_t>&)> callback = off_target_counter;
    };

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void process_kmer_to_stream(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                                const std::string& kmer,
                                const std::string& pam,
                                std::ostream& output) {
        if (sdsl::count(gi.csa, kmer.begin(), kmer.end()) > 1) return;

        std::vector<size_t> off_targets(4);
        gi.inexact_search(kmer.begin(), kmer.end(), 3, callback, off_targets);

        output << kmer << " 1-" << off_targets[1] << ":2-" << off_targets[2]
               << ":3-" << off_targets[3] << std::endl;
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

        std::ifstream is(raw_sequence_file);
        if (!is) {
            throw std::runtime_error("Sequence file not available.");
        }

        std::vector<char> kmer_buffer(k);
        size_t position = 1;

        for (; is; position++) {
            is.read(kmer_buffer.data(), k);
            is.seekg(1 - k, std::ios_base::cur);
            
            std::string kmer(kmer_buffer.data());

            if (position < start_position || position % step_size != 0)
                continue;

            if (pam_matches_foward(kmer, pam) ||
                pam_matches_reverse(kmer, pam)) {
                process_kmer_to_stream(gi, kmer, pam, output);
            }

            //if (position % 1000 == 0) std::cerr << "PROGRESS: " << position << std::endl;
            if (position > 500000) break;
        }
    }
}

#endif /* KMER_H */
