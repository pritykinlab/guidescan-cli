#include <ostream>
#include <functional>

#include "kmer.hpp"

namespace genomics {
    namespace  {
        bool pam_matches(const std::string& kmer,
                         const std::string& pam) {
            for (size_t i = 0; i < pam.length(); i++) {
                if (pam[i] == 'N') continue;
                if (kmer[kmer.length() - pam.length() + i] != pam[i]) return false;
            }

            return true;
        }

        std::string complement_strand(const std::string& kmer) {
            std::string s;
            for (int i = kmer.length() - 1; i >= 0; i--) {
                char c = kmer[i];

                switch (kmer[i]) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
                }

                s += c;
            }

            return s;
        }
    };

    kmer_producer::kmer_producer(std::istream& sequence, size_t k, const std::string &pam)
        : sequence(sequence), pam(pam), buffer_len(k + pam.length()),
          k(k), kmer_buffer(buffer_len)
    {}

    size_t kmer_producer::get_next_kmer(kmer& out_kmer) {
        if (!kmer_queue.empty()) {
            out_kmer = kmer_queue.front();
            kmer_queue.pop();
            return 1;
        }

        sequence.seekg(stream_position, std::ios_base::beg);
        sequence.read(kmer_buffer.data(), buffer_len);

        if (!sequence) return 0;

        std::string kmer_str(kmer_buffer.data());
        std::string kmer_str_c = complement_strand(kmer_str);

        if (pam_matches(kmer_str, pam)) {
            kmer kmer = {kmer_str.substr(0, k),
                         pam,
                         stream_position,
                         direction::positive};
            kmer_queue.push(kmer);
        }

        if (pam_matches(kmer_str_c, pam)) {
            kmer kmer = {kmer_str_c.substr(0, k),
                         pam,
                         stream_position,
                         direction::negative};
            kmer_queue.push(kmer);
        }

        stream_position++;
        return get_next_kmer(out_kmer); // Keep reading until EOF or error occurs.
    }
};
