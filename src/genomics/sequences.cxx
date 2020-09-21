#include "genomics/sequences.hpp"

namespace genomics {
    bool pam_matches(const std::string& kmer,
                     const std::string& pam) {
        for (size_t i = 0; i < pam.length(); i++) {
            if (pam[i] == 'N') continue;
            if (kmer[kmer.length() - pam.length() + i] != pam[i]) return false;
        }

        return true;
    }

    std::string reverse_complement(const std::string& kmer) {
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
