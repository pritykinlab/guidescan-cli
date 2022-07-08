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

    char complement(char c) {
        switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return c;
        }
    }

    std::string complement(const std::string& kmer) {
        std::string s;
        for (int64_t i = 0; i <= kmer.length(); i++) {
            char c = kmer[i];
            s += complement(c);
        }

        return s;
    }

    std::string reverse_complement(const std::string& kmer) {
        std::string s;
        for (int64_t i = kmer.length() - 1; i >= 0; i--) {
            char c = kmer[i];
            s += complement(c);
        }

        return s;
    }
};
