#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <string> 

namespace genomics {
    bool pam_matches(const std::string& kmer, const std::string& pam);

    char complement(char c);
    std::string reverse_complement(const std::string& sequence);
}

#endif /* SEQUENCES_H */
