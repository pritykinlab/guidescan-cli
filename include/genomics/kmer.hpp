/* 
   Defines functionality for finding kmers from a stream of nucleotides. 
*/

#ifndef KMER_H
#define KMER_H

#include <ostream>
#include <functional>

#include "genomics/index.hpp"
#include "genomics/structures.hpp"

namespace genomics {
    /* 
       A producer class for kmers built on top of std::istream. Returns
       the next kmer available if it can find one.

       It is assumed that the sequence represents a forward strand of
       DNA and the PAM sequence can match on either the forward or
       backward strand of it. 

       The PAM supports the wildcard 'N' matching any nucleotide.
     */
    class kmer_producer {
    private:
        std::istream& sequence;
        std::string pam;
        size_t k;

        size_t stream_position = 0;
        size_t buffer_len;
        std::vector<char> kmer_buffer;
        std::queue<kmer> kmer_queue;

    public:
        kmer_producer(std::istream& sequence, size_t k, const std::string &pam);
        kmer_producer() = delete;

        /* 
           Gets the next available kmer, returning 1 on success and 0
           if the stream is no longer valid for any reason
           (eofbit/badbit/failbit). 
        */
        size_t get_next_kmer(kmer& out_kmer);
    };
};

#endif /* KMER_H */
