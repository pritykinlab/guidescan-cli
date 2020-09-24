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
    
    class kmer_producer {
        /* 
           Gets the next available kmer, returning 1 on success and 0
	   if  no kmers are left. Can throw exceptions.
        */
    public:
	virtual size_t get_next_kmer(kmer& kmer) = 0;
    };

    /* 
       A producer class for kmers built on top of
       std::istream. Returns the next kmer available if it can find
       one.

       It is assumed that the sequence represents a forward strand of
       DNA and the PAM sequence can match on either the forward or
       backward strand of it.

       The PAM supports the wildcard 'N' matching any nucleotide.
    */
    class seq_kmer_producer : public kmer_producer {
    private:
	std::unique_ptr<std::istream> sequence;
        std::string pam;
        genome_structure gs;
        size_t k, min_chr_length;

        size_t stream_position = 0;
        std::vector<char> kmer_buffer;
        std::queue<kmer> kmer_queue;

    public:
        seq_kmer_producer(const std::string& sequence_file, genome_structure gs,
                          size_t k, const std::string &pam, size_t min_chr_length);
        seq_kmer_producer() = delete;

        /* 
           Gets the next available kmer, returning 1 on success and 0
           if the stream is no longer valid for any reason
           (eofbit/badbit/failbit). 
        */
        size_t get_next_kmer(kmer& out_kmer);
    };

    /* 
       A kmer_producer constructed on top of a kmers file that contains a
       list of kmers seperated by newlines.
    */
    class kmers_file_producer : public kmer_producer {
    private:
	std::unique_ptr<std::istream> kmers_stream;
    public:
        kmers_file_producer(const std::string& kmers_file);
        kmers_file_producer() = delete;

        /* 
           Gets the next available kmer, returning 1 on success and 0
           if the stream is no longer valid for any reason
           (eofbit/badbit/failbit). 
        */
        size_t get_next_kmer(kmer& out_kmer);
    };


};

#endif /* KMER_H */
