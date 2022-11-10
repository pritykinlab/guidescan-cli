/* 
   Defines functionality for finding kmers from a stream of nucleotides. 
*/

#ifndef KMER_H
#define KMER_H

#include <ostream>
#include <functional>

#include "csv.hpp"
#include "genomics/index.hpp"
#include "genomics/structures.hpp"

namespace genomics {
    class kmer_producer {
        /* 
           Gets the next available kmer, returning true on success and false
           if no kmers are left. Can throw exceptions.
        */
    public:
        virtual bool get_next_kmer(kmer& kmer) = 0;
    };

    /* 
       A kmer_producer constructed on top of a kmers file that contains a
       list of kmers (with 1-indexed positions) seperated by newlines.
    */
    class kmers_file_producer : public kmer_producer {
    private:
        io::CSVReader<6> kmers_stream;
    public:
        kmers_file_producer(const std::string& kmers_file);
        kmers_file_producer() = delete;

        /* 
           Gets the next available kmer, returning true on success and false
           if the stream is no longer valid for any reason
           (eofbit/badbit/failbit). 
        */
        bool get_next_kmer(kmer& out_kmer);
    };


};

#endif /* KMER_H */
