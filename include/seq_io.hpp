#ifndef SEQ_IO_H
#define SEQ_IO_H

#include <vector>
#include <tuple>
#include <istream> 
#include <ostream> 

namespace genomics::seq_io {
    typedef std::vector<std::tuple<std::string, size_t>> genome_structure;

    /* Parses a FASTA input stream to a raw sequence of uppercase
       genomic symbols and sends it to an output stream. */
    void parse_sequence(std::istream& fasta_is, std::ostream& sequence_os);

    /* Parses a FASTA input stream to a genome_structure object that
       describes how absolute coordinates map to relative coordinates
       within the genome */
    genome_structure parse_genome_structure(std::istream& fasta_is);

};

#endif /* SEQ_IO_H */
