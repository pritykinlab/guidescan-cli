#ifndef SEQ_IO_H
#define SEQ_IO_H

#include <vector>
#include <tuple>
#include <istream> 
#include <ostream> 

#include "genomics/index.hpp"

// TODO: Handle error cases
namespace genomics {
    namespace seq_io {
        /* Parses a FASTA input stream to a raw sequence of uppercase
           genomic symbols and sends it to an output stream. */
        void parse_sequence(std::istream& fasta_is, std::ostream& sequence_os);

        /* Parses a FASTA input stream to a genome_structure object that
           describes how absolute coordinates map to relative coordinates
           within the genome */
        genome_structure parse_genome_structure(std::istream& fasta_is);

	/* Parses a kmer from an input stream, returning 1 on success,
	   0 otherwise. */
	size_t parse_kmer(std::istream& kmers_stream, kmer& out_kmer);

        void write_to_file(const genome_structure& gs, const std::string& filename);
        bool load_from_file(genome_structure& gs, const std::string& filename);

	void write_to_file(const std::vector<kmer>& kmers, const std::string& filename);
	bool load_from_file(std::vector<kmer>& kmers, const std::string& filename);
    };
};

#endif /* SEQ_IO_H */
