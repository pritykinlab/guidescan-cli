#include <ostream>
#include <functional>

#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/seq_io.hpp"

namespace genomics {
    seq_kmer_producer::seq_kmer_producer(const std::string& sequence_file, genome_structure gs,
                                         size_t k, const std::string &pam, size_t min_chr_length)
        : sequence(new std::ifstream(sequence_file)),
          gs(gs),
          min_chr_length(min_chr_length),
	  pam(pam),
	  k(k),
	  kmer_buffer(k + pam.length())
    {
    }

    size_t seq_kmer_producer::get_next_kmer(kmer& out_kmer) {
        while (kmer_queue.empty()) {
            sequence->seekg(stream_position, std::ios_base::beg);
            sequence->read(kmer_buffer.data(), kmer_buffer.size());

            if (!(*sequence)) return 0;

            std::string kmer_str(kmer_buffer.data(), kmer_buffer.size());
            std::string kmer_str_c = reverse_complement(kmer_str);

            if (pam_matches(kmer_str, pam)) {
                kmer kmer = {kmer_str.substr(0, k),
                             pam,
                             stream_position,
                             direction::positive};
                coordinates c = resolve_absolute(gs, kmer.absolute_coords);
                if (c.chr.length > min_chr_length) kmer_queue.push(kmer);
            }

            if (pam_matches(kmer_str_c, pam)) {
                kmer kmer = {kmer_str_c.substr(0, k),
                             pam,
                             stream_position,
                             direction::negative};
                coordinates c = resolve_absolute(gs, kmer.absolute_coords);
                if (c.chr.length > min_chr_length) kmer_queue.push(kmer);
            }

            stream_position++;
        }

        out_kmer = kmer_queue.front();
        kmer_queue.pop();
        return 1;
    }


    kmers_file_producer::kmers_file_producer(const std::string& kmers_file)
	: kmers_stream(new std::ifstream(kmers_file))
    {}
    
    size_t kmers_file_producer::get_next_kmer(kmer& out_kmer) {
	return seq_io::parse_kmer(*kmers_stream, out_kmer);
    }
};
