#include <ostream>
#include <functional>

#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/seq_io.hpp"

namespace genomics {
    seq_kmer_producer::seq_kmer_producer(const std::string& sequence_file, size_t k, const std::string &pam)
        : sequence(new std::ifstream(sequence_file)),
	  pam(pam),
	  buffer_len(k + pam.length()),
	  k(k),
	  kmer_buffer(buffer_len)
    {
    }

    size_t seq_kmer_producer::get_next_kmer(kmer& out_kmer) {
        if (!kmer_queue.empty()) {
            out_kmer = kmer_queue.front();
            kmer_queue.pop();
            return 1;
        }

        sequence->seekg(stream_position, std::ios_base::beg);
        sequence->read(kmer_buffer.data(), buffer_len);

        if (!(*sequence)) return 0;

        std::string kmer_str(kmer_buffer.data());
        std::string kmer_str_c = reverse_complement(kmer_str);

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


    kmers_file_producer::kmers_file_producer(const std::string& kmers_file)
	: kmers_stream(new std::ifstream(kmers_file))
    {}
    
    size_t kmers_file_producer::get_next_kmer(kmer& out_kmer) {
	return seq_io::parse_kmer(*kmers_stream, out_kmer);
    }
};
