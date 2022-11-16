#include <ostream>
#include <functional>

#include "genomics/kmer.hpp"
#include "genomics/sequences.hpp"
#include "genomics/seq_io.hpp"

namespace genomics {
    kmers_file_producer::kmers_file_producer(const std::string& kmers_file)
        : kmers_stream(kmers_file)
    {
        kmers_stream.read_header(io::ignore_no_column, "id", "sequence", "pam", "chromosome", "position", "sense");
    }
    
    bool kmers_file_producer::get_next_kmer(kmer& out_kmer) {
        std::string dir;
        bool ret = kmers_stream.read_row(out_kmer.id, out_kmer.sequence, out_kmer.pam,
                                         out_kmer.chromosome, out_kmer.position, dir);
        if (ret) {
          // Convert 1-indexed position to 0-indexed position
          out_kmer.position -= 1;
          out_kmer.dir = dir == "+" ? genomics::direction::positive : genomics::direction::negative;
        }
        return ret;
    }
};
