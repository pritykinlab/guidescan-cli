#ifndef GENOMIC_STRUCTURES_H
#define GENOMIC_STRUCTURES_H

#include <string>
#include <vector>

namespace genomics {
    enum class direction {positive, negative};

    struct kmer {
        std::string id;
        std::string sequence;
        std::string pam;
        std::string chromosome;
        uint64_t position;
        direction dir;
    };

    struct chromosome {
        std::string name;
        uint64_t length;

        bool operator==(const chromosome& other) const {
            return name == other.name && length == other.length;
        }
    };

    struct coordinates {
        chromosome chr;
        uint64_t offset;
    };

    struct match {
        std::string sequence;
        uint64_t sp, ep;
        uint64_t mismatches;
        uint64_t dna_bulges;
        uint64_t rna_bulges;

        bool operator<(const match& other) const {
            return sequence < other.sequence;
        } 
    };

    typedef std::vector<chromosome> genome_structure;

    coordinates resolve_absolute(const genome_structure& gs, uint64_t absolute_coords);
    uint64_t    resolve_relative(const genome_structure& gs, coordinates coords);
};

#endif /* GENOMIC_STRUCTURES_H */
