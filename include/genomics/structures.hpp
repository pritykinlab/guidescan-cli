#ifndef GENOMIC_STRUCTURES_H
#define GENOMIC_STRUCTURES_H

#include <string>
#include <vector>

namespace genomics {
    enum class direction {positive, negative};

    struct kmer {
        std::string sequence;
        std::string pam;
        size_t absolute_coords;
        direction dir;
    };

    struct chromosome {
        std::string name;
        size_t length;

        bool operator==(const chromosome& other) const {
            return name == other.name && length == other.length;
        }
    };

    struct coordinates {
        chromosome chr;
        size_t offset;
    };

    typedef std::vector<chromosome> genome_structure;

    coordinates resolve_absolute(const genome_structure& gs, size_t absolute_coords);
    size_t      resolve_relative(const genome_structure& gs, coordinates coords);
};

#endif /* GENOMIC_STRUCTURES_H */
