#ifndef SAM_H
#define SAM_H

#include "genomics/structures.hpp"
#include "genomics/sequences.hpp"

#include <iostream>
#include <vector>

namespace genomics {
    void write_sam_header(std::ostream& os, genome_structure gs) {
	os << "@HD\tVN:1.0\tSO:unknown" << std::endl;
	for (const auto& chr : gs) {
	    os << "@SQ\tSN:" << chr.name << "\tLN:" << chr.length << std::endl;
	}
    }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    void write_sam_line(std::ostream& os, genome_index<t_wt, t_dens, t_inv_dens> gi,
			const kmer& k, const coordinates& coords,
			const std::vector<size_t>& off_targets) {
	std::string sequence(k.sequence + k.pam);
	std::string samline(sequence);

	samline += "\t";
	samline += k.dir == direction::positive ? "0" : "16";
	samline += "\t" + coords.chr.name;

	size_t position = coords.offset + 1;
	if (k.dir == direction::negative) {
	    //position -= sequence.length() - 1;
	}

	samline += "\t" + std::to_string(position);
	samline += "\t100";
	samline += "\t" + std::to_string(sequence.length()) + "M";
	samline += "\t*\t0\t0";

	if (k.dir == direction::negative) {
	    samline += "\t" + reverse_complement(sequence);
	} else {
	    samline += "\t" + sequence;
	}

	samline += "\t*\t";

	os << samline << std::endl;
    }
};

#endif /* SAM_H */
