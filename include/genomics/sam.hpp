#ifndef SAM_H
#define SAM_H

#include "genomics/structures.hpp"
#include "genomics/sequences.hpp"

#include <iostream>
#include <vector>

namespace genomics {
    namespace {
        std::string byte_to_big_endian_hex(unsigned char n) {
            std::string out("");
            int lowest1 = n % 16;
            n           = n / 16;
            int lowest2 = n % 16;

            char c1 = 'f';
            switch(lowest1) {
            case 0:  c1 = '0'; break;
            case 1:  c1 = '1'; break;
            case 2:  c1 = '2'; break;
            case 3:  c1 = '3'; break;
            case 4:  c1 = '4'; break;
            case 5:  c1 = '5'; break;
            case 6:  c1 = '6'; break;
            case 7:  c1 = '7'; break;
            case 8:  c1 = '8'; break;
            case 9:  c1 = '9'; break;
            case 10: c1 = 'a'; break;
            case 11: c1 = 'b'; break;
            case 12: c1 = 'c'; break;
            case 13: c1 = 'd'; break;
            case 14: c1 = 'e'; break;
            case 15: c1 = 'f'; break;
            };

            char c2 = 'f';
            switch(lowest2) {
            case 0:  c2 = '0'; break;
            case 1:  c2 = '1'; break;
            case 2:  c2 = '2'; break;
            case 3:  c2 = '3'; break;
            case 4:  c2 = '4'; break;
            case 5:  c2 = '5'; break;
            case 6:  c2 = '6'; break;
            case 7:  c2 = '7'; break;
            case 8:  c2 = '8'; break;
            case 9:  c2 = '9'; break;
            case 10: c2 = 'a'; break;
            case 11: c2 = 'b'; break;
            case 12: c2 = 'c'; break;
            case 13: c2 = 'd'; break;
            case 14: c2 = 'e'; break;
            case 15: c2 = 'f'; break;
            };

            out += c2;
            out += c1;

            return out;
        }

        std::string num_to_little_endian_hex(uint64_t num) {
            std::string out("");

            for (size_t i = 0; i < sizeof(num); i++) {
                unsigned char lowest = num % 256;
                num                  = num / 256;
                out += byte_to_big_endian_hex(lowest);
            }

            return out;
        }

        std::string vector_to_little_endian_hex(std::vector<uint64_t> v) {
            std::string out("");
            for (auto n : v) {
                out += num_to_little_endian_hex(n);
            }
            return out;
        }

        int64_t get_delim(genome_structure gs) {
            int64_t delim = 0;
            for (const auto& chr : gs) {
                delim += chr.length;
            } 
            return -(delim + 1);
        }

        std::string off_target_string(genome_structure gs,
                                      const std::vector<std::vector<int64_t>>& off_targets) {
            uint64_t delim = get_delim(gs);
            std::vector<uint64_t> v;

            for (uint64_t k = 0; k < off_targets.size(); k++) {
                v.insert(v.end(), off_targets[k].begin(), off_targets[k].end());
                v.push_back(k);
                v.push_back(delim);
            }

            return vector_to_little_endian_hex(v);
        }
    };

    void write_sam_header(std::ostream& os, genome_structure gs) {
        os << "@HD\tVN:1.0\tSO:unknown" << std::endl;
        for (const auto& chr : gs) {
            os << "@SQ\tSN:" << chr.name << "\tLN:" << chr.length << std::endl;
        }
    }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    std::string get_sam_line(genome_index<t_wt, t_dens, t_inv_dens> gi,
                             const kmer& k,
                             const std::vector<std::vector<int64_t>>& off_targets) {
        std::string sequence(k.sequence + k.pam);
        std::string samline(k.id);

        samline += "\t";
        samline += k.dir == direction::positive ? "0" : "16";
        samline += "\t" + k.chromosome;

        samline += "\t" + std::to_string(k.position);
        samline += "\t100";
        samline += "\t" + std::to_string(sequence.length()) + "M";
        samline += "\t*\t0\t0";

        if (k.dir == direction::negative) {
            samline += "\t" + reverse_complement(sequence);
        } else {
            samline += "\t" + sequence;
        }

        samline += "\t*";

        bool no_off_targets = true;
        for (const auto& v : off_targets) {
            if (!v.empty()) no_off_targets = false;
        }

        if (!no_off_targets) {
            std::string ots = off_target_string(gi.gs, off_targets);
            samline += "\tof:H:" + ots;
        }

        return samline;
    }
};

#endif /* SAM_H */
