#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iterator>

#include "csv.hpp"
#include "genomics/sequences.hpp"
#include "genomics/seq_io.hpp"

namespace genomics {
    namespace seq_io {
        namespace {
            inline void ltrim(std::string &s) {
                s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
                    return !std::isspace(ch);
                }));
            }

            inline void rtrim(std::string &s) {
                s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
                    return !std::isspace(ch);
                }).base(), s.end());
            }

            inline char my_toupper(char ch)
            {
                return static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            }

            template <typename Out>
            void split(const std::string &s, char delim, Out result) {
                std::istringstream iss(s);
                std::string item;
                while (std::getline(iss, item, delim)) {
                    *result++ = item;
                }
            }

            std::vector<std::string> split(const std::string &s, char delim) {
                std::vector<std::string> elems;
                split(s, delim, std::back_inserter(elems));
                return elems;
            }

            void convert_raw_sequence(std::string& seq) {
                ltrim(seq);
                rtrim(seq);

                for (size_t i = 0; i < seq.length(); i++) {
                    seq[i] = toupper(seq[i]);
                }
            }
        };

        void parse_sequence(std::istream& fasta_is, std::ostream& sequence_os) {
            for(std::string line; std::getline(fasta_is, line);) {
                if (line.length() > 0 && line[0] == '>') continue;
                convert_raw_sequence(line);
                sequence_os << line;
            }
        }

        void reverse_complement_stream(std::istream& sequence_is, std::ostream& sequence_os) {
            std::string sequence;
            while (sequence_is) {
                sequence_is >> sequence;
            }

            sequence_os << genomics::reverse_complement(sequence);
        }

        genome_structure parse_genome_structure(std::istream& fasta_is) {
            genome_structure gs;

            std::string chromosome_name, line;
            std::getline(fasta_is, line);
            if (line.size() < 1 || line[0] != '>') {
                return gs;
            }

            chromosome_name = line.substr(1);
            ltrim(chromosome_name);
            rtrim(chromosome_name);
            auto words = split(chromosome_name, ' ');
            chromosome_name = words[0];
            size_t length = 0;
        
            while(std::getline(fasta_is, line)) {
                if (line[0] == '>') {
                    chromosome c = {chromosome_name, length};
                    gs.push_back(c);
                    chromosome_name = line.substr(1);
                    ltrim(chromosome_name);
                    rtrim(chromosome_name);
                    auto words = split(chromosome_name, ' ');
                    chromosome_name = words[0];
                    length = 0;
                    continue;
                }

                length += line.length();
            }

            chromosome c = {chromosome_name, length};
            gs.push_back(c);

            return gs;
        }

        void write_to_file(const genome_structure& gs, const std::string& filename){
            std::ofstream fs;
            fs.open(filename);
        
            for (auto p : gs) {
                fs << p.name;
                fs << "\n";
                fs << p.length;
                fs << "\n";
            }
        }

        bool load_from_file(genome_structure& gs, const std::string& filename){
            std::ifstream fs(filename);

            if (!fs) return false;
        
            while (fs) {
                std::string name, length_str;
                std::getline(fs, name);
                std::getline(fs, length_str);

                if (name.length() == 0 || length_str.length() == 0) {
                    break;
                }

                size_t length = std::stoll(length_str);
                chromosome c = {name, length};
                gs.push_back(c);
            }

            return true;
        }

        /* TODO: Make function work with CSV format and new kmer format */
        bool parse_kmer(io::CSVReader<6> kmers_stream, kmer& out_kmer) {
            // returns 0 on failure, 1 on success (interesting
            // decision...)
            std::string dir;
            bool ret = kmers_stream.read_row(out_kmer.id, out_kmer.sequence, out_kmer.pam,
                                             out_kmer.chromosome, out_kmer.position, dir);
            out_kmer.dir = dir == "+" ? genomics::direction::positive : genomics::direction::negative;
            return ret;
        }
    }; 
};
