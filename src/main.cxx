#include <filesystem>
#include <istream>

#include <sdsl/suffix_arrays.hpp>

#include "genome_index.hpp"
#include "seq_io.hpp"
#include "kmer.hpp"

void collecting_callback(size_t sp, size_t ep, size_t k, std::vector<std::tuple<size_t, size_t>> &matches) {
    matches.push_back(std::make_tuple(sp, ep));
}

size_t hamming_distance(std::string s1, std::string s2) {
    size_t count = 0;
    for (auto i = 0; i < s1.length() && i < s2.length(); i++) {
        if (s1[i] != s2[i]) count++;
    }
    return count;
}

bool file_exists(const std::string& fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int main(int argc, char *argv[])
{
    using namespace std;

    if (argc < 1) {
        cout << "Usage: ./main [fasta_file]" << endl;
        return 1;
    }

    string fasta_file(argv[1]);
    string raw_sequence_file = fasta_file + ".dna";
    string genome_structure_file = fasta_file + ".gs";
    string fm_index_file = fasta_file + ".csa";

    ifstream fasta_is(fasta_file);
    if (!fasta_is) {
        cerr << "ERROR: FASTA file \"" << fasta_file
             << "\" does not exist." << endl;
        return 1;
    }

    if (!file_exists(raw_sequence_file)) {
        ofstream os(raw_sequence_file);
        if (!os) {
            cerr << "ERROR: Could not create raw sequence file." << endl;
            return 1;
        }

        cout << "No raw sequence file \"" << raw_sequence_file
             << "\". Building now..." << endl;
        genomics::seq_io::parse_sequence(fasta_is, os);
    }

    genomics::genome_structure gs;
    if (!genomics::seq_io::load_from_file(gs, genome_structure_file)) {
        cout << "No genome structure file \"" << genome_structure_file
             << "\" located. Building now..." << endl;

        fasta_is.clear();
        fasta_is.seekg(0);
        gs = genomics::seq_io::parse_genome_structure(fasta_is);
        genomics::seq_io::write_to_file(gs, genome_structure_file);
    }

    sdsl::csa_wt<sdsl::wt_huff<>, 32, 8192> fm_index;
    if (!load_from_file(fm_index, fm_index_file)) {
        cout << "No index file \"" << fm_index_file
             << "\" located. Building now..." << endl;

        construct(fm_index, raw_sequence_file, 1);
        store_to_file(fm_index, fm_index_file);
    }   


    genomics::genome_index<sdsl::wt_huff<>, 32, 8192> gi(fm_index, gs);
    cout << "Successfully loaded index." << endl;
                     
    genomics::process_kmers_to_stream(gi, raw_sequence_file, 20,
                                      string("NGG"), cout, 0, 1);

    return 0;
}
