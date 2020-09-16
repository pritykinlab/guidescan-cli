#include <cstdio>
#include <istream>
#include <sdsl/suffix_arrays.hpp>
#include "genome_index.hpp"
#include "seq_io.hpp"

void collecting_callback(size_t sp, size_t ep, std::vector<std::tuple<size_t, size_t>> &matches) {
    matches.push_back(std::make_tuple(sp, ep));
}

size_t hamming_distance(std::string s1, std::string s2) {
    size_t count = 0;
    for (auto i = 0; i < s1.length() && i < s2.length(); i++) {
        if (s1[i] != s2[i]) count++;
    }
    return count;
}

int main(int argc, char *argv[])
{
    using namespace std;

    string dna_seq("ATGCGAATCGATCGAGGAGCGATTATTCGGATCGGATCTAGCATGACGATCAG\n> asad 2\nActGGASnNtT");
    string fasta_header("> genome 1\n");
    istringstream fasta_is(fasta_header + dna_seq);

    ostringstream foo;

    auto gs = genomics::seq_io::parse_genome_structure(fasta_is);

    fasta_is.clear();
    fasta_is.seekg(0);
    genomics::seq_io::parse_sequence(fasta_is, foo);

    for (const auto &p : gs) {
        cout << get<0>(p) << " " << get<1>(p) << endl;
    }

    cout << "Sequence: " << foo.str().substr(0, 10) << "..." << endl;

    genomics::seq_io::write_to_file(gs, string("test.gs"));
    auto gs2 = genomics::seq_io::load_from_file(string("test.gs"));
    for (const auto &p : gs2) {
        cout << get<0>(p) << " " << get<1>(p) << endl;
    }

    return 0;
}
