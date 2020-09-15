#include <cstdio>
#include <istream>
#include <sdsl/suffix_arrays.hpp>

namespace genomics {
    typedef std::vector<std::tuple<std::string, size_t>> genome_structure;

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    class genome_index {
    public:
        sdsl::csa_wt<t_wt, t_dens, t_inv_dens> wt;
        genome_structure gs;

        const char* search_alphabet = "ATCG";
        size_t search_alphabet_size = strlen(search_alphabet);

        genome_index() {}
        genome_index(const genome_index& other) : wt(wt), gs(gs)
        {}

        friend void swap(genome_index& first, genome_index& second)
        {
            using std::swap;
            swap(first.wt, second.wt);
            swap(first.gs, second.gs);
        }

        genome_index& operator=(genome_index other)
        {
            swap(*this, other);
            return *this;
        }

        size_t resolve(size_t bwt_position) {
            return wt[bwt_position];
        }

        template <class t_data>
        void inexact_search(typename std::string::iterator begin, typename std::string::iterator end,
                            size_t mismatches, const std::function<void(size_t, size_t, t_data&)> callback,
                            t_data& data);

        template <class t_data>
        void inexact_search(typename std::string::iterator begin, typename std::string::iterator end,
                            size_t sp, size_t ep, size_t mismatches, const std::function<void(size_t, size_t, t_data&)> callback,
                            t_data& data);
    };

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    template <class t_data>
    void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::iterator begin,
                                                                typename std::string::iterator end,
                                                                size_t sp, size_t ep, size_t mismatches,
                                                                const std::function<void(size_t, size_t, t_data&)> callback,
                                                                t_data& data) {
        if (begin == end) {
            callback(sp, ep, data);
        }

        char c = *(end - 1);

        size_t occ_before = wt.rank_bwt(sp, c);
        size_t occ_within = wt.rank_bwt(ep + 1, c) - occ_before;

        if (occ_within > 0) {
            size_t sp_prime = wt.C[wt.char2comp[c]] + occ_before;
            size_t ep_prime = sp_prime + occ_within - 1;
            inexact_search(begin, end - 1, sp_prime, ep_prime, mismatches,
                           callback, data);
        }

        if (mismatches < 1) return;

        for (size_t i = 0; i < search_alphabet_size; i++) {
            if (search_alphabet[i] == c) continue;

            char a = search_alphabet[i];

            occ_before = wt.rank_bwt(sp, a);
            occ_within = wt.rank_bwt(ep + 1, a) - occ_before;

            if (occ_within > 0) {
                size_t sp_prime = wt.C[wt.char2comp[a]] + occ_before;
                size_t ep_prime = sp_prime + occ_within - 1;
                inexact_search(begin, end - 1, sp_prime, ep_prime, mismatches - 1,
                               callback, data);
            }
        }
    }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    template <class t_data>
    void genome_index<t_wt, t_dens, t_inv_dens>::inexact_search(typename std::string::iterator begin,
                                                                typename std::string::iterator end,
                                                                size_t mismatches,
                                                                const std::function<void(size_t, size_t, t_data&)> callback,
                                                                t_data& data) {

        inexact_search(begin, end, 0, wt.size() - 1, mismatches, callback, data);
    }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    bool load_genome_index(genome_index<t_wt, t_dens, t_inv_dens>& index, std::istream& fasta_is) {
        char* tmp_name = std::tmpnam(nullptr);

        FILE* tmpf = std::fopen(tmp_name, "a");
        if (tmpf == nullptr) return true;

        genome_structure gs;

        std::string header;
        std::getline(fasta_is, header);
        if (header.size() < 1 || header[0] != '>') {
            std::fclose(tmpf);
            return true;
        }


        header = header.substr(1);
        size_t offset = 0;
        
        for(std::string line; std::getline(fasta_is, line);) {
            if (line[0] == '>') {
                gs.push_back(std::make_tuple(header, offset));
                header = line.substr(1);
                offset = 0;
                continue;
            }

            std::fwrite(line.c_str(), sizeof line[0], line.size(), tmpf);
            offset += line.size();
        }

        std::fclose(tmpf);

        std::cout << tmp_name << std::endl;
        sdsl::construct(index.wt, std::string(tmp_name), 1);
        index.gs = gs;

        return false;
    }
};

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

    string dna_seq("ATGCGAATCGATCGAGGAGCGATTATTCGGATCGGATCTAGCATGACGATCAG");
    string fasta_header("> genome 1\n");
    istringstream fasta_is(fasta_header + dna_seq);

    genomics::genome_index<sdsl::wt_huff<>, 32, 4096> genome_index;
    bool error = genomics::load_genome_index(genome_index, fasta_is);

    if (error) {
        cout << "ERROR!" << endl;
        return 1;
    }

    string query("TCGAGG");
    vector<tuple<size_t, size_t>> matches;
    function<void(size_t, size_t, vector<tuple<size_t, size_t>>&)> callback = collecting_callback;

    genome_index.inexact_search(query.begin(), query.end(), 6, callback, matches);
    cout << dna_seq.length() - 6 << endl;

    for (auto match : matches) {
        for (auto i = get<0>(match); i <= get<1>(match); i++) {
            string match_str = dna_seq.substr(genome_index.resolve(i), 6);
            size_t dist = hamming_distance(match_str, query);
            cout << match_str << " " << dist << endl;
        }
    }

    return 0;
}
