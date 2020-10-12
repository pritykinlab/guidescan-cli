#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <random>

#define t_sa_dens 32
#define t_isa_dens 8192

typedef sdsl::enc_vector<sdsl::coder::elias_delta, 16> t_wt;

int main(int argc, char* argv[]) {
    if (argc <  2) {
        std::cout << "Usage: " << argv[0] << " text_file " << std::endl;
        return 1;
    }

    std::string index_suffix = ".csa";
    std::string index_file   = std::string(argv[1]) + index_suffix;
    sdsl::csa_sada<t_wt, t_sa_dens, t_isa_dens> fm_index;

    if (!load_from_file(fm_index, index_file)) {
        std::ifstream in(argv[1]);
        if (!in) {
            std::cout << "ERROR: File " << argv[1] << " does not exist. Exit." << std::endl;
            return 1;
        }
        std::cout << "No index "<<index_file<< " located. Building index now." << std::endl;
        sdsl::construct(fm_index, argv[1], 1);
        sdsl::store_to_file(fm_index, index_file);
    }

    std::random_device r;
    std::default_random_engine e(r());
    std::uniform_int_distribution<int> uniform_dist(0, fm_index.size());

    for (int i = 0; i < 1000000; i++) {
        size_t position = uniform_dist(e);

        size_t occ1 = fm_index.rank_bwt('A', position);
        size_t occ2 = fm_index.rank_bwt('A', position + 1);

        if (occ2 < occ1) {
            std::cout << "Oh no..." << std::endl;
        }
    }

    return 0;
}
