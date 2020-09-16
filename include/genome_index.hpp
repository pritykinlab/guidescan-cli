#ifndef GENOME_INDEX_H
#define GENOME_INDEX_H

namespace genomics {
    typedef std::vector<std::tuple<std::string, size_t>> genome_structure;

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    class genome_index {
    public:
        typedef sdsl::csa_wt<t_wt, t_dens, t_inv_dens> t_csa;

        t_csa csa;
        genome_structure gs;

        const char* search_alphabet = "ATCG";
        size_t search_alphabet_size = strlen(search_alphabet);

        genome_index() {}
        genome_index(t_csa csa, genome_structure gs) : csa(csa), gs(gs)
        {}
        genome_index(const genome_index& other) : csa(csa), gs(gs)
        {}

        friend void swap(genome_index& first, genome_index& second)
        {
            using std::swap;
            swap(first.csa, second.csa);
            swap(first.gs, second.gs);
        }

        genome_index& operator=(genome_index other)
        {
            swap(*this, other);
            return *this;
        }

        size_t resolve(size_t bwt_position) {
            return csa[bwt_position];
        }

        template <class t_data>
        void inexact_search(typename std::string::iterator begin, typename std::string::iterator end,
                            size_t mismatches,
                            const std::function<void(size_t, size_t, t_data&)> callback,
                            t_data& data);

        template <class t_data>
        void inexact_search(typename std::string::iterator begin, typename std::string::iterator end,
                            size_t sp, size_t ep, size_t mismatches,
                            const std::function<void(size_t, size_t, t_data&)> callback,
                            t_data& data);
    };

};

#endif /* GENOME_INDEX_H */
