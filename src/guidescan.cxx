#include <thread>
#include <istream>
#include <memory>

#include <sdsl/suffix_arrays.hpp>

#include "json.hpp"
#include "httplib.h"
#include "CLI/CLI.hpp"
#include "genomics/index.hpp"
#include "genomics/sam.hpp"
#include "genomics/seq_io.hpp"
#include "genomics/process.hpp"
#include "genomics/kmer.hpp"

#define t_sa_dens 64
#define t_isa_dens 8192

typedef sdsl::wt_huff<> t_wt;

struct index_cmd_options {
    std::string fasta_file;
    CLI::Option* fasta_file_opt = nullptr;

    std::string index_file_prefix;
    CLI::Option* index_file_prefix_opt = nullptr;
};

struct enumerate_cmd_options {
    int threshold;
    CLI::Option* threshold_opt = nullptr;

    size_t mismatches;
    CLI::Option* mismatches_opt = nullptr;
    
    size_t nthreads;
    CLI::Option* nthreads_opt = nullptr;

    bool enumerate_all;
    CLI::Option* enumerate_all_opt = nullptr;

    std::string index_file_prefix;
    CLI::Option* index_file_prefix_opt = nullptr;

    std::string database_file;
    CLI::Option* database_file_opt = nullptr;

    std::string kmers_file;
    CLI::Option* kmers_file_opt = nullptr;

    std::vector<std::string> alt_pams;
    CLI::Option* alt_pams_opt = nullptr;
};

struct http_server_cmd_options {
    std::string index_file_prefix;
    CLI::Option* index_file_prefix_opt = nullptr;

    size_t mismatches;
    CLI::Option* mismatches_opt = nullptr;

    int ot_limit;
    CLI::Option* ot_limit_opt = nullptr;

    size_t port;
    CLI::Option* port_opt = nullptr;
};

CLI::App* index_cmd(CLI::App &guidescan, index_cmd_options& opts) {
    auto build  = guidescan.add_subcommand("index", "Builds an genomic index over a FASTA file.");

    opts.index_file_prefix_opt = build->add_option("--index", opts.index_file_prefix, "Genomic index prefix.", true);
    opts.fasta_file_opt  = build->add_option("genome", opts.fasta_file, "Genome in FASTA format")
        ->check(CLI::ExistingFile)
        ->required();

    return build;
}

CLI::App* enumerate_cmd(CLI::App &guidescan, enumerate_cmd_options& opts) {
    auto build = guidescan.add_subcommand("enumerate", "Enumerates off-targets against a reference.");

    opts.nthreads    = std::thread::hardware_concurrency();
    opts.alt_pams    = {};
    opts.threshold   = 1;
    opts.mismatches  = 3;

    opts.nthreads_opt    = build->add_option("-n,--threads", opts.nthreads, "Number of threads to parallelize over", true);
    opts.alt_pams_opt    = build->add_option("-a,--alt-pam", opts.alt_pams, "Alternative PAMs used to find off-targets", true);
    opts.mismatches_opt  = build->add_option("-m,--mismatches", opts.mismatches, "Number of mismatches to allow when finding off-targets", true);
    opts.threshold_opt   = build->add_option("-t,--threshold", opts.threshold, "Filters gRNAs with off-targets at a distance at or below this threshold", true);
    opts.kmers_file_opt  = build->add_option("-f,--kmers-file", opts.kmers_file,
					     "File containing kmers to build gRNA database"
					     " over, if not specified, will generate the database over all kmers with the given PAM")
        ->check(CLI::ExistingFile)
        ->required();
    opts.index_file_prefix_opt  = build->add_option("index", opts.index_file_prefix, "Prefix for index files")
        ->required();
    opts.database_file_opt = build->add_option("-o, --output", opts.database_file, "Output database file.")
        ->required();
  
    return build;
}

CLI::App* http_cmd(CLI::App &guidescan, http_server_cmd_options& opts) {
    auto http = guidescan.add_subcommand("http-server",
                                         "Starts a local HTTP server to receive gRNA processing requests.");
    opts.mismatches  = 3;
    opts.port = 4500;
    opts.ot_limit = -1;

    opts.port_opt       = http->add_option("--port", opts.port, "HTTP Server Port", true);
    opts.mismatches_opt = http->add_option("-m,--mismatches", opts.mismatches, "Number of mismatches to allow when finding off-targets", true);
    opts.ot_limit_opt   = http->add_option("-l,--offtarget-limit", opts.ot_limit, "Max number of off-targets to enumerate", true);
    opts.index_file_prefix_opt = http->add_option("genome", opts.index_file_prefix, "Genome in FASTA format")
        ->check(CLI::ExistingFile)
        ->required();

    return http;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

bool file_exists(const std::string& fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int do_index_cmd(const index_cmd_options& opts) {
    using namespace std;

    string index_file_prefix = opts.fasta_file + ".index";
    if (!(*opts.index_file_prefix_opt)) {
        cout << "Index file prefix not specified. Using default: " << index_file_prefix << endl;
    } else {
        index_file_prefix = opts.index_file_prefix;
    }

    string forward_raw_sequence_file = opts.fasta_file + ".forward.dna";
    string reverse_raw_sequence_file = opts.fasta_file + ".reverse.dna";

    string forward_fm_index_file = index_file_prefix + ".forward";
    string reverse_fm_index_file = index_file_prefix + ".reverse";
    string genome_structure_file = index_file_prefix + ".gs";

    ifstream fasta_is(opts.fasta_file);
    if (!fasta_is) {
        cerr << "ERROR: FASTA file \"" << opts.fasta_file
             << "\" does not exist." << endl;
        return 1;
    }

    cout << "Attempting to read raw sequence file (if constructed)..." << endl;
    if (!file_exists(forward_raw_sequence_file)) {
        ofstream os(forward_raw_sequence_file);
        if (!os) {
            cerr << "ERROR: Could not create forward raw sequence file." << endl;
            return 1;
        }

        cout << "No raw sequence file \"" << forward_raw_sequence_file
             << "\". Building now..." << endl;
        genomics::seq_io::parse_sequence(fasta_is, os);
    }

    if (!file_exists(reverse_raw_sequence_file)) {
        ofstream os(reverse_raw_sequence_file);
        if (!os) {
            cerr << "ERROR: Could not create reverse raw sequence file." << endl;
            return 1;
        }

        cout << "No raw sequence file \"" << reverse_raw_sequence_file
             << "\". Building now..." << endl;
        ifstream is(forward_raw_sequence_file);
        genomics::seq_io::reverse_complement_stream(is, os);
    }

    cout << "Constructing genome index..." << endl;
    genomics::genome_structure gs;
    gs = genomics::seq_io::parse_genome_structure(fasta_is);
    genomics::seq_io::write_to_file(gs, genome_structure_file);

    cout << "Constructing forward genomic index." << endl;
    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> forward_fm_index;
    construct(forward_fm_index, forward_raw_sequence_file, 1);
    store_to_file(forward_fm_index, forward_fm_index_file);

    cout << "Constructing reverse genomic index." << endl;
    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> reverse_fm_index;
    construct(reverse_fm_index, reverse_raw_sequence_file, 1);
    store_to_file(reverse_fm_index, reverse_fm_index_file);

    cout << "Index construction complete." << endl;
    return 0;
}

int do_enumerate_cmd(const enumerate_cmd_options& opts) {
    using namespace std;

    string forward_fm_index_file = opts.index_file_prefix + ".forward";
    string reverse_fm_index_file = opts.index_file_prefix + ".reverse";
    string genome_structure_file = opts.index_file_prefix + ".gs";

    cout << "Loading genome index..." << endl;

    genomics::genome_structure gs;
    if (!genomics::seq_io::load_from_file(gs, genome_structure_file)) {
        cerr << "No genome structure file \"" << genome_structure_file
             << "\" located. Abort." << endl;
        return 1;
    }

    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> forward_fm_index;
    if (!load_from_file(forward_fm_index, forward_fm_index_file)) {
        cerr << "No forward index file \"" << forward_fm_index_file
             << "\" located. Abort." << endl;
        return 1;
    }   

    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> reverse_fm_index;
    if (!load_from_file(reverse_fm_index, reverse_fm_index_file)) {
        cerr << "No reverse index file \"" << reverse_fm_index_file
             << "\" located. Abort." << endl;
        return 1;
    }   

    genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_forward(forward_fm_index, gs);
    genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_reverse(reverse_fm_index, gs);
    cout << "Successfully loaded indices." << endl;

    ofstream output(opts.database_file);
    genomics::write_sam_header(output, gi_forward.gs);

    genomics::kmers_file_producer kmer_p(opts.kmers_file);

    /* Require KMER file as input */

    std::mutex output_mtx;
    std::mutex kmer_mtx;

    std::vector<std::string> pams = opts.alt_pams;

    vector<thread> threads;
    for (size_t i = 0; i < opts.nthreads; i++) {
        thread t(genomics::process_kmers_to_stream<t_wt, t_sa_dens, t_isa_dens>,
                 cref(gi_forward), cref(gi_reverse),
                 cref(pams), opts.mismatches, opts.threshold,
		 ref(kmer_p), ref(kmer_mtx),
		 ref(output), ref(output_mtx));
        threads.push_back(move(t));
    }

    for (auto &thread : threads) {
        thread.join();
    }
 
    return 0;
}

int do_http_server_cmd(const http_server_cmd_options& opts) {
    using namespace std;
    using json = nlohmann::json;

    string forward_fm_index_file = opts.index_file_prefix + ".forward";
    string reverse_fm_index_file = opts.index_file_prefix + ".reverse";
    string genome_structure_file = opts.index_file_prefix + ".gs";

    cout << "Loading genome index..." << endl;

    genomics::genome_structure gs;
    if (!genomics::seq_io::load_from_file(gs, genome_structure_file)) {
        cerr << "No genome structure file \"" << genome_structure_file
             << "\" located. Abort." << endl;
        return 1;
    }

    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> forward_fm_index;
    if (!load_from_file(forward_fm_index, forward_fm_index_file)) {
        cerr << "No forward index file \"" << forward_fm_index_file
             << "\" located. Abort." << endl;
        return 1;
    }   

    sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> reverse_fm_index;
    if (!load_from_file(reverse_fm_index, reverse_fm_index_file)) {
        cerr << "No reverse index file \"" << reverse_fm_index_file
             << "\" located. Abort." << endl;
        return 1;
    }   

    genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_forward(forward_fm_index, gs);
    genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_reverse(reverse_fm_index, gs);
    cout << "Successfully loaded indices." << endl;

    httplib::Server svr;
    svr.Get("/search", [&gi_forward, &gi_reverse, &opts](const httplib::Request& req, httplib::Response& res){
        if (!req.has_param("sequence")) {
            return;
        }

        auto sequence = req.get_param_value("sequence");
        if (sequence.length() == 0) return;

        json result = search_kmer(gi_forward, gi_reverse, sequence, opts.mismatches, opts.ot_limit);
        res.set_content(result.dump(), "application/json");
    });

    cout << "Successfully started local server." << endl;
    svr.listen("localhost", opts.port);
        
    return 0;
}

int main(int argc, char *argv[])
{
    CLI::App guidescan("Guidescan all-in-one interface.\n");
    guidescan.require_subcommand(1);
    guidescan.failure_message(CLI::FailureMessage::help);

    enumerate_cmd_options enumerate_opts;
    index_cmd_options index_opts;
    http_server_cmd_options http_opts;

    auto enumerate = enumerate_cmd(guidescan, enumerate_opts);
    auto index = index_cmd(guidescan, index_opts);
    auto http  = http_cmd(guidescan, http_opts);

    (void) enumerate; (void) index; (void) http; // supress unused variable warnings

    try {
        guidescan.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return guidescan.exit(e);
    }
    
    if (guidescan.got_subcommand("index")) {
        return do_index_cmd(index_opts);
    }

    if (guidescan.got_subcommand("enumerate")) {
        return do_enumerate_cmd(enumerate_opts);
    }

    if (guidescan.got_subcommand("http-server")) {
        return do_http_server_cmd(http_opts);
    }

    return 1;
}
