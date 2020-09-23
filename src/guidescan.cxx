#include <thread>
#include <istream>
#include <memory>

#include <sdsl/suffix_arrays.hpp>

#include "CLI/CLI.hpp"
#include "genomics/index.hpp"
#include "genomics/sam.hpp"
#include "genomics/seq_io.hpp"
#include "genomics/process.hpp"
#include "genomics/kmer.hpp"

struct build_cmd_options {
    size_t kmer_length;
    CLI::Option* kmer_length_opt;

    size_t mismatches;
    CLI::Option* mismatches_opt;
    
    size_t nthreads;
    CLI::Option* nthreads_opt;

    std::string fasta_file;
    CLI::Option* fasta_file_opt;

    std::string database_file;
    CLI::Option* database_file_opt;

    std::string kmers_file;
    CLI::Option* kmers_file_opt;

    std::string pam;
    CLI::Option* pam_opt;

    std::vector<std::string> alt_pams;
    CLI::Option* alt_pams_opt;
};

struct kmer_cmd_options {
    size_t kmer_length;
    CLI::Option* kmer_length_opt;

    std::string kmers_file;
    CLI::Option* kmers_file_opt;

    std::string fasta_file;
    CLI::Option* fasta_file_opt;

    std::string pam;
    CLI::Option* pam_opt;
};

CLI::App* build_cmd(CLI::App &guidescan, build_cmd_options& opts) {
    auto build = guidescan.add_subcommand("build", "Builds a gRNA database over the given genome.");

    opts.kmer_length = 20;
    opts.nthreads    = std::thread::hardware_concurrency();
    opts.pam         = std::string("NGG");
    opts.alt_pams    = {std::string("NAG")};
    opts.mismatches  = 3;

    opts.kmer_length_opt = build->add_option("-k,--kmer-length", opts.kmer_length, "Length of kmers excluding the PAM", true);
    opts.nthreads_opt    = build->add_option("-n,--threads", opts.nthreads, "Number of threads to parallelize over", true);
    opts.pam_opt         = build->add_option("-p,--pam", opts.pam, "Main PAM to generate gRNAs and find off-targets", true);
    opts.alt_pams_opt    = build->add_option("-a,--alt-pam", opts.alt_pams, "Alternative PAMs used to find off-targets", true);
    opts.mismatches_opt  = build->add_option("-m,--mismatches", opts.mismatches, "Number of mismatches to allow when finding off-targets", true);
    opts.kmers_file_opt  = build->add_option("-f,--kmers-file", opts.kmers_file,
					     "File containing kmers to build gRNA database"
					     " over, if not specified, will generate the database over all kmers with the given PAM")
	->check(CLI::ExistingFile);
    opts.fasta_file_opt  = build->add_option("genome", opts.fasta_file, "Genome in FASTA format")
	->check(CLI::ExistingFile)
	->required();
    opts.database_file_opt = build->add_option("-o, --output", opts.database_file, "Output database file.")
	->required();
  
    return build;
}

CLI::App* kmer_cmd(CLI::App &guidescan, kmer_cmd_options& opts) {
    auto kmers = guidescan.add_subcommand("kmers",
					  "Generates a list of kmers for a specific PAM written and"
					  " writes them to stdout.");
    opts.kmer_length = 20;
    opts.pam         = std::string("NGG");

    opts.kmer_length_opt = kmers->add_option("-k,--kmer-length", opts.kmer_length, "Length of kmers excluding the PAM", true);
    opts.pam_opt         = kmers->add_option("-p,--pam", opts.pam, "PAM to generate kmers for", true);
    opts.fasta_file_opt  = kmers->add_option("genome", opts.fasta_file, "Genome in FASTA format")
	->check(CLI::ExistingFile)
	->required();
    opts.kmers_file_opt  = kmers->add_option("-o, --output", opts.kmers_file, "Output kmers file.")
	->required();

    return kmers;
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

int do_build_cmd(const build_cmd_options& opts) {
    using namespace std;

    string raw_sequence_file = opts.fasta_file + ".dna";
    string genome_structure_file = opts.fasta_file + ".gs";
    string fm_index_file = opts.fasta_file + ".csa";
    
    ifstream fasta_is(opts.fasta_file);
    if (!fasta_is) {
        cerr << "ERROR: FASTA file \"" << opts.fasta_file
             << "\" does not exist." << endl;
        return 1;
    }

    cout << "Reading sequence file..." << endl;
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

    cout << "Loading genome index..." << endl;
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

    ofstream output(opts.database_file);
    genomics::write_sam_header(output, gi.gs);

    std::unique_ptr<genomics::kmer_producer> kmer_p;

    if (opts.kmers_file_opt->count() > 0) {
	kmer_p = make_unique<genomics::kmers_file_producer>(opts.kmers_file);
    } else {
	kmer_p = make_unique<genomics::seq_kmer_producer>(raw_sequence_file, opts.kmer_length, opts.pam);
    }

    std::mutex output_mtx;
    std::mutex kmer_mtx;

    std::vector<std::string> pams = opts.alt_pams;
    pams.push_back(opts.pam);

    vector<thread> threads;
    for (int i = 0; i < opts.nthreads; i++) {
        thread t(genomics::process_kmers_to_stream<sdsl::wt_huff<>, 32, 8192>,
                 cref(gi), cref(pams), opts.mismatches,
		 ref(kmer_p), ref(kmer_mtx),
		 ref(output), ref(output_mtx));
        threads.push_back(move(t));
    }

    for (auto &thread : threads) {
        thread.join();
    }
 
    return 0;
}

int do_kmers_cmd(const kmer_cmd_options& opts) {
    using namespace std;

    string raw_sequence_file = opts.fasta_file + ".dna";

    ifstream fasta_is(opts.fasta_file);
    if (!fasta_is) {
        cerr << "ERROR: FASTA file \"" << opts.fasta_file
             << "\" does not exist." << endl;
        return 1;
    }

    cout << "Reading sequence file..." << endl;
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

    cout << "Loading kmers from sequence..." << endl;
    genomics::seq_kmer_producer kmer_p(raw_sequence_file, opts.kmer_length, opts.pam);

    genomics::kmer k;
    vector<genomics::kmer> kmers;

    while(kmer_p.get_next_kmer(k)) {
	kmers.push_back(k);
    }
    
    cout << "Writing kmers to file..." << endl;
    genomics::seq_io::write_to_file(kmers, opts.kmers_file);

    return 0;
}
int main(int argc, char *argv[])
{
    CLI::App guidescan("Guidescan all-in-one interface.\n");
    guidescan.require_subcommand(1);
    guidescan.failure_message(CLI::FailureMessage::help);

    build_cmd_options build_opts;
    kmer_cmd_options kmer_opts;

    build_cmd(guidescan, build_opts);
    kmer_cmd(guidescan, kmer_opts);

    try {
	guidescan.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
	return guidescan.exit(e);
    }
    
    if (guidescan.got_subcommand("kmers")) {
	return do_kmers_cmd(kmer_opts);
    }

    if (guidescan.got_subcommand("build")) {
	return do_build_cmd(build_opts);
    }

    return 1;
}
