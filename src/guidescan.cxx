#include <thread>
#include <istream>
#include <memory>
#include <atomic>
#include <chrono>

#include <sdsl/suffix_arrays.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include "genomics/index.hpp"
#include "genomics/printer.hpp"
#include "genomics/seq_io.hpp"
#include "genomics/process.hpp"
#include "genomics/kmer.hpp"
#include "guidescan.hpp"
#include "version.hpp"

#define t_sa_dens 64
#define t_isa_dens 8192

typedef sdsl::wt_huff<> t_wt;
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
  opts.rna_bulges  = 0;
  opts.dna_bulges  = 0;
  opts.start       = false;
  opts.out_format  = "csv";
  opts.out_mode    = "complete";
  opts.max_off_targets = -1;

  opts.start_opt       = build->add_flag("--start", opts.start, "Match PAM at start of kmer instead at end (default).");
  opts.max_off_targets_opt = build->add_option("--max-off-targets", opts.max_off_targets, "Maximum number of off-targets to store for each number of mismatches.", true);
  opts.nthreads_opt    = build->add_option("-n,--threads", opts.nthreads, "Number of threads to parallelize over", true);
  opts.alt_pams_opt    = build->add_option("-a,--alt-pam", opts.alt_pams, "Alternative PAMs used to find off-targets", true);
  opts.mismatches_opt  = build->add_option("-m,--mismatches", opts.mismatches, "Number of mismatches to allow when finding off-targets", true);
  opts.mismatches_opt  = build->add_option("--rna-bulges", opts.rna_bulges, "Max number of RNA bulges to allow when finding off-targets", true);
  opts.mismatches_opt  = build->add_option("--dna-bulges", opts.dna_bulges, "Number of DNA bulges to allow when finding off-targets", true);
  opts.threshold_opt   = build->add_option("-t,--threshold", opts.threshold, "Filters gRNAs with off-targets at a distance at or below this threshold", true);
  opts.out_format_opt  = build->add_option("--format", opts.out_format, "File format for output. Choices are ['csv', 'sam'].")
    ->transform(CLI::IsMember({"csv", "sam"}, CLI::ignore_case));
  opts.out_mode_opt  = build->add_option("--mode", opts.out_mode, "Information to output. Choices are ['succinct', 'complete'].")
    ->transform(CLI::IsMember({"succinct", "complete"}, CLI::ignore_case));
  opts.kmers_file_opt  = build->add_option("-f,--kmers-file", opts.kmers_file,
                                           "File containing kmers to build gRNA database"
                                           " over, if not specified, will generate the database over all kmers with the given PAM")
    ->check(CLI::ExistingFile)
    ->required();
  opts.index_file_prefix_opt  = build->add_option("index", opts.index_file_prefix, "Prefix for index files")
    ->required();
  opts.database_file_opt = build->add_option("-o, --output", opts.database_file, "Output file.")
    ->required();
  
  return build;
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

  spdlog::info("Constructing genomic index.");

  genomics::genome_structure gs;
  fasta_is.clear();
  fasta_is.seekg(0);
  gs = genomics::seq_io::parse_genome_structure(fasta_is);
  genomics::seq_io::write_to_file(gs, genome_structure_file);

  spdlog::info("Constructing forward genomic index.");
  sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> forward_fm_index;
  construct(forward_fm_index, forward_raw_sequence_file, 1);
  store_to_file(forward_fm_index, forward_fm_index_file);

  spdlog::info("Constructing reverse genomic index.");
  sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> reverse_fm_index;
  construct(reverse_fm_index, reverse_raw_sequence_file, 1);
  store_to_file(reverse_fm_index, reverse_fm_index_file);

  spdlog::info("Index construction complete.");
  return 0;
}

int do_enumerate_cmd(const enumerate_cmd_options& opts) {
  using namespace std;

  const auto& error_log = spdlog::get("error");

  string forward_fm_index_file = opts.index_file_prefix + ".forward";
  string reverse_fm_index_file = opts.index_file_prefix + ".reverse";
  string genome_structure_file = opts.index_file_prefix + ".gs";

  spdlog::info("Loading genome index at \"{}\".", opts.index_file_prefix);

  std::chrono::steady_clock::time_point _start;
  std::chrono::steady_clock::time_point _end;

  _start = std::chrono::steady_clock::now();

  genomics::genome_structure gs;
  if (!genomics::seq_io::load_from_file(gs, genome_structure_file)) {
    error_log->error("No genome structure file {} located.", genome_structure_file);
    return 1;
  }

  sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> forward_fm_index;
  if (!load_from_file(forward_fm_index, forward_fm_index_file)) {
    error_log->error("No forward index file {} located.", forward_fm_index_file);
    return 1;
  }   

  sdsl::csa_wt<t_wt, t_sa_dens, t_isa_dens> reverse_fm_index;
  if (!load_from_file(reverse_fm_index, reverse_fm_index_file)) {
    error_log->error("No forward index file {} located.", reverse_fm_index_file);
    return 1;
  }
  _end = std::chrono::steady_clock::now();
  auto _elapsed = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start).count();
  spdlog::info("Loaded gs/forward/reverse in {} microseconds.", _elapsed);

  _start = std::chrono::steady_clock::now();
  genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_forward(forward_fm_index, gs);
  genomics::genome_index<t_wt, t_sa_dens, t_isa_dens> gi_reverse(reverse_fm_index, gs);
  spdlog::info("Successfully loaded genome index.");
  _end = std::chrono::steady_clock::now();
  _elapsed = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start).count();
  spdlog::info("Constructed structs in {} microseconds.", _elapsed);

  ofstream output(opts.database_file);
  bool complete = opts.out_mode == "complete";
  if (opts.out_format == "sam") {
    genomics::write_sam_header(output, gi_forward.gs);
  } else {
    genomics::write_csv_header(output, complete);
  }

  spdlog::info("Loading kmers.");
  genomics::kmers_file_producer kmer_p(opts.kmers_file);

  // split kmers across threads.
  vector<vector<genomics::kmer>> kmers(opts.nthreads, vector<genomics::kmer>());
  genomics::kmer out_kmer;
  int kmer_count = 0;
  for (; kmer_p.get_next_kmer(out_kmer); kmer_count++) {
    kmers[kmer_count % opts.nthreads].push_back(out_kmer);
  }

  spdlog::info("Read in {} kmer(s).", kmer_count);

  std::mutex output_mtx;
  std::vector<std::string> pams = opts.alt_pams;

  std::atomic<uint64_t> num_kmers_processed(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  vector<thread> threads;
  for (size_t i = 0; i < opts.nthreads; i++) {
    thread t(genomics::process_kmers_to_stream<t_wt, t_sa_dens, t_isa_dens>,
             cref(gi_forward), cref(gi_reverse),
             cref(opts), cref(kmers[i]),
             ref(output), ref(output_mtx), ref(num_kmers_processed), kmer_count, start_time, complete);
    threads.push_back(move(t));
  }

  for (auto &thread : threads) {
    thread.join();
  }

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
  spdlog::info("Processed {} kmers in {} seconds.", num_kmers_processed, elapsed_time);
 
  return 0;
}

int main(int argc, char *argv[])
{
  CLI::App guidescan("Guidescan all-in-one interface.\n");
  guidescan.set_version_flag("--version", std::string(GUIDESCAN_VERSION));

  guidescan.require_subcommand(1);
  guidescan.failure_message(CLI::FailureMessage::help);

  auto console_logger = spdlog::stdout_color_mt("guidescan2");
  spdlog::set_default_logger(console_logger);

  auto error_logger = spdlog::stderr_color_mt("error");

  enumerate_cmd_options enumerate_opts;
  index_cmd_options index_opts;

  auto enumerate = enumerate_cmd(guidescan, enumerate_opts);
  auto index = index_cmd(guidescan, index_opts);

  (void) enumerate; (void) index;; // supress unused variable warnings

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

  return 1;
}
