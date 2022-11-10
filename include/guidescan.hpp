#ifndef GUIDESCAN_H
#define GUIDESCAN_H

#include "CLI/CLI.hpp"

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

  size_t rna_bulges;
  CLI::Option* rna_bulges_opt = nullptr;
    
  size_t dna_bulges;
  CLI::Option* dna_bulges_opt = nullptr;

  std::string out_format;
  CLI::Option* out_format_opt = nullptr;

  std::string out_mode;
  CLI::Option* out_mode_opt = nullptr;

  size_t nthreads;
  CLI::Option* nthreads_opt = nullptr;

  int64_t max_off_targets;
  CLI::Option* max_off_targets_opt = nullptr;

  bool start;
  CLI::Option* start_opt = nullptr;

  std::string index_file_prefix;
  CLI::Option* index_file_prefix_opt = nullptr;

  std::string database_file;
  CLI::Option* database_file_opt = nullptr;

  std::string kmers_file;
  CLI::Option* kmers_file_opt = nullptr;

  std::vector<std::string> alt_pams;
  CLI::Option* alt_pams_opt = nullptr;
};

#endif
