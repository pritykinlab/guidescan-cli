#ifndef PRINTER_H
#define PRINTER_H

#include "genomics/structures.hpp"
#include "genomics/sequences.hpp"
#include "genomics/doench.hpp"
#include "spdlog/fmt/fmt.h"
#include "version.hpp"

#include <algorithm>
#include <random>
#include <tuple>
#include <iostream>
#include <vector>

namespace genomics {
  namespace {
    std::string byte_to_big_endian_hex(unsigned char n) {
      std::string out("");
      int lowest1 = n % 16;
      n           = n / 16;
      int lowest2 = n % 16;

      char c1 = 'f';
      switch(lowest1) {
      case 0:  c1 = '0'; break;
      case 1:  c1 = '1'; break;
      case 2:  c1 = '2'; break;
      case 3:  c1 = '3'; break;
      case 4:  c1 = '4'; break;
      case 5:  c1 = '5'; break;
      case 6:  c1 = '6'; break;
      case 7:  c1 = '7'; break;
      case 8:  c1 = '8'; break;
      case 9:  c1 = '9'; break;
      case 10: c1 = 'a'; break;
      case 11: c1 = 'b'; break;
      case 12: c1 = 'c'; break;
      case 13: c1 = 'd'; break;
      case 14: c1 = 'e'; break;
      case 15: c1 = 'f'; break;
      };

      char c2 = 'f';
      switch(lowest2) {
      case 0:  c2 = '0'; break;
      case 1:  c2 = '1'; break;
      case 2:  c2 = '2'; break;
      case 3:  c2 = '3'; break;
      case 4:  c2 = '4'; break;
      case 5:  c2 = '5'; break;
      case 6:  c2 = '6'; break;
      case 7:  c2 = '7'; break;
      case 8:  c2 = '8'; break;
      case 9:  c2 = '9'; break;
      case 10: c2 = 'a'; break;
      case 11: c2 = 'b'; break;
      case 12: c2 = 'c'; break;
      case 13: c2 = 'd'; break;
      case 14: c2 = 'e'; break;
      case 15: c2 = 'f'; break;
      };

      out += c2;
      out += c1;

      return out;
    }

    std::string num_to_little_endian_hex(uint64_t num) {
      std::string out("");

      for (size_t i = 0; i < sizeof(num); i++) {
        unsigned char lowest = num % 256;
        num                  = num / 256;
        out += byte_to_big_endian_hex(lowest);
      }

      return out;
    }

    std::string list_to_little_endian_hex(std::list<int64_t> v) {
      std::string out("");
      for (const auto& n : v) {
        out += num_to_little_endian_hex(n);
      }
      return out;
    }

    int64_t get_delim(genome_structure gs) {
      int64_t delim = 0;
      for (const auto& chr : gs) {
        delim += chr.length;
      } 
      return -(delim + 1);
    }

    float calculate_cfd(std::string sgRNA, std::string sequence, std::string PAM) {
      if (PAM.length() != 3) return 1.0; // CFD only defined for length 3 PAMs.

      float cfd = 1.0;
      for (int i=0; i<sgRNA.length(); i++) {
        char _sgRNA = sgRNA.at(i);
        char _sequence = sequence.at(i);
        if (_sgRNA != _sequence) {
          if (_sgRNA == 'T') _sgRNA = 'U';
          auto key = fmt::format("r{}:d{},{}", _sgRNA, static_cast<char>(toupper(complement(_sequence))), i+1);
          cfd *= genomics::mm_scores[key];
        }
      }
      cfd *= genomics::pam_scores[PAM.substr(1,2)];
      return cfd;
    }

    std::tuple<std::string, float> off_target_fields(const kmer& k,
                                                     genome_structure gs,
                                                     std::vector<std::vector<std::tuple<int64_t, match>>>& off_targets,
                                                     int64_t max_off_targets) {
      uint64_t delim = get_delim(gs);
      std::string offtarget_hex("");
      float cfd_sum = 0.0;
      bool perfect_match = false;

      for (uint64_t i = 0; i < off_targets.size(); i++) {
          std::list<int64_t> v;

          int64_t n_off_targets = 0;
          for (const auto& t : off_targets[i]) {
              if ((max_off_targets != -1) && (n_off_targets >= max_off_targets))
                  break;
              auto pos = std::get<0>(t);
              auto match = std::get<1>(t);
              std::string match_sequence = complement(match.sequence);
              // Note that we cannot simply take k.pam as the PAM
              // because the match may well have been found using
              // an alternate PAM. The relevant PAM to consider for
              // CFD calculations is found at the end of the match.
              std::string pam = match_sequence.substr(20, 3);

              if ((match.mismatches == 0) && (pam.substr(1, 2) == "GG"))
                perfect_match = true;

              coordinates c;
              std::string strand;
              tie(c, strand) = resolve_absolute(gs, pos, k);
              if (c.chr.name == "") continue;
              v.push_back(pos);

              cfd_sum += calculate_cfd(k.sequence, match_sequence, pam);

              n_off_targets++;
          }

          v.push_back(i);
          v.push_back(delim);

          offtarget_hex += list_to_little_endian_hex(v);
      }

      float specificity = 0.0;
      if (!perfect_match) cfd_sum += 1;
      if (cfd_sum > 0) specificity = 1/cfd_sum;

      return std::make_tuple(offtarget_hex, specificity);
    }
  };

  void write_sam_header(std::ostream& os, genome_structure gs) {
    os << "@HD\tVN:1.0\tSO:unknown" << std::endl;
    os << "@PG\tID:Guidescan\tVN:" << std::string(GUIDESCAN_VERSION) << std::endl;
    for (const auto& chr : gs) {
      os << "@SQ\tSN:" << chr.name << "\tLN:" << chr.length << std::endl;
    }
  }

  void write_csv_header(std::ostream& os, bool complete) {
    os << "id,sequence,match_chrm,match_position,match_strand,match_distance";
    if (complete) {
      os << ",match_sequence,rna_bulges,dna_bulges";
    }
    os << ",specificity" << std::endl;
  }

    template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
    std::string get_csv_line_nomatch(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                             const kmer& k, bool start, bool complete) {

      std::string csvline(k.id);
      std::string sequence = start ? k.pam + k.sequence : k.sequence + k.pam;
      csvline += "," + sequence + ",NA,NA,NA,0,1.0";  // sequence, chrom, pos, strand, distance, specificity
      if (complete) csvline += ",NA,NA,NA";  // match sequence, rna bulges, dna bulges;
      return csvline;
    }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  std::string get_csv_line(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                           const kmer& k, bool start, match m,
                           int64_t off_target_abs_coords, bool complete) {

    coordinates c;
    std::string strand;
    tie(c, strand) = resolve_absolute(gi.gs, off_target_abs_coords, k);

    if (c.chr.name=="") return std::string("");  // Handle sentinel value for boundary conditions

    std::string sequence = start ? k.pam + k.sequence : k.sequence + k.pam;
    std::string csvline(k.id);

    csvline += ",";
    csvline += sequence;

    csvline += ",";
    csvline += c.chr.name;

    csvline += ",";
    csvline += std::to_string(c.offset);

    csvline += ",";
    csvline += strand;

    csvline += ",";
    csvline += std::to_string(m.mismatches);

    if (complete) {
      csvline += ",";
      csvline += complement(m.sequence);

      csvline += ",";
      csvline += std::to_string(m.rna_bulges);

      csvline += ",";
      csvline += std::to_string(m.dna_bulges);
    }

    return csvline;
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  std::string get_csv_lines(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                            const kmer& k, bool start,
                            int64_t max_off_targets,
                            std::vector<std::vector<std::tuple<int64_t, match>>>& off_targets, bool complete) {
    std::string csvlines;

    float cfd_sum = 0.0;
    bool perfect_match = false;
    bool no_match_no_offtargets = true;
    std::vector<std::string> off_targets_lines;

    for (uint64_t d = 0; d < off_targets.size(); d++) {
      for (int64_t i = 0; i < off_targets[d].size(); i++) {
        no_match_no_offtargets = false;
        if ((max_off_targets != -1) && (i >= max_off_targets)) break;

        const auto& off_target = off_targets[d][i];
        int64_t off_target_abs_coords = std::get<0>(off_target);
        match match = std::get<1>(off_target);
        std::string match_sequence = complement(match.sequence);
        // Note that we cannot simply take k.pam as the PAM
        // because the match may well have been found using
        // an alternate PAM. The relevant PAM to consider for
        // CFD calculations is found at the end of the match.
        std::string pam = match_sequence.substr(20, 3);

        if ((match.mismatches == 0) && (pam.substr(1, 2) == "GG"))
          perfect_match = true;

        std::string csvline = get_csv_line(gi, k, start, match, off_target_abs_coords, complete);
        if (csvline != "") {
          off_targets_lines.push_back(csvline);
          cfd_sum += calculate_cfd(k.sequence, match_sequence, pam);
        }
      }
    }

    if (no_match_no_offtargets) {
      csvlines += get_csv_line_nomatch(gi, k, start, complete) + "\n";
    } else {
      // Each off-target will have the same specificity value - add at the end of each line
      float specificity = 0.0;
      if (!perfect_match) cfd_sum += 1;
      if (cfd_sum > 0) specificity = 1 / cfd_sum;
      for (auto line: off_targets_lines) {
        csvlines += line + "," + std::to_string(specificity) + "\n";
      }
    }

    return csvlines;
  }

  template <class t_wt, uint32_t t_dens, uint32_t t_inv_dens>
  std::string get_sam_line(const genome_index<t_wt, t_dens, t_inv_dens>& gi,
                           const kmer& k, bool start, int64_t max_off_targets,
                           std::vector<std::vector<std::tuple<int64_t, match>>>& off_targets,
                           bool complete) {
    std::string sequence = start ? k.pam + k.sequence : k.sequence + k.pam;
    std::string samline(k.id);                                     // SAM:QNAME

    samline += "\t";
    samline += k.dir == direction::positive ? "0" : "16";          // SAM:FLAG (16 => SEQ being reverse complemented)
    samline += "\t" + k.chromosome;                                // SAM:RNAME

    samline += "\t" + std::to_string(k.position + 1);              // SAM:POS (1-indexed)
    samline += "\t100";                                            // SAM:MAPQ
    samline += "\t" + std::to_string(sequence.length()) + "M";     // SAM:CIGAR
    samline += "\t*\t0\t0";                                        // SAM:RNEXT, SAM:PNEXT, SAM:TLEN

    if (k.dir == direction::negative) {
      samline += "\t" + reverse_complement(sequence);              // SAM:SEQ
    } else {
      samline += "\t" + sequence;                                  // SAM:SEQ
    }

    samline += "\t*";                                              // SAM:QUAL
    
    // store off-target counts in tag
    for (uint64_t k = 0; k < off_targets.size(); k++) {
        samline += "\tk" + std::to_string(k) + ":i:" + std::to_string(off_targets[k].size());
    }

    std::string offtarget_hex;
    float specificity;
    tie(offtarget_hex, specificity) = off_target_fields(k, gi.gs, off_targets, max_off_targets);

    if (complete) {
        samline += "\tof:H:" + offtarget_hex;
    }
    samline += "\tsp:f:" + std::to_string(specificity);

    return samline;
  }
};

#endif /* PRINTER_H */
