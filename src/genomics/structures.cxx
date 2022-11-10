#include <assert.h>
#include <tuple>

#include "genomics/structures.hpp"

namespace genomics {
    std::tuple<coordinates, std::string> resolve_absolute(const genome_structure& gs, int64_t off_target_abs_coords, const kmer& k) {

      std::string strand("+");
      if (off_target_abs_coords < 0) {
        off_target_abs_coords = -off_target_abs_coords;
        strand = "-";
      }

      chromosome c = {"", 0};
      for (const chromosome& chr : gs) {
        if (off_target_abs_coords <= (int64_t)(chr.length-1)) {
          c = chr;
          break;
        } else {
          off_target_abs_coords -= chr.length;
        }
      }

      assert(c.name != "" && "Absolute coords longer than genome.");

      /* Start and end positions of match. start_position may turn out -ve so these are int64_t */
      int64_t start_position, end_position;
      if (strand == "+") {
        /* For + strand, absolute_coords denotes the 0-indexed endpoint (inclusive) of the match
           Convert this to 1-indexed (inclusive) to match samfile notation (1-indexed, both ends inclusive) */
        end_position = off_target_abs_coords + 1;
        /* The start index (1-indexed, inclusive) is end - len(seq_including_pam) + 1
         * Important to typecast operands here to preserve type! */
        start_position = end_position - (int64_t)k.sequence.length() - (int64_t)k.pam.length() + 1;
      } else {
        /* For - strand, absolute_coords denotes the 0-indexed startpoint (inclusive) of the match
           Convert this to 1-indexed (inclusive) to match samfile notation (1-indexed, both ends inclusive) */
        start_position = off_target_abs_coords + 1;
        /* The end index (1-indexed, inclusive) is start + len(seq_including_pam) - 1 */
        end_position = start_position + k.sequence.length() + k.pam.length() - 1;
      }

      /* If start_position < 0 or end_position > chr_length, we have gone beyond the chromosome boundary.
       * Return sentinel chromosome to indicate these (rare) cases. */
      if ((start_position < 0) || (end_position > (int64_t)c.length)) {
        return std::make_tuple(coordinates{chromosome{"", 0}, 0}, "");
      }

      /* For both + and - strand matches, we return the (lower in value) start_position, following convention */
      return std::make_tuple(coordinates{c, (uint64_t)start_position}, strand);
  }
}
