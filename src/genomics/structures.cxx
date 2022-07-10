#include <assert.h>

#include "genomics/structures.hpp"

namespace genomics {
  coordinates resolve_absolute(const genome_structure& gs, uint64_t absolute_coords) {
    for (const chromosome& chr : gs) {
      if (absolute_coords < chr.length) {
        coordinates c = { chr, absolute_coords };
        return c;
      }

      absolute_coords -= chr.length;
    }

    assert(0 && "Absolute coords longer than genome.");

    return *((coordinates *) ((void *) nullptr)); // to suppress warning
  }


  uint64_t resolve_relative(const genome_structure& gs, coordinates coords) {
    uint64_t absolute_coords = 0;
    for (const chromosome& chr : gs) {
      if (chr == coords.chr) break;
      absolute_coords += chr.length;
    }

    return absolute_coords;
  }
}
