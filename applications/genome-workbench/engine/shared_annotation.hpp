// The annotation of a sequence from its letters alone, as the Genome view (genome.annotate) and
// the toolbox (tool.annot.*) both compute it: genes by the domain's self-trained gene finder under
// a genetic code; CpG islands by the two-state HMM of Durbin et al., kept from 200 bases (the
// length of Gardiner-Garden and Frommer's criterion, 1987); maximal repeats on both strands by the
// suffix array and its LCP.

#pragma once

#include "genomes.hpp"

#include <genomics/annotation/genes.hpp>
#include <genomics/annotation/repeats.hpp>
#include <genomics/hmm/cpg.hpp>

#include <cstddef>
#include <vector>

namespace workbench::detail {

inline std::vector<genomics::annotation::Gene> predicted_genes(const seq::Codes& x, int table, genomics::annotation::GeneReport* report = nullptr) {
  genomics::annotation::GeneOptions o;
  o.table = table;
  return genomics::annotation::find_genes(x, o, report);
}

inline std::vector<genomics::hmm::Island> cpg_islands_of(const seq::Codes& x) {
  genomics::hmm::IslandOptions o;
  o.min_length = 200;
  return genomics::hmm::cpg_islands(x, o);
}

inline std::vector<genomics::annotation::Repeat> repeats_of(const seq::Codes& x, std::size_t min_length) {
  genomics::annotation::RepeatOptions o;
  o.min_length = min_length;
  o.both_strands = true;
  return genomics::annotation::maximal_repeats(x, o);
}

}  // namespace workbench::detail
