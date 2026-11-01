// A test of prediction: how well the indices of codon usage of Escherichia coli K-12's genes rank
// their proteins' measured abundances (PaxDb 6.0's integrated dataset, data/registry/paxdb-ecoli.json).
// The thresholds are ours, stated with the measured values: codon usage is one of several things
// that set a protein's abundance (transcription, translation initiation, decay), so a substantial
// rank correlation, not a perfect one, is what it can show.

#include "common.hpp"

#include <genomics/codon/abundance.hpp>
#include <genomics/codon/adaptation.hpp>
#include <genomics/codon/trna.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <cstdio>
#include <string>

using namespace genomics;
using namespace genomics::codon;

namespace {
std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  // The reader, on a small text.
  {
    const auto t = read_paxdb("#name: test\n#coverage: 50\n#gene_name\tstring_external_id\tabundance\nuspA\t511145.b3495\t30087\nrplL\t511145.b3986\t11544.5\nbad\t511145.b0001\tx\n");
    s.require("PaxDb's format read: two rows, the locus after the taxon, a bad number reported with its line",
              t.rows.size() == 2 && t.by_locus.at("b3986") == 11544.5 && t.name == "test" && t.problems.size() == 1 && t.problems[0].line == 6);
  }
  const auto ecoli = check::genome("NC_000913.3");
  const auto text = check::read_file(check::cache() / "paxdb-ecoli" / "511145-WHOLE_ORGANISM-integrated.txt");
  if (!ecoli || !text) {
    s.skip("E. coli's genome or PaxDb's abundances are not in the data cache (ncbi-genomes, paxdb-ecoli): the indices are not tested against measured abundance");
    return s.done();
  }
  const auto table = read_paxdb(*text);
  s.require("PaxDb 6.0's E. coli table read without problems", table.problems.empty() && table.rows.size() > 3000, std::to_string(table.rows.size()) + " proteins, coverage " + table.coverage);
  const auto& code = sequence::genetic_code(11);
  const auto set = coding_genes(*ecoli);
  const auto rp = genes_by_product(set.genes, RIBOSOMAL_PROTEINS);
  std::vector<CodonCounts> counts;
  for (const auto& g : set.genes) counts.push_back(codon_counts(g));
  const auto w = relative_adaptiveness(pooled_counts(set.genes, rp), code);
  const auto pool = trna_pool(trna_genes(*ecoli));
  const auto ta = trna_adaptiveness(pool, code);
  const auto optimal = optimal_by_adaptiveness(w);
  const auto reference = pooled_counts(set.genes, rp), all = pooled_counts(set.genes);
  ExpressionStandards st;
  st.all = all;
  st.ribosomal = reference;
  st.chaperones = pooled_counts(set.genes, genes_by_name(set.genes, ecoli_chaperones()));
  st.factors = pooled_counts(set.genes, genes_by_name(set.genes, ecoli_translation_factors()));
  std::vector<double> v_cai, v_tai, v_fop, v_nc, v_melp, v_e, v_dcbs;
  for (const auto& c : counts) {
    v_cai.push_back(cai(c, w));
    v_tai.push_back(tai(c, ta));
    v_fop.push_back(fop(c, optimal, code));
    v_nc.push_back(enc_wright(c, code));
    v_melp.push_back(milc(c, all, code) / std::max(1e-12, milc(c, reference, code)));
    v_e.push_back(predicted_expression(c, st, code).e);
    v_dcbs.push_back(dcbs(c, code));
  }
  const auto cai_a = agreement(set.genes, v_cai, table);
  const auto tai_a = agreement(set.genes, v_tai, table);
  const auto fop_a = agreement(set.genes, v_fop, table);
  const auto nc_a = agreement(set.genes, v_nc, table);
  const auto melp_a = agreement(set.genes, v_melp, table);
  const auto e_a = agreement(set.genes, v_e, table);
  const auto dcbs_a = agreement(set.genes, v_dcbs, table);
  s.require("genes matched to PaxDb by locus tag: " + std::to_string(cai_a.matched) + " of " + std::to_string(set.genes.size()), cai_a.matched > 3500);
  s.require("CAI (ribosomal proteins) ranks protein abundance: Spearman ρ = " + fixed(cai_a.spearman) + " (Pearson r with log abundance " + fixed(cai_a.pearson_log) + "), above 0.4",
            cai_a.spearman > 0.4);
  s.require("tAI (E. coli's tRNA genes) ranks it too: ρ = " + fixed(tai_a.spearman) + ", above 0.25", tai_a.spearman > 0.25);
  s.require("Fop (CAI's optimal codons): ρ = " + fixed(fop_a.spearman) + "; MELP: ρ = " + fixed(melp_a.spearman) + "; Karlin and Mrázek's E: ρ = " + fixed(e_a.spearman) +
                "; DCBS: ρ = " + fixed(dcbs_a.spearman) + "; all positive",
            fop_a.spearman > 0.2 && melp_a.spearman > 0.2 && e_a.spearman > 0.2 && dcbs_a.spearman > 0);
  s.require("Wright's N_c, which falls as bias rises, is anti-correlated with abundance: ρ = " + fixed(nc_a.spearman), nc_a.spearman < -0.1);
  // A sanity check of the data: the ribosomal proteins are among the most abundant proteins.
  {
    std::vector<double> ab;
    for (const auto& [k, v] : table.by_locus) ab.push_back(v);
    const double median = tyche::descriptive::median(ab);
    std::size_t above = 0, n = 0;
    for (std::size_t i : rp) {
      const auto it = table.by_locus.find(set.genes[i].locus_tag);
      if (it == table.by_locus.end()) continue;
      n++;
      above += it->second > 10 * median;
    }
    s.require("the ribosomal proteins measured are, nearly all, over ten times the median protein's abundance", above * 10 >= n * 9, std::to_string(above) + " of " + std::to_string(n));
  }
  return s.done();
}
