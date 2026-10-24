// Transcriptomes and RNA-seq fragments: transcripts with the genes they belong to, read from
// FASTA; the distribution of fragment lengths and the effective length of a transcript; and a
// simulator of transcriptomes (genes of exons, isoforms sharing them, paralogous genes) and of
// RNA-seq reads drawn from known abundances, with their truth.
//
// - Effective length (Pachter 2011, "Models for transcript quantification from RNA-Seq",
//   arXiv:1104.3889; Bray, Pimentel, Melsted and Pachter 2016, Nature Biotechnology 34,
//   525–527, Online Methods): a transcript of length l holds l − f + 1 fragments of length f, so
//   with fragment lengths of law p the expected number of places a fragment can come from is
//   Σ_{f ≤ l} p(f)(l − f + 1) = P(F ≤ l) · (l − μ_l + 1), μ_l the mean of p truncated to [1, l].
//   kallisto's effective length is l̃ = l − μ_l + 1, and l itself when that falls below 1 (or
//   when no fragment fits).
// - The model of a sample (Pachter 2011): ρ_t is the fraction of the molecules that are
//   copies of t (TPM = 10⁶ ρ); a fragment comes from t with probability α_t = ρ_t l̃_t /
//   Σ_u ρ_u l̃_u, has a length f of law p truncated to [1, l_t], and starts uniformly in
//   [0, l_t − f]. The simulator draws exactly this (it leaves out the factor P(F ≤ l_t), as the
//   likelihood of quantify.hpp does; it is 1 for transcripts longer than every fragment). A read
//   is the first `length` bases of the fragment on its own strand, or of its reverse complement
//   (unstranded, each with probability one half); a pair is both, the second mate from the other
//   end on the other strand. Each base is miscalled with probability `error`. Fragment i comes from
//   its own stream, so a shard of a sample is two numbers.
// - Simulated genes: a gene is a run of exons (random sequence); its first isoform has every
//   exon, the others random subsets (an internal exon skipped with probability 0.3, the first or
//   last with 0.1, at least two exons, every isoform distinct), so that isoforms share exons and
//   the k-mers of shared exons are ambiguous between them. A fraction of the genes are paralogs:
//   copies of an earlier gene's exons with each base substituted with probability
//   `paralog_divergence`, so that some k-mers are shared across genes too.
//
// Everything is in memory; a simulated transcriptome of n genes of mean length L takes O(n L).

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::expression {

using sequence::Code;
using sequence::Codes;
using sequence::Problem;

namespace detail {

// A base of GC content gc, and a base other than b, uniformly (as simulate/genome.hpp draws them).
inline Code draw_base(calx::random::Stream& r, double gc) {
  const double u = r.uniform();
  return u < gc / 2 ? sequence::G : u < gc ? sequence::C : u < (1 + gc) / 2 ? sequence::A : sequence::T;
}
inline Code substitute(calx::random::Stream& r, Code b) {
  const Code x = static_cast<Code>(r.below(3));
  return b > 3 ? x : static_cast<Code>(x >= b ? x + 1 : x);
}

}  // namespace detail

// The purposes of the streams of this module.
inline constexpr std::uint32_t TRANSCRIPTOME = 0x5452414e, ABUNDANCE = 0x4142554e, FRAGMENTS = 0x46524147, BOOTSTRAP = 0x424f4f54, COUNTS = 0x434f554e;

struct Transcript {
  std::string name, gene;
  Codes sequence;
};

struct Transcriptome {
  std::vector<Transcript> transcripts;
  std::vector<std::string> genes;        // the distinct genes, in order of first appearance
  std::vector<std::uint32_t> gene_of;    // per transcript, its gene's index

  std::size_t size() const { return transcripts.size(); }
  std::vector<std::size_t> lengths() const {
    std::vector<std::size_t> out;
    out.reserve(size());
    for (const auto& t : transcripts) out.push_back(t.sequence.size());
    return out;
  }
};

// The transcriptome of a list of transcripts: genes numbered in order of first appearance.
inline Transcriptome make_transcriptome(std::vector<Transcript> transcripts) {
  Transcriptome out;
  out.transcripts = std::move(transcripts);
  for (const auto& t : out.transcripts) {
    std::size_t g = 0;
    while (g < out.genes.size() && out.genes[g] != t.gene) g++;
    if (g == out.genes.size()) out.genes.push_back(t.gene);
    out.gene_of.push_back(static_cast<std::uint32_t>(g));
  }
  return out;
}

// Transcripts from FASTA: the name is the record's identifier, the gene the value of a
// "gene=" or "gene:" field of the description (GENCODE's and Ensembl's headers), else the name.
inline Transcriptome read_transcriptome(std::string_view fasta, std::vector<Problem>* problems = nullptr) {
  const auto records = sequence::read_fasta(fasta, problems);
  std::vector<Transcript> ts;
  for (const auto& r : records) {
    Transcript t;
    t.name = r.id;
    t.sequence = sequence::encode_dna(r.sequence);
    t.gene = r.id;
    for (const char* key : {"gene=", "gene:"}) {
      const std::size_t at = r.description.find(key);
      if (at == std::string::npos) continue;
      const std::size_t b = at + 5, e = r.description.find_first_of(" \t;", b);
      t.gene = r.description.substr(b, e == std::string::npos ? std::string::npos : e - b);
      break;
    }
    ts.push_back(std::move(t));
  }
  return make_transcriptome(std::move(ts));
}

// ---- fragment lengths and effective lengths ------------------------------------------------------------------------------------
struct FragmentLengths {
  std::vector<double> p;  // p[f] for f = 0 .. max, Σ = 1, p[0] = 0

  std::size_t max() const { return p.empty() ? 0 : p.size() - 1; }
  double mean() const {
    double m = 0;
    for (std::size_t f = 0; f < p.size(); f++) m += static_cast<double>(f) * p[f];
    return m;
  }
  // A normal law of the given mean and standard deviation, on the integers 1 .. max (by
  // default the mean plus six deviations), renormalised.
  static FragmentLengths normal(double mean, double sd, std::size_t max = 0) {
    if (!(mean > 0) || !(sd > 0)) throw std::invalid_argument("genomics::expression::FragmentLengths::normal: a positive mean and deviation");
    if (max == 0) max = static_cast<std::size_t>(std::ceil(mean + 6 * sd));
    FragmentLengths d;
    d.p.assign(max + 1, 0.0);
    double total = 0;
    for (std::size_t f = 1; f <= max; f++) {
      const double z = (static_cast<double>(f) - mean) / sd;
      total += d.p[f] = std::exp(-0.5 * z * z);
    }
    if (!(total > 0)) throw std::invalid_argument("genomics::expression::FragmentLengths::normal: no mass on 1 .. max");
    for (double& x : d.p) x /= total;
    return d;
  }
};

// kallisto's effective lengths, l − μ_l + 1, or l where that is below 1.
inline std::vector<double> effective_lengths(const std::vector<std::size_t>& lengths, const FragmentLengths& d) {
  // Cumulative mass and first moment, so that each transcript costs O(1).
  std::vector<double> mass(d.p.size() + 1, 0.0), moment(d.p.size() + 1, 0.0);
  for (std::size_t f = 0; f < d.p.size(); f++) {
    mass[f + 1] = mass[f] + d.p[f];
    moment[f + 1] = moment[f] + static_cast<double>(f) * d.p[f];
  }
  std::vector<double> out(lengths.size());
  for (std::size_t t = 0; t < lengths.size(); t++) {
    const std::size_t l = lengths[t], top = std::min(l, d.max());
    const double m = mass[top + 1];
    const double eff = m > 0 ? static_cast<double>(l) - moment[top + 1] / m + 1 : 0;
    out[t] = eff >= 1 ? eff : static_cast<double>(l);
  }
  return out;
}

// ---- simulated transcriptomes ------------------------------------------------------------------------------------
struct TranscriptomeSpec {
  std::size_t genes = 200;
  std::size_t min_exons = 2, max_exons = 10;
  std::size_t min_exon = 60, max_exon = 400;  // exon lengths, uniform
  std::size_t max_isoforms = 4;               // per gene, uniform from 1
  double gc = 0.5;
  double paralogs = 0.0;                      // the fraction of genes copied from an earlier one
  double paralog_divergence = 0.05;
};

struct SimulatedTranscriptome {
  Transcriptome transcriptome;
  std::vector<std::vector<std::uint32_t>> exons;  // per transcript, the exons of its gene it holds
};

inline SimulatedTranscriptome simulate_transcriptome(const TranscriptomeSpec& spec, std::uint64_t seed) {
  if (spec.min_exons < 1 || spec.max_exons < spec.min_exons || spec.min_exon < 1 || spec.max_exon < spec.min_exon || spec.max_isoforms < 1)
    throw std::invalid_argument("genomics::expression::simulate_transcriptome: inconsistent ranges");
  SimulatedTranscriptome out;
  std::vector<Transcript> ts;
  std::vector<std::vector<Codes>> gene_exons;
  for (std::size_t g = 0; g < spec.genes; g++) {
    calx::random::Stream r(seed, TRANSCRIPTOME, static_cast<std::uint32_t>(g), 0);
    std::vector<Codes> exons;
    if (g > 0 && r.uniform() < spec.paralogs) {
      exons = gene_exons[r.below(g)];
      for (auto& e : exons)
        for (auto& b : e)
          if (r.uniform() < spec.paralog_divergence) b = detail::substitute(r, b);
    } else {
      const std::size_t E = spec.min_exons + r.below(spec.max_exons - spec.min_exons + 1);
      for (std::size_t e = 0; e < E; e++) {
        Codes x(spec.min_exon + r.below(spec.max_exon - spec.min_exon + 1));
        for (auto& b : x) b = detail::draw_base(r, spec.gc);
        exons.push_back(std::move(x));
      }
    }
    const std::size_t E = exons.size();
    const std::size_t want = 1 + r.below(spec.max_isoforms);
    std::vector<std::vector<std::uint32_t>> isoforms;
    std::vector<std::uint32_t> all(E);
    for (std::size_t e = 0; e < E; e++) all[e] = static_cast<std::uint32_t>(e);
    isoforms.push_back(all);
    for (int attempt = 0; attempt < 50 && isoforms.size() < want && E >= 3; attempt++) {
      std::vector<std::uint32_t> keep;
      for (std::size_t e = 0; e < E; e++) {
        const bool end = e == 0 || e + 1 == E;
        if (r.uniform() >= (end ? 0.1 : 0.3)) keep.push_back(static_cast<std::uint32_t>(e));
      }
      if (keep.size() < 2 || std::find(isoforms.begin(), isoforms.end(), keep) != isoforms.end()) continue;
      isoforms.push_back(keep);
    }
    const std::string gene = "g" + std::to_string(g + 1);
    for (std::size_t k = 0; k < isoforms.size(); k++) {
      Transcript t;
      t.gene = gene;
      t.name = gene + ".t" + std::to_string(k + 1);
      for (std::uint32_t e : isoforms[k]) t.sequence.insert(t.sequence.end(), exons[e].begin(), exons[e].end());
      ts.push_back(std::move(t));
      out.exons.push_back(isoforms[k]);
    }
    gene_exons.push_back(std::move(exons));
  }
  out.transcriptome = make_transcriptome(std::move(ts));
  return out;
}

// Molar fractions ρ (summing to 1) from a log-normal law of standard deviation `sigma`, with a
// fraction `zero` of the transcripts not expressed.
inline std::vector<double> simulate_abundances(std::size_t transcripts, double sigma, double zero, std::uint64_t seed) {
  std::vector<double> rho(transcripts, 0.0);
  double total = 0;
  for (std::size_t t = 0; t < transcripts; t++) {
    calx::random::Stream r(seed, ABUNDANCE, static_cast<std::uint32_t>(t), 0);
    if (r.uniform() < zero) continue;
    total += rho[t] = r.lognormal(0, sigma);
  }
  if (total > 0)
    for (double& x : rho) x /= total;
  return rho;
}

// ---- simulated reads ------------------------------------------------------------------------------------
struct ReadModel {
  std::size_t length = 75;
  bool paired = false;
  double error = 0.0;  // per base, a substitution
};

struct Fragment {
  std::uint32_t transcript = 0;
  std::uint32_t begin = 0, length = 0;  // on the transcript
  bool reverse = false;                 // the first read from the reverse strand
  Codes read1, read2;                   // read2 empty unless paired
};

// The fractions of fragments α_t = ρ_t l̃_t / Σ ρ_u l̃_u of the model above.
inline std::vector<double> fragment_fractions(const std::vector<double>& rho, const std::vector<double>& effective) {
  std::vector<double> a(rho.size(), 0.0);
  double total = 0;
  for (std::size_t t = 0; t < rho.size(); t++) total += a[t] = rho[t] * effective[t];
  if (total > 0)
    for (double& x : a) x /= total;
  return a;
}

// Fragments first .. first + count − 1 of a sample of the transcriptome at molar fractions ρ.
inline std::vector<Fragment> simulate_fragments(const Transcriptome& tx, const std::vector<double>& rho, const FragmentLengths& d, const ReadModel& model,
                                                std::size_t count, std::uint64_t seed, std::size_t first = 0) {
  if (rho.size() != tx.size()) throw std::invalid_argument("genomics::expression::simulate_fragments: one abundance per transcript");
  const auto lengths = tx.lengths();
  const auto eff = effective_lengths(lengths, d);
  // Transcripts no fragment fits in are never drawn.
  std::vector<double> weight(tx.size()), cumulative(tx.size());
  double total = 0;
  for (std::size_t t = 0; t < tx.size(); t++) {
    const bool fits = lengths[t] >= 1 && std::any_of(d.p.begin() + 1, d.p.begin() + static_cast<std::ptrdiff_t>(std::min(lengths[t], d.max()) + 1), [](double x) { return x > 0; });
    weight[t] = fits ? rho[t] * eff[t] : 0.0;
    cumulative[t] = total += weight[t];
  }
  if (!(total > 0)) throw std::invalid_argument("genomics::expression::simulate_fragments: nothing to draw from");
  std::vector<double> mass(d.p.size() + 1, 0.0);
  for (std::size_t f = 0; f < d.p.size(); f++) mass[f + 1] = mass[f] + d.p[f];
  std::vector<Fragment> out;
  out.reserve(count);
  for (std::size_t i = first; i < first + count; i++) {
    calx::random::Stream r(seed, FRAGMENTS, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i >> 32));
    Fragment fr;
    const double u = r.uniform() * total;
    // The first transcript whose cumulative weight exceeds u, which has a weight of its own.
    fr.transcript = static_cast<std::uint32_t>(std::min<std::size_t>(tx.size() - 1, static_cast<std::size_t>(std::upper_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin())));
    const std::size_t l = lengths[fr.transcript], top = std::min(l, d.max());
    // A length from p truncated to [1, l], by inversion.
    const double v = r.uniform() * mass[top + 1];
    std::size_t f = static_cast<std::size_t>(std::upper_bound(mass.begin() + 1, mass.begin() + static_cast<std::ptrdiff_t>(top + 2), v) - mass.begin()) - 1;
    f = std::clamp<std::size_t>(f, 1, top);
    fr.length = static_cast<std::uint32_t>(f);
    fr.begin = static_cast<std::uint32_t>(r.below(l - f + 1));
    fr.reverse = r.bernoulli(0.5);
    const Codes& s = tx.transcripts[fr.transcript].sequence;
    const Codes piece(s.begin() + fr.begin, s.begin() + fr.begin + f);
    const Codes rc = sequence::reverse_complement(piece);
    const std::size_t L = std::min(model.length, f);
    const Codes& one = fr.reverse ? rc : piece;
    const Codes& two = fr.reverse ? piece : rc;
    fr.read1.assign(one.begin(), one.begin() + static_cast<std::ptrdiff_t>(L));
    if (model.paired) fr.read2.assign(two.begin(), two.begin() + static_cast<std::ptrdiff_t>(L));
    for (Codes* read : {&fr.read1, &fr.read2})
      for (auto& b : *read)
        if (model.error > 0 && r.uniform() < model.error) b = detail::substitute(r, b);
    out.push_back(std::move(fr));
  }
  return out;
}

}  // namespace genomics::expression
