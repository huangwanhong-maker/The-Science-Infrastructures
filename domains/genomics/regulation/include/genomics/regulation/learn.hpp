// Expression learnt from regulatory sequence: measured promoter activities read, families of near
// copies kept together, a penalised linear model on k-mer and gapped k-mer counts, a convolutional
// network trained on one-hot sequences, its first layer's filters read as position weight matrices,
// and those matched to known motifs.
//
// - Data: a measured activity per sequence (a massively parallel reporter assay, such as the
//   Random Promoter DREAM Challenge's; Rafi et al. 2025, Nature Biotechnology 43, 1373–1383; de Boer
//   et al. 2020, Nature Biotechnology 38, 56–65). Splits for held-out evaluation keep sequences that
//   share a key (here, two stretches of the variable region) on one side, so that a variant and the
//   sequence it varies are never one in training and one in test: random splits of such sets reward
//   memory (homology-induced leakage).
// - The linear baseline: counts of k-mers and gapped k-mers (encoding.hpp) as the columns of the
//   elastic net (tyche::regression's coordinate descent, Friedman, Hastie and Tibshirani 2010), λ by
//   K-fold cross-validation (the least mean error). Its coefficients name the words that raise or
//   lower activity.
// - The convolutional network (tyche::learning): one-hot input, a layer of filters (each a position
//   weight matrix in effect; Alipanahi et al. 2015), rectification, max pooling, dropout, a dense
//   layer and a linear output; trained by Adam on the squared error of the standardised activity.
//   The architecture is DeepBind's and Basset's in small (Kelley, Snoek and Rinn 2016, Genome Research
//   26, 990–999), far from the state of the art (attention and recurrent layers, millions of training
//   sequences: the DREAM Challenge's winners, Enformer, Borzoi, AlphaGenome), which the theory note
//   places it against.
// - Filters as motifs (Kelley et al. 2016, §"Filter interpretation"; Alipanahi et al. 2015): for each
//   filter, the windows of a set of sequences where its activation (the convolution before
//   rectification) exceeds half its largest, aligned and counted into a position frequency matrix,
//   with the filter's information content.
// - Matching to known motifs: a query matrix against each known one (both strands, every offset
//   overlapping at least `min_overlap` columns) by the mean Pearson correlation of overlapping
//   columns' probabilities (Pietrokovski 1996, Nucleic Acids Research 24, 3836–3845; one of Tomtom's
//   column scores, Gupta, Stamatoyannopoulos, Bailey and Noble 2007, Genome Biology 8, R24); its
//   significance empirical, against the best matches of the query with its columns shuffled.
// - Known motifs in MEME's format (JASPAR's downloads; the MEME Suite's motif format, version 4).
//
// Limits: one output (an activity); the network sees the sequence on its given strand only.

#pragma once

#include <genomics/annotation/motifs.hpp>
#include <genomics/regulation/encoding.hpp>
#include <genomics/sequence/fasta.hpp>

#include <tyche/descriptive/rank.hpp>
#include <tyche/descriptive/summary.hpp>
#include <tyche/learning/train.hpp>
#include <tyche/regression/penalized.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace genomics::regulation {

inline constexpr std::uint32_t GROUP_SPLIT = 0x47535054;  // "GSPT": the split of families
inline constexpr std::uint32_t SHUFFLE_COLUMNS = 0x53434f4c;  // "SCOL": the null of motif matching

// ---- measured activities ---------------------------------------------------------------------------------
struct Measured {
  std::vector<Codes> sequences;
  std::vector<double> values;
  std::vector<sequence::Problem> problems;
  std::size_t size() const { return values.size(); }
};

// Lines "sequence<TAB>value".
inline Measured read_measured(std::string_view text) {
  Measured m;
  sequence::detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (line.empty() || line[0] == '#') return;
    const auto tab = line.find('\t');
    if (tab == std::string_view::npos) { m.problems.push_back({number, "a line without a tab"}); return; }
    const std::string v(sequence::detail::trim(line.substr(tab + 1)));
    char* end = nullptr;
    const double value = std::strtod(v.c_str(), &end);
    if (v.empty() || *end || !std::isfinite(value)) { m.problems.push_back({number, "the value '" + v + "' is not a number"}); return; }
    m.sequences.push_back(sequence::encode_dna(line.substr(0, tab)));
    m.values.push_back(value);
  });
  return m;
}

inline Measured subset(const Measured& m, const std::vector<std::size_t>& rows) {
  Measured s;
  for (std::size_t r : rows) {
    s.sequences.push_back(m.sequences[r]);
    s.values.push_back(m.values[r]);
  }
  return s;
}

// A split that keeps sequences of one key together: the keys are shuffled by a stream and their
// sequences assigned to test, validation and training in turn until each part has its fraction.
struct Key {
  std::size_t from = 0, length = 0;
};
inline tyche::learning::Split grouped_split(const Measured& m, const std::vector<Key>& key, double validation, double test, std::uint64_t seed) {
  std::unordered_map<std::string, std::vector<std::size_t>> groups;
  std::vector<std::string> order;
  for (std::size_t i = 0; i < m.size(); i++) {
    std::string k;
    for (const auto& part : key)
      for (std::size_t t = part.from; t < std::min(m.sequences[i].size(), part.from + part.length); t++) k += sequence::dna_letter(m.sequences[i][t]);
    auto [it, fresh] = groups.try_emplace(k);
    if (fresh) order.push_back(k);
    it->second.push_back(i);
  }
  calx::random::Stream r(seed, GROUP_SPLIT, 0, 0);
  for (std::size_t i = order.size(); i > 1; i--) std::swap(order[i - 1], order[r.below(i)]);
  tyche::learning::Split s;
  const double n = static_cast<double>(m.size());
  for (const auto& k : order) {
    auto& part = static_cast<double>(s.test.size()) < test * n ? s.test : static_cast<double>(s.validation.size()) < validation * n ? s.validation : s.train;
    const auto& rows = groups[k];
    part.insert(part.end(), rows.begin(), rows.end());
  }
  for (auto* p : {&s.train, &s.validation, &s.test}) std::sort(p->begin(), p->end());
  return s;
}

struct Agreement {
  double pearson = 0, spearman = 0, r2 = 0;  // r2 = 1 − SSE/SST on the held-out values
  std::size_t n = 0;
};
inline Agreement agreement(const std::vector<double>& predicted, const std::vector<double>& measured) {
  Agreement a;
  a.n = measured.size();
  if (a.n < 3) return a;
  a.pearson = tyche::descriptive::correlation(predicted, measured);
  a.spearman = tyche::descriptive::spearman(predicted, measured);
  const double mean = tyche::descriptive::mean(measured);
  double sse = 0, sst = 0;
  for (std::size_t i = 0; i < a.n; i++) {
    sse += (predicted[i] - measured[i]) * (predicted[i] - measured[i]);
    sst += (measured[i] - mean) * (measured[i] - mean);
  }
  a.r2 = sst > 0 ? 1 - sse / sst : 0;
  return a;
}

// ---- the k-mer baseline ----------------------------------------------------------------------------------
struct KmerBaselineOptions {
  std::size_t k_min = 1, k_max = 4;
  std::size_t gapped_l = 6, gapped_k = 4;  // gapped_k = 0: none
  bool canonical = false;                  // promoters of a reporter assay have a strand
  double alpha = 0.5;                      // the elastic net's mixing
  std::size_t lambda_count = 20, folds = 5;
  double lambda_min_ratio = 1e-2;          // the path ends at this fraction of λ_max
  double tolerance = 1e-5;                 // of coordinate descent, relative to the null deviance
  std::size_t from = 0, length = 0;        // the region counted (length 0: the whole sequence)
  std::uint64_t seed = 1;
};

struct KmerBaseline {
  KmerBaselineOptions options;
  std::vector<std::string> names;          // the columns
  tyche::regression::PenalizedFit fit;
  double lambda = 0;
  std::vector<double> lambdas, cv_error;   // the cross-validation curve
  std::size_t nonzero = 0;
};

namespace detail {
inline calx::linalg::Rectangular kmer_matrix(const std::vector<Codes>& xs, const KmerBaselineOptions& o, std::vector<std::string>* names = nullptr) {
  const KmerFeatures km(o.k_min, o.k_max, o.canonical);
  std::unique_ptr<GappedKmerFeatures> gk;
  if (o.gapped_k) gk = std::make_unique<GappedKmerFeatures>(o.gapped_l, o.gapped_k, o.canonical);
  const std::size_t p = km.columns() + (gk ? gk->columns() : 0);
  calx::linalg::Rectangular X(xs.size(), p);
  for (std::size_t i = 0; i < xs.size(); i++) {
    const Codes& full = xs[i];
    const std::size_t a = std::min(o.from, full.size()), b = o.length ? std::min(full.size(), o.from + o.length) : full.size();
    const Codes x(full.begin() + static_cast<std::ptrdiff_t>(a), full.begin() + static_cast<std::ptrdiff_t>(b));
    double* row = &X.a[i * p];
    km.count(x, row);
    if (gk) gk->count(x, row + km.columns());
  }
  if (names) {
    *names = km.names();
    if (gk) names->insert(names->end(), gk->names().begin(), gk->names().end());
  }
  return X;
}
}  // namespace detail

inline KmerBaseline fit_kmer_baseline(const Measured& train, const KmerBaselineOptions& o = {}) {
  KmerBaseline b;
  b.options = o;
  const auto X = detail::kmer_matrix(train.sequences, o, &b.names);
  tyche::regression::PenalizedOptions po;
  po.alpha = o.alpha;
  po.lambda_count = o.lambda_count;
  po.lambda_min_ratio = o.lambda_min_ratio;
  po.tolerance = o.tolerance;
  tyche::regression::CvOptions cv;
  cv.folds = o.folds;
  cv.seed = o.seed;
  const auto r = tyche::regression::cross_validate(X, train.values, po, cv);
  b.fit = r.path.fits[r.index_min];
  b.lambda = r.lambda_min;
  b.lambdas = r.path.lambda;
  b.cv_error = r.mean_error;
  for (double v : b.fit.beta) b.nonzero += v != 0;
  return b;
}

inline std::vector<double> predict(const KmerBaseline& b, const std::vector<Codes>& xs) {
  return tyche::regression::predict(b.fit, detail::kmer_matrix(xs, b.options));
}

// ---- the convolutional network --------------------------------------------------------------------------
struct CnnOptions {
  std::size_t filters = 32, width = 15, pool = 4, hidden = 32;
  double dropout = 0.2;
  std::size_t epochs = 10, batch = 64, patience = 3;
  double rate = 1e-3, weight_decay = 0;
  std::uint64_t seed = 1;
};

struct SequenceCnn {
  CnnOptions options;
  std::size_t length = 0;
  double mean = 0, sd = 1;                      // the standardisation of the activity
  tyche::learning::Network<float> net{{4, 1}};
  tyche::learning::Training training;
};

namespace detail {
inline tyche::learning::Data<float> one_hot_data(const std::vector<Codes>& xs, const std::vector<double>& ys, std::size_t L, double mean, double sd) {
  tyche::learning::Data<float> d;
  d.input = {4, L};
  d.outputs = 1;
  d.x.resize(xs.size() * 4 * L);
  d.y.resize(xs.size());
  for (std::size_t i = 0; i < xs.size(); i++) {
    if (xs[i].size() != L) throw std::invalid_argument("genomics::regulation: sequences of one length");
    one_hot_into(xs[i], d.x.data() + i * 4 * L);
    d.y[i] = ys.empty() ? 0.0f : static_cast<float>((ys[i] - mean) / sd);
  }
  return d;
}
}  // namespace detail

inline SequenceCnn train_cnn(const Measured& train, const Measured& validation, const CnnOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  using namespace tyche::learning;
  if (train.size() == 0) throw std::invalid_argument("genomics::regulation::train_cnn: no training sequences");
  SequenceCnn m;
  m.options = o;
  m.length = train.sequences.front().size();
  m.mean = tyche::descriptive::mean(train.values);
  m.sd = tyche::descriptive::standard_deviation(train.values);
  if (!(m.sd > 0)) m.sd = 1;
  m.net = Network<float>({4, m.length});
  m.net.add<Conv1d<float>>(o.filters, o.width).add<Relu<float>>().add<MaxPool1d<float>>(o.pool);
  if (o.dropout > 0) m.net.add<Dropout<float>>(o.dropout);
  m.net.add<Dense<float>>(o.hidden).add<Relu<float>>().add<Dense<float>>(1, Init::glorot);
  m.net.initialise(o.seed);
  const auto tr = detail::one_hot_data(train.sequences, train.values, m.length, m.mean, m.sd);
  const auto va = detail::one_hot_data(validation.sequences, validation.values, m.length, m.mean, m.sd);
  TrainOptions to;
  to.epochs = o.epochs;
  to.batch = o.batch;
  to.patience = o.patience;
  to.adam.rate = o.rate;
  to.adam.weight_decay = o.weight_decay;
  to.seed = o.seed;
  m.training = tyche::learning::train(m.net, tr, va.size() ? &va : nullptr, to, pool);
  return m;
}

inline std::vector<double> predict(const SequenceCnn& m, const std::vector<Codes>& xs, calx::scheduler::ThreadPool* pool = nullptr) {
  const auto d = detail::one_hot_data(xs, {}, m.length, 0, 1);
  const auto y = tyche::learning::predict(m.net, d, pool);
  std::vector<double> out(y.size());
  for (std::size_t i = 0; i < y.size(); i++) out[i] = m.mean + m.sd * static_cast<double>(y[i]);
  return out;
}

// ---- filters as motifs --------------------------------------------------------------------------------------
struct FilterMotif {
  std::size_t filter = 0;
  annotation::Pwm pwm;                 // from the activating windows (counts, pseudocount 0.5)
  std::size_t sites = 0;               // windows above half the largest activation
  double information = 0;              // Σ of the columns' information, bits
  std::vector<double> weights;         // the filter itself, width × 4 (base fastest)
};

// The filters of the network's first layer (a Conv1d on one-hot input) as motifs, from the windows
// of `xs` where each fires above half its largest activation.
inline std::vector<FilterMotif> filter_motifs(const SequenceCnn& m, const std::vector<Codes>& xs) {
  const auto* conv = dynamic_cast<const tyche::learning::Conv1d<float>*>(&m.net.layer(0));
  if (!conv) throw std::invalid_argument("genomics::regulation::filter_motifs: the first layer is not a convolution");
  const std::size_t F = conv->filters(), W = conv->width();
  const auto theta = m.net.slice(0);
  const float* Wt = theta.data();
  const float* bias = theta.data() + F * 4 * W;
  std::vector<FilterMotif> out(F);
  // Largest activation of each filter over the windows.
  std::vector<float> top(F, -std::numeric_limits<float>::infinity());
  const auto activation = [&](std::size_t f, const Codes& x, std::size_t t) {
    float s = bias[f];
    for (std::size_t k = 0; k < W; k++) {
      const Code c = x[t + k];
      if (c < 4) s += Wt[(f * 4 + c) * W + k];
      else
        for (int b = 0; b < 4; b++) s += 0.25f * Wt[(f * 4 + static_cast<std::size_t>(b)) * W + k];
    }
    return s;
  };
  for (const auto& x : xs)
    for (std::size_t t = 0; t + W <= x.size(); t++)
      for (std::size_t f = 0; f < F; f++) top[f] = std::max(top[f], activation(f, x, t));
  std::vector<std::vector<double>> counts(F, std::vector<double>(4 * W, 0.0));
  for (const auto& x : xs)
    for (std::size_t t = 0; t + W <= x.size(); t++)
      for (std::size_t f = 0; f < F; f++) {
        if (!(top[f] > 0) || activation(f, x, t) <= 0.5f * top[f]) continue;
        out[f].sites++;
        for (std::size_t k = 0; k < W; k++)
          if (x[t + k] < 4) counts[f][4 * k + x[t + k]] += 1;
      }
  for (std::size_t f = 0; f < F; f++) {
    out[f].filter = f;
    out[f].pwm = annotation::pwm_from_counts(counts[f], 0.5);
    if (out[f].sites)
      for (double v : annotation::information_content(out[f].pwm)) out[f].information += v;
    out[f].weights.resize(4 * W);
    for (std::size_t k = 0; k < W; k++)
      for (std::size_t b = 0; b < 4; b++) out[f].weights[4 * k + b] = static_cast<double>(Wt[(f * 4 + b) * W + k]);
  }
  return out;
}

// ---- known motifs and matching ---------------------------------------------------------------------------------
struct KnownMotif {
  std::string id, name;
  std::size_t width = 0;
  std::vector<double> probability;  // width × 4, base fastest
};

// MEME's motif format (version 4): "MOTIF id name", "letter-probability matrix: … w= W …", W rows.
inline std::vector<KnownMotif> read_meme_motifs(std::string_view text, std::vector<sequence::Problem>* problems = nullptr) {
  std::vector<KnownMotif> out;
  KnownMotif* cur = nullptr;
  std::size_t rows = 0;
  sequence::detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    const auto t = sequence::detail::trim(line);
    if (t.substr(0, 6) == "MOTIF ") {
      out.emplace_back();
      cur = &out.back();
      const auto rest = sequence::detail::trim(t.substr(6));
      const auto sp = rest.find(' ');
      cur->id = std::string(rest.substr(0, sp));
      cur->name = sp == std::string_view::npos ? cur->id : std::string(sequence::detail::trim(rest.substr(sp + 1)));
      rows = 0;
      return;
    }
    if (!cur) return;
    if (t.substr(0, 25) == "letter-probability matrix") {
      const auto w = t.find("w=");
      if (w != std::string_view::npos) cur->width = static_cast<std::size_t>(std::strtoul(std::string(t.substr(w + 2)).c_str(), nullptr, 10));
      rows = cur->width;
      return;
    }
    if (rows > 0 && !t.empty()) {
      double v[4];
      const std::string s(t);
      char* p = const_cast<char*>(s.c_str());
      for (double& x : v) x = std::strtod(p, &p);
      double total = v[0] + v[1] + v[2] + v[3];
      if (!(total > 0.9 && total < 1.1)) {
        if (problems) problems->push_back({number, "a row of " + cur->id + " does not sum to one"});
        total = total > 0 ? total : 1;
      }
      for (double x : v) cur->probability.push_back(x / total);
      rows--;
    }
  });
  out.erase(std::remove_if(out.begin(), out.end(), [](const KnownMotif& m) { return m.width == 0 || m.probability.size() != 4 * m.width; }), out.end());
  return out;
}

struct MotifMatch {
  std::size_t motif = 0;           // index among the known motifs
  double similarity = -1;          // mean column Pearson correlation over the overlap
  long offset = 0;                 // the known motif's first column against the query's
  bool reverse = false;            // the known motif's reverse complement matched
  std::size_t overlap = 0;
  double p_value = 1;              // empirical, against queries of shuffled columns
};

namespace detail {

inline double column_pcc(const double* a, const double* b) {
  double ma = 0, mb = 0;
  for (int i = 0; i < 4; i++) { ma += a[i]; mb += b[i]; }
  ma /= 4;
  mb /= 4;
  double sab = 0, saa = 0, sbb = 0;
  for (int i = 0; i < 4; i++) {
    sab += (a[i] - ma) * (b[i] - mb);
    saa += (a[i] - ma) * (a[i] - ma);
    sbb += (b[i] - mb) * (b[i] - mb);
  }
  return saa > 0 && sbb > 0 ? sab / std::sqrt(saa * sbb) : 0;
}

inline std::vector<double> reverse_complement_matrix(const std::vector<double>& p, std::size_t w) {
  std::vector<double> r(4 * w);
  for (std::size_t i = 0; i < w; i++)
    for (int b = 0; b < 4; b++) r[4 * (w - 1 - i) + static_cast<std::size_t>(3 - b)] = p[4 * i + static_cast<std::size_t>(b)];
  return r;
}

// The best mean column PCC of query q (wq) against target t (wt) over offsets with overlap ≥ min_overlap.
inline MotifMatch best_alignment(const std::vector<double>& q, std::size_t wq, const std::vector<double>& t, std::size_t wt, std::size_t min_overlap) {
  MotifMatch best;
  for (long off = -static_cast<long>(wt) + 1; off < static_cast<long>(wq); off++) {
    // target column j aligns with query column j + off
    double s = 0;
    std::size_t n = 0;
    for (std::size_t j = 0; j < wt; j++) {
      const long i = static_cast<long>(j) + off;
      if (i < 0 || i >= static_cast<long>(wq)) continue;
      s += column_pcc(&q[4 * static_cast<std::size_t>(i)], &t[4 * j]);
      n++;
    }
    if (n < min_overlap) continue;
    const double v = s / static_cast<double>(n);
    if (v > best.similarity) { best.similarity = v; best.offset = off; best.overlap = n; }
  }
  return best;
}

}  // namespace detail

// The best match of a query (W × 4 probabilities) among known motifs, with an empirical p-value
// from `shuffles` queries whose columns are permuted (each against every motif, its best).
inline MotifMatch match_motif(const std::vector<double>& query, std::size_t width, const std::vector<KnownMotif>& known, std::size_t min_overlap = 5,
                              std::size_t shuffles = 100, std::uint64_t seed = 1) {
  const auto best_of = [&](const std::vector<double>& q) {
    MotifMatch best;
    for (std::size_t k = 0; k < known.size(); k++) {
      const auto& m = known[k];
      const std::size_t mo = std::min({min_overlap, width, m.width});
      for (int strand = 0; strand < 2; strand++) {
        const auto target = strand ? detail::reverse_complement_matrix(m.probability, m.width) : m.probability;
        auto a = detail::best_alignment(q, width, target, m.width, mo);
        if (a.similarity > best.similarity) {
          a.motif = k;
          a.reverse = strand == 1;
          best = a;
        }
      }
    }
    return best;
  };
  MotifMatch out = best_of(query);
  if (shuffles == 0) return out;
  calx::random::Stream r(seed, SHUFFLE_COLUMNS, 0, 0);
  std::size_t as_good = 0;
  std::vector<std::size_t> cols(width);
  for (std::size_t s = 0; s < shuffles; s++) {
    std::iota(cols.begin(), cols.end(), std::size_t{0});
    for (std::size_t i = width; i > 1; i--) std::swap(cols[i - 1], cols[r.below(i)]);
    std::vector<double> q(4 * width);
    for (std::size_t i = 0; i < width; i++)
      for (std::size_t b = 0; b < 4; b++) q[4 * i + b] = query[4 * cols[i] + b];
    as_good += best_of(q).similarity >= out.similarity;
  }
  out.p_value = (static_cast<double>(as_good) + 1) / (static_cast<double>(shuffles) + 1);
  return out;
}

// A filter's matrix trimmed to its informative core: the longest run of columns from the first to
// the last whose information exceeds `bits`.
inline std::pair<std::size_t, std::size_t> informative_core(const annotation::Pwm& p, double bits = 0.3) {
  const auto ic = annotation::information_content(p);
  std::size_t a = 0, b = ic.size();
  while (a < b && ic[a] < bits) a++;
  while (b > a && ic[b - 1] < bits) b--;
  return {a, b};
}

}  // namespace genomics::regulation
