// The genomes a worker has loaded, and its threads: the state that lives between requests.
//
// A request names its genome in one of three ways, and the worker keeps the last few it was
// asked for, so that the many requests of one view (tracks at every zoom, a feature, a
// search) parse a file once:
//
//   {"file": "/path/NC_000913.3.gb"}              a GenBank or FASTA file the server points to
//   {"text": ">x\nACGT...", "name": "pasted"}     a sequence given inline
//   {"simulated": {"length": 200000, "gc": 0.5, "seed": 1, "repeats": {...}, "tandem": {...}}}
//
// A simulated genome is regenerated from its spec on every worker that is asked for it, so the
// cluster's workers on other machines need nothing but the spec.

#pragma once

#include "protocol.hpp"

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/cpu/target.hpp>
#include <calx/scheduler/parallel_for.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace workbench {

namespace seq = genomics::sequence;
namespace sim = genomics::simulate;

struct Contig {
  std::string name, accession, definition;
  seq::Codes codes;
  bool circular = false;
  std::vector<seq::Feature> features;  // GenBank's, source excluded
};

struct Genome {
  std::string key, title, organism, source;  // source: "file", "text" or "simulated"
  std::vector<std::string> lineage;
  std::vector<Contig> contigs;
  std::vector<sim::RepeatCopy> repeats;      // a simulated genome's truth
  int table = 1;                             // the genetic code its genes use
  std::vector<seq::Problem> problems;
  std::size_t length() const {
    std::size_t n = 0;
    for (const auto& c : contigs) n += c.codes.size();
    return n;
  }
  const Contig& contig(long long i) const {
    if (i < 0 || static_cast<std::size_t>(i) >= contigs.size()) throw Fault("no contig " + std::to_string(i) + " in " + title);
    return contigs[static_cast<std::size_t>(i)];
  }
};

namespace detail {

inline std::string read_text(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw Fault("cannot read " + path);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

// The genetic code most of the CDSs name, or the lineage's.
inline int code_of(const Genome& g) {
  std::map<int, int> votes;
  for (const auto& c : g.contigs)
    for (const auto& f : c.features)
      if (f.key == "CDS") {
        const auto t = f.qualifier("transl_table");
        votes[t.empty() ? 1 : std::atoi(t.c_str())]++;
      }
  if (!votes.empty()) return std::max_element(votes.begin(), votes.end(), [](auto& a, auto& b) { return a.second < b.second; })->first;
  for (const auto& l : g.lineage)
    if (l == "Bacteria" || l == "Archaea" || l == "Caudoviricetes" || l == "Duplodnaviria") return 11;
  return 1;
}

inline void parse_into(Genome& g, const std::string& text, const std::string& name) {
  std::size_t i = text.find_first_not_of(" \t\r\n");
  if (i == std::string::npos) throw Fault("the genome " + name + " is empty");
  if (text.compare(i, 5, "LOCUS") == 0) {
    for (auto& r : seq::read_genbank(text, &g.problems)) {
      Contig c;
      c.name = r.locus;
      c.accession = r.version.empty() ? r.accession : r.version;
      c.definition = r.definition;
      c.codes = seq::encode_dna(r.sequence);
      c.circular = r.circular();
      for (auto& f : r.features)
        if (f.key != "source") c.features.push_back(std::move(f));
      if (g.title.empty()) {
        g.title = r.definition;
        g.organism = r.organism;
        g.lineage = r.lineage;
      }
      g.contigs.push_back(std::move(c));
    }
  } else if (text[i] == '>') {
    for (auto& r : seq::read_fasta(text, &g.problems)) {
      Contig c;
      c.name = r.id;
      c.definition = r.description;
      c.codes = seq::encode_dna(r.sequence);
      if (g.title.empty()) g.title = r.description.empty() ? r.id : r.description;
      g.contigs.push_back(std::move(c));
    }
  } else {
    // Bare letters: one contig.
    Contig c;
    c.name = name;
    std::string letters;
    for (char ch : text)
      if (std::isalpha(static_cast<unsigned char>(ch))) letters += ch;
    c.codes = seq::encode_dna(letters);
    g.title = name;
    g.contigs.push_back(std::move(c));
  }
  if (g.contigs.empty()) throw Fault("no sequence could be read from " + name);
  g.table = code_of(g);
}

}  // namespace detail

struct Worker {
  calx::scheduler::ThreadPool pool{std::max(1u, std::thread::hardware_concurrency())};
  std::map<std::string, std::shared_ptr<Genome>> genomes;
  std::vector<std::string> recent;  // most recent last
  static constexpr std::size_t KEEP = 6;

  std::shared_ptr<Genome> genome(const Json& request) {
    const Json& spec = request["genome"];
    if (!spec.is_object()) throw Fault("the request names no genome");
    std::string key;
    if (spec.has("file")) key = "file:" + spec["file"].text();
    else if (spec.has("text")) key = "text:" + std::to_string(std::hash<std::string>{}(spec["text"].text())) + ":" + spec["name"].text("pasted");
    else if (spec.has("simulated")) key = "sim:" + spec["simulated"].dump();
    else throw Fault("a genome is a file, a text or a simulated spec");
    if (auto it = genomes.find(key); it != genomes.end()) {
      touch(key);
      return it->second;
    }
    auto g = std::make_shared<Genome>();
    g->key = key;
    if (spec.has("file")) {
      g->source = "file";
      detail::parse_into(*g, detail::read_text(spec["file"].text()), spec["file"].text());
    } else if (spec.has("text")) {
      g->source = "text";
      detail::parse_into(*g, spec["text"].text(), spec["name"].text("pasted"));
    } else {
      const Json& s = spec["simulated"];
      sim::GenomeSpec gs;
      gs.length = static_cast<std::size_t>(bounded(s, "length", 100000, 100, 20'000'000));
      gs.gc = num(s, "gc", 0.5);
      if (s["repeats"].is_object()) {
        const Json& r = s["repeats"];
        gs.repeats = {static_cast<std::size_t>(integer(r, "elements", 0)), static_cast<std::size_t>(integer(r, "copies", 0)),
                      static_cast<std::size_t>(integer(r, "length", 0)), num(r, "divergence", 0)};
      }
      if (s["tandem"].is_object()) {
        const Json& t = s["tandem"];
        gs.tandem = {static_cast<std::size_t>(integer(t, "arrays", 0)), static_cast<std::size_t>(integer(t, "unit", 0)),
                     static_cast<std::size_t>(integer(t, "copies", 0)), num(t, "divergence", 0)};
      }
      auto made = sim::make_genome(gs, static_cast<std::uint64_t>(integer(s, "seed", 1)));
      g->source = "simulated";
      g->title = s["name"].text("A simulated genome");
      Contig c;
      c.name = "sim";
      c.definition = g->title;
      c.codes = std::move(made.sequence);
      g->repeats = std::move(made.repeats);
      g->contigs.push_back(std::move(c));
      g->table = 11;
    }
    genomes[key] = g;
    touch(key);
    while (recent.size() > KEEP) {
      genomes.erase(recent.front());
      recent.erase(recent.begin());
    }
    return g;
  }

  // A structure derived from a genome (an index, a sketch), made once per worker and kept while
  // the genome is: derived<FmdIndex>(genome, "fmd", [&] { return FmdIndex(...); }).
  template <class T, class Make>
  std::shared_ptr<T> derived(const Genome& g, const std::string& kind, Make&& make, double* seconds = nullptr) {
    const std::string key = g.key + "\x1f" + kind;
    if (auto it = derived_.find(key); it != derived_.end()) {
      if (seconds) *seconds = built_in_[key];  // the time it took when it was made
      return std::static_pointer_cast<T>(it->second);
    }
    const auto t0 = std::chrono::steady_clock::now();
    auto made = std::make_shared<T>(make());
    if (seconds) *seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // Only the structures of genomes still held are kept.
    for (auto it = derived_.begin(); it != derived_.end();) {
      const auto sep = it->first.find('\x1f');
      it = genomes.count(it->first.substr(0, sep)) ? std::next(it) : derived_.erase(it);
    }
    derived_[key] = made;
    if (seconds) built_in_[key] = *seconds;
    return made;
  }

 private:
  std::map<std::string, std::shared_ptr<void>> derived_;
  std::map<std::string, double> built_in_;
  void touch(const std::string& key) {
    recent.erase(std::remove(recent.begin(), recent.end(), key), recent.end());
    recent.push_back(key);
  }
};

// Reads 0, …, count − 1 of a short-read sample, drawn in chunks on the worker's threads. Each
// read draws from its own stream (seed, index), so the reads are those of one call to
// simulate::short_reads, whatever the threads.
inline std::vector<sim::Read> short_reads(Worker& worker, const sim::Source& source, std::size_t count, const sim::ShortReadModel& model, std::uint64_t seed) {
  constexpr std::size_t CHUNK = 16384;
  const std::size_t chunks = (count + CHUNK - 1) / CHUNK;
  std::vector<std::vector<sim::Read>> parts(chunks);
  calx::scheduler::parallel_for(worker.pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t c = b; c < e; c++) parts[c] = sim::short_reads(source, std::min(CHUNK, count - c * CHUNK), model, seed, c * CHUNK);
  });
  std::vector<sim::Read> out;
  out.reserve(count);
  for (auto& part : parts) std::move(part.begin(), part.end(), std::back_inserter(out));
  return out;
}

}  // namespace workbench
