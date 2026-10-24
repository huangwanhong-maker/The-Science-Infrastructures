// The inputs of the toolbox's operations: sequences given inline (`text`: FASTA, FASTQ, GenBank,
// or bare letters as one record, DNA or protein) or as a file the server points to
// (`<name>_file`: fasta_file, reads_file, transcripts_file), and variants as a VCF file
// (`vcf_file`). A text's format is told from its first letters, never from its file's name:
// GenBank by LOCUS, FASTA by '>', FASTQ by '@'.
//
// Everything is read by the domain's readers (genomics::sequence's FASTA, FASTQ and GenBank, a
// FASTA or FASTQ file a line at a time, so that a file of gigabytes is never held), what they could
// not read is reported with its line, and every input is held within bounds, so that a pasted text
// or a large file cannot hang a worker. A tool says how much it reads (Bounds); where reading stops
// early, the result says so (Reading). A GenBank file is read whole, within the bound.

#pragma once

#include "protocol.hpp"
#include "genomes.hpp"

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genbank.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <fstream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

namespace workbench {

namespace detail {

// What an input held and how much of it was read.
struct Reading {
  std::string source = "text";  // "text", or the file's name
  std::string format;           // "fasta", "fastq", "genbank" or "letters"
  std::size_t records = 0, letters = 0, bytes = 0, file_bytes = 0;
  bool truncated = false;       // reading stopped before the input's end
  bool cut = false;             // ... at the bound of bytes, within a file
  std::string stopped;          // why
  std::vector<seq::Problem> problems;

  Json json() const {
    Json p = Json::array();
    for (std::size_t i = 0; i < problems.size() && i < 20; i++) p.push(Json::object({{"line", problems[i].line}, {"message", problems[i].message}}));
    Json out = Json::object({{"source", source}, {"format", format}, {"records", records}, {"letters", letters}, {"bytes", bytes}, {"truncated", truncated},
                             {"problems", p}, {"problem_count", problems.size()}});
    if (file_bytes) out.set("file_bytes", file_bytes);
    if (!stopped.empty()) out.set("stopped", stopped);
    return out;
  }
};

// How much of an input a tool reads: the bytes of a file, the records, the characters of a text
// given inline. `whole`: a tool whose answer would be wrong for part of its input refuses an input
// beyond the bounds rather than read part of it.
struct Bounds {
  std::size_t bytes = std::size_t{256} << 20;
  std::size_t records = 2'000'000;
  std::size_t text = std::size_t{16} << 20;
  bool whole = true;
};

inline double round_to(double x, double scale) { return std::round(x * scale) / scale; }

inline std::string megabytes(std::size_t b) {
  char buf[32];
  std::snprintf(buf, sizeof buf, b >= (std::size_t{10} << 20) ? "%.0f MB" : "%.1f MB", static_cast<double>(b) / 1e6);
  return buf;
}

inline std::string file_name(const std::string& path) {
  const auto at = path.find_last_of("/\\");
  return at == std::string::npos ? path : path.substr(at + 1);
}

// The format of a text from its first characters; VCF, SAM and Newick are named so that a tool
// given one where it wants sequences can say so.
inline std::string sniff(std::string_view head) {
  const std::size_t i = head.find_first_not_of(" \t\r\n");
  if (i == std::string_view::npos) return "empty";
  if (head.compare(i, 5, "LOCUS") == 0) return "genbank";
  if (head.compare(i, 16, "##fileformat=VCF") == 0) return "VCF";
  if (head.compare(i, 4, "@HD\t") == 0 || head.compare(i, 4, "@SQ\t") == 0) return "SAM";
  if (head[i] == '(') return "Newick";
  if (head[i] == '>') return "fasta";
  if (head[i] == '@') return "fastq";
  return "letters";
}

// Whether letters read as protein rather than DNA, judged on the first 5000: more than one in a
// hundred of them outside IUPAC's letters for bases (E, F, I, L, P, Q, …), or fewer than half of
// them A, C, G, T, U or N (a peptide spelled in the letters both alphabets share, MKV).
inline bool looks_protein(std::string_view s) {
  std::size_t n = 0, nucleotides = 0, other = 0;
  for (char c : s.substr(0, 5000)) {
    if (!std::isalpha(static_cast<unsigned char>(c))) continue;
    n++;
    nucleotides += std::string_view("ACGTUNacgtun").find(c) != std::string_view::npos;
    other += seq::iupac_mask(c) == 0;
  }
  return n > 0 && (100 * other > n || 2 * nucleotides < n);
}

// The message of a domain's exception without the name of the function that threw it
// ("genomics::index::iupac_words: …"), for the page.
inline std::string message_of(const std::exception& e) {
  const std::string_view w = e.what();
  const auto colon = w.find(": ");
  if (colon != std::string_view::npos && w.substr(0, colon).find("::") != std::string_view::npos && w.substr(0, colon).find(' ') == std::string_view::npos)
    return std::string(w.substr(colon + 2));
  return std::string(w);
}

// f(record) for a callback returning void or bool (false: no more records wanted).
template <class F, class R>
bool hand_on(F& f, R&& r) {
  if constexpr (std::is_same_v<std::invoke_result_t<F&, R&&>, bool>) return f(std::forward<R>(r));
  else {
    f(std::forward<R>(r));
    return true;
  }
}

// The lines of a text, f(line) for each, without their \n.
template <class F>
void text_lines(std::string_view text, F&& f) {
  std::size_t i = 0;
  while (i < text.size()) {
    std::size_t j = text.find('\n', i);
    if (j == std::string_view::npos) j = text.size();
    if (!f(text.substr(i, j - i))) return;
    i = j + 1;
  }
}

// The records of GenBank text, as records of their accession (with its version) and definition,
// their letters in capitals.
inline std::vector<seq::Record> genbank_records(std::string_view text, std::vector<seq::Problem>* problems) {
  std::vector<seq::Record> out;
  for (auto& g : seq::read_genbank(text, problems)) {
    seq::Record r;
    r.id = !g.version.empty() ? g.version : !g.accession.empty() ? g.accession : g.locus;
    r.description = g.definition;
    r.sequence = std::move(g.sequence);
    for (char& c : r.sequence) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    out.push_back(std::move(r));
  }
  return out;
}

}  // namespace detail

// Every record of a request's sequences, handed to f(seq::Record&&) as it is read: the file of
// `<name>_file` when the request names one, else its `text`. f may return false to stop.
template <class F>
detail::Reading each_record(const Json& request, const std::string& name, const detail::Bounds& b, F&& f) {
  detail::Reading reading;
  bool more = true;
  const auto give = [&](seq::Record&& r) {
    if (!more) return false;
    if (reading.records >= b.records) {
      if (b.whole) throw Fault("the input has more than " + std::to_string(b.records) + " records, more than this tool takes");
      reading.truncated = true;
      reading.stopped = "the first " + std::to_string(b.records) + " records were read";
      return more = false;
    }
    reading.records++;
    reading.letters += r.sequence.size();
    if (!detail::hand_on(f, std::move(r))) {
      reading.truncated = true;
      if (reading.stopped.empty()) reading.stopped = "the tool read the records it needed";
      more = false;
    }
    return more;
  };
  const std::string key = name + "_file";
  const std::string path = str(request, key.c_str());
  if (!path.empty()) {
    std::vector<char> buffer(std::size_t{1} << 20);
    std::ifstream in;
    in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    in.open(path, std::ios::binary);
    if (!in) throw Fault("cannot read " + detail::file_name(path));
    reading.source = detail::file_name(path);
    std::error_code ec;
    reading.file_bytes = static_cast<std::size_t>(std::filesystem::file_size(path, ec));
    if (b.whole && reading.file_bytes > b.bytes)
      throw Fault(reading.source + " has " + detail::megabytes(reading.file_bytes) + "; this tool reads files of up to " + detail::megabytes(b.bytes));
    std::string head(65536, '\0');
    in.read(head.data(), static_cast<std::streamsize>(head.size()));
    head.resize(static_cast<std::size_t>(in.gcount()));
    in.clear();
    in.seekg(0);
    reading.format = detail::sniff(head);
    if (reading.format == "empty") throw Fault(reading.source + " is empty");
    if (reading.format == "VCF" || reading.format == "SAM" || reading.format == "Newick") throw Fault(reading.source + " is " + reading.format + ", not sequences");
    if (reading.format == "genbank") {
      if (reading.file_bytes > b.bytes) throw Fault(reading.source + " is a GenBank file of " + detail::megabytes(reading.file_bytes) + "; this tool reads GenBank whole, up to " + detail::megabytes(b.bytes));
      const std::string text = detail::read_text(path);
      reading.bytes = text.size();
      for (auto& r : detail::genbank_records(text, &reading.problems))
        if (!give(std::move(r))) break;
      return reading;
    }
    // FASTA, FASTQ or bare letters, a line at a time, until the bound.
    seq::FastaReader fasta(&reading.problems);
    seq::FastqReader fastq(&reading.problems);
    seq::Record bare, done;
    bare.id = reading.source;
    std::string line;
    bool cut = false;
    while (more && std::getline(in, line)) {
      reading.bytes += line.size() + 1;
      if (reading.format == "fasta") {
        if (fasta.line(line, done) && !give(std::move(done))) break;
      } else if (reading.format == "fastq") {
        if (fastq.line(line, done) && !give(std::move(done))) break;
      } else {
        for (char c : line)
          if (std::isalpha(static_cast<unsigned char>(c)) || c == '*') bare.sequence += c;
      }
      if (reading.bytes > b.bytes && in.peek() != std::char_traits<char>::eof()) {
        if (b.whole) throw Fault(reading.source + " has " + detail::megabytes(reading.file_bytes) + "; this tool reads files of up to " + detail::megabytes(b.bytes));
        cut = true;
        break;
      }
    }
    if (more) {
      if (reading.format == "fasta" && fasta.finish(done)) give(std::move(done));  // when cut, a prefix of the last record
      else if (reading.format == "fastq" && !cut) fastq.finish();
      else if (reading.format == "letters") give(std::move(bare));
    }
    if (cut) {
      reading.truncated = reading.cut = true;
      reading.stopped = "the first " + detail::megabytes(reading.bytes) + " of " + detail::megabytes(reading.file_bytes) + " were read";
    }
    return reading;
  }
  const Json& t = request["text"];
  if (!t.is_string() || t.as_string().find_first_not_of(" \t\r\n") == std::string::npos) throw Fault("no sequences were given");
  const std::string& text = t.as_string();
  if (text.size() > b.text) throw Fault("a pasted text has at most " + detail::megabytes(b.text) + " here: give a file instead");
  reading.bytes = text.size();
  reading.format = detail::sniff(text);
  if (reading.format == "VCF" || reading.format == "SAM" || reading.format == "Newick") throw Fault("the text is " + reading.format + ", not sequences");
  if (reading.format == "genbank") {
    for (auto& r : detail::genbank_records(text, &reading.problems))
      if (!give(std::move(r))) break;
  } else if (reading.format == "fasta" || reading.format == "fastq") {
    seq::FastaReader fasta(&reading.problems);
    seq::FastqReader fastq(&reading.problems);
    seq::Record done;
    const bool fa = reading.format == "fasta";
    detail::text_lines(text, [&](std::string_view line) { return (fa ? fasta.line(line, done) : fastq.line(line, done)) ? give(std::move(done)) : true; });
    if (more) {
      if (fa && fasta.finish(done)) give(std::move(done));
      if (!fa) fastq.finish();
    }
  } else {
    seq::Record r;
    r.id = "sequence";
    for (char c : text)
      if (std::isalpha(static_cast<unsigned char>(c)) || c == '*') r.sequence += c;
    give(std::move(r));
  }
  return reading;
}

// Every record, held.
inline std::vector<seq::Record> all_records(const Json& request, const std::string& name, const detail::Bounds& b, detail::Reading& reading) {
  std::vector<seq::Record> out;
  reading = each_record(request, name, b, [&](seq::Record&& r) { out.push_back(std::move(r)); });
  if (out.empty()) throw Fault("no sequence could be read from " + reading.source);
  return out;
}

// The first record, and whether the input holds more (reading stops at the second).
struct FirstRecord {
  seq::Record record;
  bool more = false;
};
inline FirstRecord first_record(const Json& request, const std::string& name, const detail::Bounds& b, detail::Reading& reading) {
  FirstRecord out;
  std::size_t seen = 0;
  reading = each_record(request, name, b, [&](seq::Record&& r) {
    if (seen++ == 0) {
      out.record = std::move(r);
      return true;
    }
    out.more = true;
    return false;
  });
  if (seen == 0) throw Fault("no sequence could be read from " + reading.source);
  if (reading.cut && !out.more) throw Fault("the first record of " + reading.source + " is longer than this tool reads (" + detail::megabytes(b.bytes) + ")");
  if (out.more) {
    reading.records = 1;
    reading.letters = out.record.sequence.size();
    reading.stopped = "the first record, " + out.record.id + ", of several";
  }
  return out;
}

// A read as the mapper, the caller and the counters take it: its name, bases and qualities.
struct ReadRecord {
  std::string name;
  seq::Codes bases;
  std::string quality;  // Phred + 33, or empty for reads from FASTA
};

// Every read of the request's `reads_file` (or text), handed to f(ReadRecord&&).
template <class F>
detail::Reading each_read(const Json& request, const detail::Bounds& b, F&& f) {
  return each_record(request, "reads", b, [&](seq::Record&& r) {
    ReadRecord x{std::move(r.id), seq::encode_dna(r.sequence), std::move(r.quality)};
    return detail::hand_on(f, std::move(x));
  });
}

// Names made distinct: a second "x" becomes "x (2)".
inline std::vector<std::string> distinct_names(const std::vector<seq::Record>& records) {
  std::vector<std::string> out;
  std::set<std::string> seen;
  for (const auto& r : records) {
    std::string name = r.id.empty() ? "sequence" : r.id;
    for (int k = 2; seen.count(name); k++) name = (r.id.empty() ? "sequence" : r.id) + " (" + std::to_string(k) + ")";
    seen.insert(name);
    out.push_back(std::move(name));
  }
  return out;
}

// The genetic code a request names by NCBI's number (its `table`).
inline const seq::GeneticCode& genetic_code_of(const Json& request, int fallback) {
  const auto id = integer(request, "table", fallback);
  try {
    return seq::genetic_code(static_cast<int>(id));
  } catch (const std::invalid_argument&) {
    throw Fault("NCBI has no genetic code " + std::to_string(id));
  }
}

// A fault unless the letters read as DNA.
inline void require_dna(const seq::Record& r, const char* what) {
  if (detail::looks_protein(r.sequence)) throw Fault(r.id + " reads as protein; " + what + " is of DNA");
}

// The text of the request's VCF file: at most `bytes` of it, cut after its last whole line.
inline std::string vcf_text(const Json& request, std::size_t bytes, detail::Reading& reading) {
  const std::string path = str(request, "vcf_file");
  if (path.empty()) throw Fault("choose a VCF file");
  std::ifstream in(path, std::ios::binary);
  if (!in) throw Fault("cannot read " + detail::file_name(path));
  reading.source = detail::file_name(path);
  reading.format = "vcf";
  std::error_code ec;
  reading.file_bytes = static_cast<std::size_t>(std::filesystem::file_size(path, ec));
  std::string text(std::min(bytes, reading.file_bytes ? reading.file_bytes : bytes), '\0');
  in.read(text.data(), static_cast<std::streamsize>(text.size()));
  text.resize(static_cast<std::size_t>(in.gcount()));
  if (in.peek() != std::char_traits<char>::eof()) {
    const auto last = text.rfind('\n');
    text.resize(last == std::string::npos ? 0 : last + 1);
    reading.truncated = true;
    reading.stopped = "the first " + detail::megabytes(text.size()) + " of " + detail::megabytes(reading.file_bytes) + " were read";
  }
  reading.bytes = text.size();
  return text;
}

}  // namespace workbench
