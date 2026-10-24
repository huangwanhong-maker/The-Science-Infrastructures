// FASTA and FASTQ, read from text and written to text, and the Phred scale of qualities.
//
// The readers take the whole text (the caller reads the file whose path it was given) and
// never throw on content: what cannot be read becomes a problem, with its line and why, and
// the rest is read on. For files too large to hold, FastaReader and FastqReader take the same
// grammars a line at a time and hand on each record as it completes; read_fasta and read_fastq
// are those readers run over a whole text, so the two ways of reading never differ.
//
// - FASTA (Pearson and Lipman 1988): records begin with '>'; the identifier runs to the first
//   white space, the description is the rest of the line; sequence lines are joined with
//   white space removed; lines beginning with ';' are comments; line endings may be \n or
//   \r\n. Letters are kept as they are (codes come from alphabet.hpp).
// - FASTQ (Cock et al. 2010, Nucleic Acids Research 38, 1767–1771): '@' header, sequence lines
//   up to a '+' line (whose optional title must repeat the header's), then quality lines until
//   the quality is as long as the sequence. Qualities are Sanger's, Phred + 33.
// - Phred: Q = −10 log10 p (Ewing and Green 1998, Genome Research 8, 186–194).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::sequence {

// A fault found while reading, with where it is (a line, counted from 1) and why.
struct Problem {
  std::size_t line = 0;
  std::string message;
};

struct Record {
  std::string id, description;
  std::string sequence;  // letters, as in the file
  std::string quality;   // FASTQ only: one character per base, Phred + 33
};

namespace detail {

// Splits text into lines without their endings, calling f(line, number).
template <class F>
void for_each_line(std::string_view text, F&& f) {
  std::size_t number = 0, i = 0;
  while (i < text.size()) {
    std::size_t j = text.find('\n', i);
    if (j == std::string_view::npos) j = text.size();
    std::string_view line = text.substr(i, j - i);
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    f(line, ++number);
    i = j + 1;
  }
}

inline std::string_view trim(std::string_view s) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.remove_prefix(1);
  while (!s.empty() && (s.back() == ' ' || s.back() == '\t')) s.remove_suffix(1);
  return s;
}

inline void header(std::string_view line, Record& r) {
  line = trim(line);
  const std::size_t sp = line.find_first_of(" \t");
  r.id = std::string(line.substr(0, sp));
  r.description = sp == std::string_view::npos ? "" : std::string(trim(line.substr(sp)));
}

}  // namespace detail

// ---- a line at a time ------------------------------------------------------------------------
// Each reader takes the lines of its input in order, without their endings (a final '\r' is
// dropped, so lines of a \r\n file may be passed as std::getline gives them), and counts them
// from 1 for its problems. line() is true when the line completes a record, which it then moves
// into `out`; finish() ends the input. A FASTA record completes when the next header arrives, so
// the last one comes from finish(); a FASTQ record completes with its quality, and finish() notes
// one left incomplete.

class FastaReader {
 public:
  explicit FastaReader(std::vector<Problem>* problems = nullptr) : problems_(problems) {}

  bool line(std::string_view line, Record& out) {
    number_++;
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    if (line.empty() || line.front() == ';') return false;
    if (line.front() == '>') {
      const bool done = open_;
      if (done) out = std::move(current_);
      current_ = Record{};
      detail::header(line.substr(1), current_);
      open_ = true;
      if (current_.id.empty() && problems_) problems_->push_back({number_, "a record without an identifier"});
      return done;
    }
    if (!open_) {
      if (!warned_ && problems_) problems_->push_back({number_, "sequence before the first '>' header: skipped"});
      warned_ = true;
      return false;
    }
    for (char c : line)
      if (c != ' ' && c != '\t') current_.sequence += c;
    return false;
  }
  // The last record, when there is one.
  bool finish(Record& out) {
    if (!open_) return false;
    out = std::move(current_);
    open_ = false;
    return true;
  }
  std::size_t lines() const { return number_; }

 private:
  std::vector<Problem>* problems_;
  Record current_;
  bool open_ = false, warned_ = false;
  std::size_t number_ = 0;
};

class FastqReader {
 public:
  explicit FastqReader(std::vector<Problem>* problems = nullptr) : problems_(problems) {}

  bool line(std::string_view line, Record& out) {
    number_++;
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    switch (state_) {
      case State::header:
        if (line.empty()) return false;
        if (line.front() != '@') { problem(number_, "expected a FASTQ header beginning with '@'"); return false; }
        r_ = Record{};
        detail::header(line.substr(1), r_);
        started_ = number_;
        state_ = State::sequence;
        return false;
      case State::sequence:
        if (!line.empty() && line.front() == '+') {
          const auto title = detail::trim(line.substr(1));
          if (!title.empty() && title.substr(0, r_.id.size()) != r_.id) problem(number_, "the '+' line names another record than its header");
          if (r_.sequence.empty()) {  // no bases: no quality lines follow
            state_ = State::header;
            out = std::move(r_);
            return true;
          }
          state_ = State::quality;
          return false;
        }
        for (char c : line)
          if (c != ' ' && c != '\t') r_.sequence += c;
        return false;
      case State::quality:
        r_.quality += line;
        if (r_.quality.size() >= r_.sequence.size()) {
          if (r_.quality.size() > r_.sequence.size()) {
            problem(number_, "record " + r_.id + ": the quality is longer than the sequence, and is cut to it");
            r_.quality.resize(r_.sequence.size());
          }
          state_ = State::header;
          out = std::move(r_);
          return true;
        }
        return false;
    }
    return false;
  }
  // The end of the input: a record begun and not completed is a problem on its header's line.
  void finish() {
    if (state_ != State::header) problem(started_, "record " + r_.id + " ends before its quality is complete");
    state_ = State::header;
  }
  // Whether a record is begun and not yet complete.
  bool within_record() const { return state_ != State::header; }
  std::size_t lines() const { return number_; }

 private:
  enum class State { header, sequence, quality };
  void problem(std::size_t line, std::string m) {
    if (problems_) problems_->push_back({line, std::move(m)});
  }
  std::vector<Problem>* problems_;
  State state_ = State::header;
  Record r_;
  std::size_t started_ = 0, number_ = 0;
};

inline std::vector<Record> read_fasta(std::string_view text, std::vector<Problem>* problems = nullptr) {
  std::vector<Record> out;
  FastaReader reader(problems);
  Record r;
  detail::for_each_line(text, [&](std::string_view line, std::size_t) {
    if (reader.line(line, r)) out.push_back(std::move(r));
  });
  if (reader.finish(r)) out.push_back(std::move(r));
  return out;
}

inline std::vector<Record> read_fastq(std::string_view text, std::vector<Problem>* problems = nullptr) {
  std::vector<Record> out;
  FastqReader reader(problems);
  Record r;
  detail::for_each_line(text, [&](std::string_view line, std::size_t) {
    if (reader.line(line, r)) out.push_back(std::move(r));
  });
  reader.finish();
  return out;
}

inline std::string write_fasta(const std::vector<Record>& records, std::size_t width = 70) {
  std::string out;
  for (const auto& r : records) {
    out += '>' + r.id + (r.description.empty() ? "" : " " + r.description) + '\n';
    for (std::size_t i = 0; i < r.sequence.size(); i += width) out += r.sequence.substr(i, width) + '\n';
  }
  return out;
}

inline std::string write_fastq(const std::vector<Record>& records) {
  std::string out;
  for (const auto& r : records) {
    out += '@' + r.id + (r.description.empty() ? "" : " " + r.description) + '\n' + r.sequence + "\n+\n";
    out += r.quality.size() == r.sequence.size() ? r.quality : std::string(r.sequence.size(), 'I');
    out += '\n';
  }
  return out;
}

// ---- the Phred scale ------------------------------------------------------------------------
constexpr int phred_of(char q) { return static_cast<int>(static_cast<unsigned char>(q)) - 33; }
constexpr char phred_char(int q) { return static_cast<char>(std::clamp(q, 0, 93) + 33); }
inline double error_probability(int q) { return std::pow(10.0, -q / 10.0); }
// The Phred score of an error probability, rounded, at most 93.
inline int phred_from_probability(double p) {
  if (!(p > 0)) return 93;
  return std::clamp(static_cast<int>(std::lround(-10 * std::log10(p))), 0, 93);
}

}  // namespace genomics::sequence
