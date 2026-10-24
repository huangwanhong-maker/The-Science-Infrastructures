// FASTA and FASTQ, read from text and written to text, and the Phred scale of qualities.
//
// The readers take the whole text (the caller reads the file whose path it was given) and
// never throw on content: what cannot be read becomes a problem, with its line and why, and
// the rest is read on.
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

inline std::vector<Record> read_fasta(std::string_view text, std::vector<Problem>* problems = nullptr) {
  std::vector<Record> out;
  bool warned = false;
  detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (line.empty() || line.front() == ';') return;
    if (line.front() == '>') {
      out.emplace_back();
      detail::header(line.substr(1), out.back());
      if (out.back().id.empty() && problems) problems->push_back({number, "a record without an identifier"});
      return;
    }
    if (out.empty()) {
      if (!warned && problems) problems->push_back({number, "sequence before the first '>' header: skipped"});
      warned = true;
      return;
    }
    for (char c : line)
      if (c != ' ' && c != '\t') out.back().sequence += c;
  });
  return out;
}

inline std::vector<Record> read_fastq(std::string_view text, std::vector<Problem>* problems = nullptr) {
  std::vector<Record> out;
  enum class State { header, sequence, quality } state = State::header;
  Record r;
  std::size_t started = 0;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    switch (state) {
      case State::header:
        if (line.empty()) return;
        if (line.front() != '@') { problem(number, "expected a FASTQ header beginning with '@'"); return; }
        r = Record{};
        detail::header(line.substr(1), r);
        started = number;
        state = State::sequence;
        return;
      case State::sequence:
        if (!line.empty() && line.front() == '+') {
          const auto title = detail::trim(line.substr(1));
          if (!title.empty() && title.substr(0, r.id.size()) != r.id) problem(number, "the '+' line names another record than its header");
          state = State::quality;
          if (r.sequence.empty()) { out.push_back(r); state = State::header; }
          return;
        }
        for (char c : line)
          if (c != ' ' && c != '\t') r.sequence += c;
        return;
      case State::quality:
        r.quality += line;
        if (r.quality.size() >= r.sequence.size()) {
          if (r.quality.size() > r.sequence.size()) {
            problem(number, "record " + r.id + ": the quality is longer than the sequence, and is cut to it");
            r.quality.resize(r.sequence.size());
          }
          out.push_back(std::move(r));
          state = State::header;
        }
        return;
    }
  });
  if (state != State::header) problem(started, "record " + r.id + " ends before its quality is complete");
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
