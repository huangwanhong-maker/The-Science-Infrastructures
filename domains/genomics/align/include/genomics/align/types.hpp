// The common language of every alignment in the domain: modes, scorings and penalties, CIGAR
// strings, and the alignment of a query to a target with its score and its path.
//
// - The query is the sequence being placed (a read), the target what it is placed on (a
//   reference); both are codes (sequence/alphabet.hpp).
// - Modes, by which ends of the edit graph's paths are free:
//     global      both sequences end to end (Needleman and Wunsch 1970);
//     local       any substring of each (Smith and Waterman 1981);
//     semiglobal  the query end to end, the target's ends free ("glocal": a read on a reference);
//     overlap     free end gaps on both: a suffix of one with a prefix of the other, or one
//                 contained in the other (the overlaps of assembly).
// - A Scoring is maximised: a substitution matrix s(a, b) on codes and affine gaps, a gap of
//   length L scoring −(open + L · extend), so that open = 0 is a linear gap (the convention of
//   BLAST and minimap2). Penalties are minimised, as the wavefront algorithm minimises them:
//   0 for a match, `mismatch` for a mismatch, open + L · extend for a gap.
// - A CIGAR is the path as runs of operations, in the SAM specification's letters (Li et al.
//   2009, Bioinformatics 25, 2078–2079, and The SAM/BAM Format Specification, §1.4): M an
//   aligned pair, = and X a match and a mismatch, I a base of the query against a gap
//   (an insertion to the target), D a base of the target against a gap (a deletion from it),
//   S soft-clipped query bases, H hard-clipped, N a skipped stretch of the target.
// - rescore() computes the score of any CIGAR under a Scoring, and penalty() under Penalties,
//   from the two sequences: every algorithm's claimed score is held to the score of the path it
//   returns, which is how the fast kernels are checked against one another.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::align {

using sequence::Code;
using sequence::Codes;

enum class Mode { global, local, semiglobal, overlap };

inline const char* mode_name(Mode m) {
  switch (m) {
    case Mode::global: return "global";
    case Mode::local: return "local";
    case Mode::semiglobal: return "semiglobal";
    case Mode::overlap: return "overlap";
  }
  return "?";
}

// A substitution matrix on the codes 0 .. size − 1.
struct Matrix {
  std::string name;
  int size = 0;
  std::vector<int> score;  // row-major, size × size

  int operator()(Code a, Code b) const { return score[static_cast<std::size_t>(a) * static_cast<std::size_t>(size) + b]; }
  int max_score() const { return score.empty() ? 0 : *std::max_element(score.begin(), score.end()); }
  int min_score() const { return score.empty() ? 0 : *std::min_element(score.begin(), score.end()); }

  // DNA: match on A, C, G, T; mismatch otherwise; N against anything scores `n`.
  static Matrix dna(int match = 2, int mismatch = -4, int n = -1) {
    Matrix m;
    m.name = "dna(" + std::to_string(match) + "," + std::to_string(mismatch) + ")";
    m.size = 5;
    m.score.assign(25, mismatch);
    for (int a = 0; a < 5; a++)
      for (int b = 0; b < 5; b++)
        m.score[static_cast<std::size_t>(a * 5 + b)] = (a == 4 || b == 4) ? n : a == b ? match : mismatch;
    return m;
  }
};

struct Scoring {
  Matrix matrix = Matrix::dna();
  int gap_open = 4, gap_extend = 2;  // a gap of length L scores −(gap_open + L · gap_extend)
  int gap(std::size_t length) const { return length ? -(gap_open + static_cast<int>(length) * gap_extend) : 0; }
};

struct Penalties {
  int mismatch = 4, gap_open = 6, gap_extend = 2;  // a match costs 0; a gap of length L, open + L · extend
};

// ---- CIGAR -------------------------------------------------------------------------------------
struct CigarOp {
  char op = 'M';
  std::uint32_t length = 0;
  friend bool operator==(const CigarOp&, const CigarOp&) = default;
};
using Cigar = std::vector<CigarOp>;

constexpr bool consumes_query(char op) { return op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X'; }
constexpr bool consumes_target(char op) { return op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X'; }

// Appends an operation, merging it with the last when they are the same.
inline void push(Cigar& c, char op, std::uint32_t length = 1) {
  if (length == 0) return;
  if (!c.empty() && c.back().op == op) c.back().length += length;
  else c.push_back({op, length});
}
inline std::string cigar_string(const Cigar& c) {
  std::string s;
  for (const auto& o : c) s += std::to_string(o.length) + o.op;
  return s.empty() ? "*" : s;
}
inline Cigar parse_cigar(std::string_view s) {
  Cigar c;
  if (s == "*") return c;
  std::uint32_t n = 0;
  bool digits = false;
  for (char ch : s) {
    if (ch >= '0' && ch <= '9') { n = n * 10 + static_cast<std::uint32_t>(ch - '0'); digits = true; continue; }
    if (!digits || std::string_view("MIDNSHP=X").find(ch) == std::string_view::npos) throw std::invalid_argument("genomics::align::parse_cigar: not a CIGAR: " + std::string(s));
    push(c, ch, n);
    n = 0;
    digits = false;
  }
  if (digits) throw std::invalid_argument("genomics::align::parse_cigar: a length without its operation: " + std::string(s));
  return c;
}
inline std::size_t query_length(const Cigar& c) {
  std::size_t n = 0;
  for (const auto& o : c) n += consumes_query(o.op) ? o.length : 0;
  return n;
}
inline std::size_t target_length(const Cigar& c) {
  std::size_t n = 0;
  for (const auto& o : c) n += consumes_target(o.op) ? o.length : 0;
  return n;
}
// The CIGAR reversed, as for the alignment of the reversed sequences.
inline Cigar reversed(Cigar c) {
  std::reverse(c.begin(), c.end());
  return c;
}

// ---- alignments --------------------------------------------------------------------------------------
struct Alignment {
  std::int64_t score = 0;
  std::size_t query_begin = 0, query_end = 0;    // the query's aligned part, [begin, end)
  std::size_t target_begin = 0, target_end = 0;  // the target's
  Cigar cigar;                                   // from (query_begin, target_begin), without clipping
};

// Two codes match when they are equal and below `limit`: 4 for DNA (N matches nothing), 20 for
// protein (B, Z, X and * match nothing).
constexpr bool same(Code a, Code b, Code limit = 4) { return a == b && a < limit; }

// The CIGAR with M split into = and X by the sequences.
inline Cigar explicit_matches(const Cigar& c, const Codes& query, const Codes& target, std::size_t qb = 0, std::size_t tb = 0, Code limit = 4) {
  Cigar out;
  std::size_t i = qb, j = tb;
  for (const auto& o : c) {
    if (o.op == 'M' || o.op == '=' || o.op == 'X') {
      for (std::uint32_t k = 0; k < o.length; k++, i++, j++) push(out, same(query[i], target[j], limit) ? '=' : 'X');
    } else {
      push(out, o.op, o.length);
      if (consumes_query(o.op)) i += o.length;
      if (consumes_target(o.op)) j += o.length;
    }
  }
  return out;
}

// The score of a path under a Scoring: the sum of its pairs and the costs of its gaps.
// Throws if the path leaves either sequence.
inline std::int64_t rescore(const Alignment& a, const Codes& query, const Codes& target, const Scoring& s) {
  std::int64_t score = 0;
  std::size_t i = a.query_begin, j = a.target_begin;
  for (const auto& o : a.cigar) {
    switch (o.op) {
      case 'M': case '=': case 'X':
        if (i + o.length > query.size() || j + o.length > target.size()) throw std::out_of_range("genomics::align::rescore: the path leaves a sequence");
        for (std::uint32_t k = 0; k < o.length; k++) score += s.matrix(query[i + k], target[j + k]);
        i += o.length;
        j += o.length;
        break;
      case 'I': score += s.gap(o.length); i += o.length; break;
      case 'D': score += s.gap(o.length); j += o.length; break;
      default: break;  // clipping and skips are not scored
    }
  }
  if (i != a.query_end || j != a.target_end) throw std::out_of_range("genomics::align::rescore: the path does not end where the alignment says");
  return score;
}

// The penalty of a path under Penalties (matches 0).
inline std::int64_t penalty(const Alignment& a, const Codes& query, const Codes& target, const Penalties& p, Code limit = 4) {
  std::int64_t cost = 0;
  std::size_t i = a.query_begin, j = a.target_begin;
  for (const auto& o : a.cigar) {
    switch (o.op) {
      case 'M': case '=': case 'X':
        for (std::uint32_t k = 0; k < o.length; k++) cost += same(query[i + k], target[j + k], limit) ? 0 : p.mismatch;
        i += o.length;
        j += o.length;
        break;
      case 'I': cost += p.gap_open + static_cast<std::int64_t>(o.length) * p.gap_extend; i += o.length; break;
      case 'D': cost += p.gap_open + static_cast<std::int64_t>(o.length) * p.gap_extend; j += o.length; break;
      default: break;
    }
  }
  return cost;
}

// The edit distance of a path: mismatches plus gap columns (SAM's NM).
inline std::size_t edit_distance(const Alignment& a, const Codes& query, const Codes& target, Code limit = 4) {
  std::size_t nm = 0, i = a.query_begin, j = a.target_begin;
  for (const auto& o : a.cigar) {
    if (o.op == 'M' || o.op == '=' || o.op == 'X') {
      for (std::uint32_t k = 0; k < o.length; k++) nm += !same(query[i + k], target[j + k], limit);
      i += o.length;
      j += o.length;
    } else if (o.op == 'I') { nm += o.length; i += o.length; }
    else if (o.op == 'D') { nm += o.length; j += o.length; }
  }
  return nm;
}

// Matches over aligned columns (gaps counted as columns): the BLAST identity.
inline double identity(const Alignment& a, const Codes& query, const Codes& target, Code limit = 4) {
  std::size_t columns = 0;
  for (const auto& o : a.cigar)
    if (o.op == 'M' || o.op == '=' || o.op == 'X' || o.op == 'I' || o.op == 'D') columns += o.length;
  return columns ? 1.0 - static_cast<double>(edit_distance(a, query, target, limit)) / static_cast<double>(columns) : 0;
}

// The two rows of the alignment and the line between them, for reading: '|' a match, '.' a
// mismatch, ' ' a gap; letters by `letter` (sequence::dna_letter by default).
struct Rendered {
  std::string query, middle, target;
};
inline Rendered render(const Alignment& a, const Codes& query, const Codes& target, char (*letter)(Code) = sequence::dna_letter) {
  Rendered r;
  std::size_t i = a.query_begin, j = a.target_begin;
  for (const auto& o : a.cigar)
    for (std::uint32_t k = 0; k < o.length; k++) {
      if (o.op == 'M' || o.op == '=' || o.op == 'X') {
        r.query += letter(query[i]);
        r.target += letter(target[j]);
        r.middle += query[i] == target[j] ? '|' : '.';
        i++;
        j++;
      } else if (o.op == 'I') {
        r.query += letter(query[i++]);
        r.target += '-';
        r.middle += ' ';
      } else if (o.op == 'D') {
        r.query += '-';
        r.target += letter(target[j++]);
        r.middle += ' ';
      }
    }
  return r;
}

}  // namespace genomics::align
