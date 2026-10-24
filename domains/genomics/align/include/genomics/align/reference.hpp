// The reference every alignment of the domain is held to: morphe's generic sum over the paths of
// the edit graph (morphe/formal/alignment.hpp), with the domain's scoring as weights in the
// tropical semiring (max, +) on three states.
//
// A Scoring (types.hpp) is the affine scheme of morphe::formal::affine: the pair (a, b) weighs
// s(a, b), entering a gap weighs −open and each letter of a gap −extend, so a gap of length L
// weighs −(open + L · extend), and the modes are the same four. The generic programme is slow
// (sixteen semiring operations per vertex, the whole table kept for traceback) and obviously
// correct; the schedules of this module (gotoh.hpp, hirschberg.hpp, banded.hpp, and the kernels
// on lanes) are checked against it on random and homologous pairs.
//
// The same scheme in another semiring gives other readings of one alignment problem:
// reference_count counts the optimal alignments exactly (the tropical semiring with
// multiplicities), reference_alignments counts all of them (Counting).
//
// Complexity: O(mn) semiring operations; reference_alignment keeps O(mn) values.

#pragma once

#include <genomics/align/types.hpp>

#include <morphe/arithmetic/integer.hpp>
#include <morphe/formal/alignment.hpp>
#include <morphe/formal/semiring.hpp>

#include <cstdint>

namespace genomics::align {

inline morphe::formal::Mode formal_mode(Mode m) {
  switch (m) {
    case Mode::global: return morphe::formal::Mode::global;
    case Mode::local: return morphe::formal::Mode::local;
    case Mode::semiglobal: return morphe::formal::Mode::semiglobal;
    case Mode::overlap: return morphe::formal::Mode::overlap;
  }
  return morphe::formal::Mode::global;
}

// The scoring as morphe's affine scheme in (max, +). The scheme refers to `s`, which must outlive it.
inline auto formal_weights(const Scoring& s) {
  using T = morphe::formal::MaxPlus<std::int64_t>;
  const std::int64_t e = -static_cast<std::int64_t>(s.gap_extend);
  return morphe::formal::affine<T>([&s](Code a, Code b) { return static_cast<std::int64_t>(s.matrix(a, b)); }, -static_cast<std::int64_t>(s.gap_open),
                                   [e](Code) { return e; }, [e](Code) { return e; });
}

// The optimal score by the generic programme.
inline std::int64_t reference_score(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  return morphe::formal::path_sum(query, target, formal_weights(s), formal_mode(m));
}

// An optimal alignment by the generic programme's traceback.
inline Alignment reference_alignment(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  const auto best = morphe::formal::optimal_path(query, target, formal_weights(s), formal_mode(m));
  Alignment a;
  a.score = best.value;
  a.query_begin = best.path.x_begin;
  a.target_begin = best.path.y_begin;
  a.query_end = best.path.x_end();
  a.target_end = best.path.y_end();
  for (const auto& r : best.path.runs()) push(a.cigar, r.op, static_cast<std::uint32_t>(r.length));
  return a;
}

// The optimal score and the number of optimal alignments, exactly.
inline morphe::formal::ScoreCount reference_count(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  using T = morphe::formal::MaxPlusCount;
  const auto e = T::single(-static_cast<std::int64_t>(s.gap_extend));
  const auto w = morphe::formal::affine<T>([&s](Code a, Code b) { return T::single(s.matrix(a, b)); }, T::single(-static_cast<std::int64_t>(s.gap_open)),
                                           [e](Code) { return e; }, [e](Code) { return e; });
  return morphe::formal::path_sum(query, target, w, formal_mode(m));
}

// The number of alignments the mode allows (a function of the lengths alone; global: Delannoy).
inline morphe::arithmetic::Integer reference_alignments(const Codes& query, const Codes& target, Mode m) {
  using C = morphe::formal::Counting;
  const auto one = [](auto...) { return morphe::arithmetic::Integer(1); };
  return morphe::formal::path_sum(query, target, morphe::formal::linear<C>(one, one, one), formal_mode(m));
}

}  // namespace genomics::align
