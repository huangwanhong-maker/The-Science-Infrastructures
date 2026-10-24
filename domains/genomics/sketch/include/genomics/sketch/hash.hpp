// The hash of canonical k-mers that the sketches sample: a k-mer (k ≤ 32, kmer.hpp) and its reverse
// complement have one canonical integer, which is mixed into 64 bits by SplitMix64's finaliser
// (Steele, Lea and Flood 2014, "Fast splittable pseudorandom number generators", OOPSLA, 453–472,
// the function mix64; Stafford's variant 13 of MurmurHash3's finaliser), after the seed is
// folded in.
//
// Every step of the mix (xor with a right shift, multiplication by an odd constant) is invertible
// modulo 2^64, so for a fixed seed the hash is a bijection of 64-bit words: distinct k-mers never
// collide, and a sketch's error is the sampling error alone. The sketches assume the hashes of
// distinct k-mers behave as independent uniform draws; the mix is not a cryptographic hash, but
// it passes the statistical tests the estimators' checks make (Jaccard, containment and
// cardinality within their theoretical errors). hash64 of kmer.hpp is weaker at k ≤ 15 and
// covers only [0, 4^k), so it is not used here: the FracMinHash threshold and HyperLogLog's
// registers need all 64 bits uniform.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::sketch {

using sequence::Code;
using sequence::Codes;
using sequence::Kmer;

constexpr std::uint64_t mix64(std::uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// The hash of a canonical k-mer under a seed.
constexpr std::uint64_t kmer_hash(Kmer canonical, std::uint64_t seed) { return mix64(canonical ^ mix64(seed + 0x9e3779b97f4a7c15ULL)); }

// f(hash) for every canonical k-mer free of N of a sequence, in order (with repeats).
template <class F>
void for_each_kmer_hash(const Codes& s, int k, std::uint64_t seed, F&& f) {
  sequence::for_each_kmer(s, k, [&](std::size_t, Kmer forward, Kmer reverse) { f(kmer_hash(forward < reverse ? forward : reverse, seed)); });
}
template <class F>
void for_each_kmer_hash(const std::vector<Codes>& contigs, int k, std::uint64_t seed, F&& f) {
  for (const auto& c : contigs) for_each_kmer_hash(c, k, seed, f);
}

}  // namespace genomics::sketch
