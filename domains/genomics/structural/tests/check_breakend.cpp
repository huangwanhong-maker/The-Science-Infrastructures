// Checks of breakends and adjacencies: the six breakend records of the VCF 4.2 specification's
// example of every orientation (§5.4, Figure 1) and its example of inserted bases (Figure 2)
// written back exactly; the microhomology of a junction against the brute-force enumeration of
// every placement that leaves the derivative unchanged, for deletions, tandem duplications and
// translocations in random sequence with planted repeats; the canonical form is the leftmost such
// placement, the same for every member of the class; reversal is an involution and the signature
// names each rearrangement.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/vcf.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <set>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::structural;
using sequence::Codes;

namespace {

std::string letters(const Codes& x, std::size_t b, std::size_t e) {
  std::string s;
  for (std::size_t i = b; i < e; i++) s += sequence::dna_letter(x[i]);
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- VCF 4.2 §5.4: the breakends of Figure 1, and Figure 2's inserted bases --------------------
  {
    // Contigs "2", "13", "17" carrying the bases the records name at their positions (1-based).
    std::vector<Codes> contigs(3, Codes(330000, sequence::N));
    const std::vector<std::string> names{"2", "13", "17"};
    const auto put = [&](std::size_t c, std::size_t pos1, char b) { contigs[c][pos1 - 1] = sequence::dna_code(b); };
    put(0, 321681, 'G');
    put(0, 321682, 'T');
    put(1, 123456, 'C');
    put(1, 123457, 'A');
    put(2, 198982, 'A');
    put(2, 198983, 'C');
    const auto B = [](std::uint32_t c, std::size_t pos1, bool plus) { return Breakend{c, pos1 - 1, plus}; };
    // The three adjacencies of the 3-break: W–Y, V–U, X–Z.
    const Adjacency wy{B(0, 321681, true), B(2, 198982, true), {}};
    const Adjacency vu{B(0, 321682, false), B(1, 123456, true), {}};
    const Adjacency xz{B(1, 123457, false), B(2, 198983, false), {}};
    const auto alt = [&](const Adjacency& j) { return detail::breakend_alt(contigs, names, j); };
    s.require("bnd_W: G]17:198982]", alt(wy) == "G]17:198982]", alt(wy));
    s.require("bnd_Y: A]2:321681]", alt(reversed(wy)) == "A]2:321681]", alt(reversed(wy)));
    s.require("bnd_V: ]13:123456]T", alt(vu) == "]13:123456]T", alt(vu));
    s.require("bnd_U: C[2:321682[", alt(reversed(vu)) == "C[2:321682[", alt(reversed(vu)));
    s.require("bnd_X: [17:198983[A", alt(xz) == "[17:198983[A", alt(xz));
    s.require("bnd_Z: [13:123457[C", alt(reversed(xz)) == "[13:123457[C", alt(reversed(xz)));
    // Figure 2: the derivative reads 13:…C AGTNNNNNCA T…2:321682, leaving 13 at 123456 (+).
    const Adjacency ins{B(1, 123456, true), B(0, 321682, false), "AGTNNNNNCA"};
    s.require("inserted bases: bnd_U CAGTNNNNNCA[2:321682[", alt(ins) == "CAGTNNNNNCA[2:321682[", alt(ins));
    s.require("and its mate bnd_V ]13:123456]AGTNNNNNCAT", alt(reversed(ins)) == "]13:123456]AGTNNNNNCAT", alt(reversed(ins)));
    s.require("reversal is an involution, also with novel bases", reversed(reversed(ins)).inserted == ins.inserted && reversed(reversed(ins)).a == ins.a);
  }

  // ---- microhomology against brute force --------------------------------------------------------
  {
    calx::random::Stream r(5, 0x42524b31, 0, 0);
    std::size_t trials = 0, agree = 0, canonical_agree = 0, members_agree = 0;
    for (int t = 0; t < 600; t++) {
      // A random contig with a repeat planted at the two ends of the event, so that the junction slides.
      const std::size_t n = 400;
      Codes x(n);
      for (auto& b : x) b = static_cast<sequence::Code>(r.below(4));
      const std::size_t a = 100 + r.below(50), L = 60 + r.below(100), h = r.below(9);
      for (std::size_t k = 0; k < h && a + L + k < n; k++) x[a + L + k] = x[a + k];  // ref[a, a + h) = ref[a + L, a + L + h)
      const std::vector<Codes> contigs{x};
      const int kind = t % 2;  // 0 deletion of [a, a + L), 1 tandem duplication of [a, a + L)
      std::string derivative;
      if (kind == 0) derivative = letters(x, 0, a) + letters(x, a + L, n);
      else derivative = letters(x, 0, a + L) + letters(x, a, n);
      // Every placement of the same kind and length that makes the same word.
      std::set<std::size_t> same;
      for (std::size_t b = 1; b + L < n; b++) {
        const std::string w = kind == 0 ? letters(x, 0, b) + letters(x, b + L, n) : letters(x, 0, b + L) + letters(x, b, n);
        if (w == derivative) same.insert(b);
      }
      const Adjacency j = kind == 0 ? Adjacency{{0, a - 1, true}, {0, a + L, false}, {}} : Adjacency{{0, a + L - 1, true}, {0, a, false}, {}};
      const Homology hm = homology(contigs, j);
      trials++;
      agree += hm.length() + 1 == same.size();
      const Adjacency c = canonical(contigs, j);
      const std::size_t lowest = *same.begin();
      const Adjacency expected = kind == 0 ? ordered(Adjacency{{0, lowest - 1, true}, {0, lowest + L, false}, {}}) : ordered(Adjacency{{0, lowest + L - 1, true}, {0, lowest, false}, {}});
      canonical_agree += c.a == expected.a && c.b == expected.b;
      // Every member of the class has the same canonical form.
      bool all = true;
      for (const std::size_t b : same) {
        const Adjacency m = kind == 0 ? Adjacency{{0, b - 1, true}, {0, b + L, false}, {}} : Adjacency{{0, b + L - 1, true}, {0, b, false}, {}};
        const Adjacency cm = canonical(contigs, m);
        all = all && cm.a == c.a && cm.b == c.b;
      }
      members_agree += all;
    }
    s.require("a deletion's or tandem duplication's microhomology is the number of its equivalent placements less one (600 random junctions, brute force)", agree == trials, std::to_string(agree) + " of " + std::to_string(trials));
    s.require("its canonical form is the leftmost placement in the deletion's reading (VCF's left alignment), and the tandem duplication's in its own", canonical_agree == trials, std::to_string(canonical_agree));
    s.require("every equivalent placement has the same canonical form", members_agree == trials, std::to_string(members_agree));
  }
  {
    // A translocation's junction: der(A) = A[0, p) B[q, |B|); its placements by brute force.
    calx::random::Stream r(6, 0x42524b32, 0, 0);
    std::size_t trials = 0, agree = 0;
    for (int t = 0; t < 300; t++) {
      Codes A(300), Bc(300);
      for (auto& b : A) b = static_cast<sequence::Code>(r.below(4));
      for (auto& b : Bc) b = static_cast<sequence::Code>(r.below(4));
      const std::size_t p = 100 + r.below(100), q = 100 + r.below(100), h = r.below(7);
      for (std::size_t k = 0; k < h; k++) Bc[q - 1 - k] = A[p - 1 - k];  // the bases before both cuts agree: the junction slides left
      const std::vector<Codes> contigs{A, Bc};
      const std::string der = letters(A, 0, p) + letters(Bc, q, Bc.size());
      std::size_t same = 0;
      for (std::size_t m = 0; m + 1 < std::min(p, q); m++)  // the cut moved m to the left on both
        same += letters(A, 0, p - m) + letters(Bc, q - m, Bc.size()) == der;
      for (std::size_t m = 1; p + m < A.size() && q + m < Bc.size(); m++)  // or to the right
        same += letters(A, 0, p + m) + letters(Bc, q + m, Bc.size()) == der;
      const Homology hm = homology(contigs, {{0, p - 1, true}, {1, q, false}, {}});
      trials++;
      agree += hm.length() + 1 == same;
    }
    s.require("a translocation junction's microhomology is its equivalent placements less one (300, brute force)", agree == trials, std::to_string(agree) + " of " + std::to_string(trials));
  }

  // ---- signatures -------------------------------------------------------------------------------
  {
    const Adjacency del{{0, 99, true}, {0, 500, false}, {}}, dup{{0, 100, false}, {0, 499, true}, {}}, inv1{{0, 99, true}, {0, 499, true}, {}}, inv2{{0, 100, false}, {0, 500, false}, {}},
        tra{{0, 99, true}, {1, 7, false}, {}};
    s.require("signatures: (+,−) a deletion, (−,+) a tandem duplication, (+,+) and (−,−) an inversion's junctions, two contigs a translocation",
              signature(del) == Signature::deletion && signature(dup) == Signature::duplication && signature(inv1) == Signature::inversion_plus &&
                  signature(inv2) == Signature::inversion_minus && signature(tra) == Signature::translocation);
    s.require("ordering writes a ≤ b, and is idempotent", ordered(reversed(del)).a == del.a && ordered(ordered(dup)).b == dup.b);
  }
  return s.done();
}
