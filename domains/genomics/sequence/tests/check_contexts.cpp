// Checks of the classes of substitutions (contexts.hpp): SBS-96 as the orbits of the strand group on
// the 192 directed trinucleotide substitutions, each orbit of two members and every class reached;
// names and the lexicographic order; DBS-78 as the 78 orbits of the 144 directed doublets, the list
// held to SigProfilerMatrixGenerator 1.4.0's (the fixture somatic-classes.txt, made from its source
// by make_somatic_fixtures.py).

#include <genomics/sequence/contexts.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::sequence;

int main() {
  calx::verification::Suite s;

  // ---- SBS-96 ---------------------------------------------------------------------------------------
  std::vector<int> preimages(SBS96, 0);
  bool strand = true, pyrimidine = true;
  for (Code l = 0; l < 4; l++)
    for (Code r = 0; r < 4; r++)
      for (Code a = 0; a < 4; a++)
        for (Code x = 0; x < 4; x++) {
          const int k = sbs96_index(l, r, a, x);
          if (a == r) { strand = strand && k == -1; continue; }
          if (k < 0) { strand = false; continue; }
          preimages[static_cast<std::size_t>(k)]++;
          // The same event read on the other strand.
          strand = strand && sbs96_index(complement(x), complement(r), complement(a), complement(l)) == k;
          const SubstitutionContext c = sbs96_context(static_cast<std::size_t>(k));
          pyrimidine = pyrimidine && (c.ref == C || c.ref == T);
          if (r == C || r == T) pyrimidine = pyrimidine && c.left == l && c.ref == r && c.alt == a && c.right == x;
        }
  s.require("SBS-96: a class for every directed substitution, the same on both strands", strand);
  s.require("SBS-96: every class has exactly two directed preimages (the strand group acts freely)", std::all_of(preimages.begin(), preimages.end(), [](int n) { return n == 2; }));
  s.require("SBS-96: classes are named on the pyrimidine strand", pyrimidine);
  s.require("SBS-96: the first, the 17th and the last classes are A[C>A]A, A[C>G]A, T[T>G]T",
            sbs96_name(0) == "A[C>A]A" && sbs96_name(16) == "A[C>G]A" && sbs96_name(95) == "T[T>G]T");
  {
    std::vector<std::string> names;
    for (std::size_t k = 0; k < SBS96; k++) names.push_back(sbs96_name(k));
    std::vector<std::string> sorted = names;
    std::sort(sorted.begin(), sorted.end());
    bool ok = true;
    std::set<std::size_t> seen;
    for (std::size_t i = 0; i < SBS96; i++) {
      ok = ok && names[sbs96_lexicographic(i)] == sorted[i];
      seen.insert(sbs96_lexicographic(i));
    }
    s.require("sbs96_lexicographic is the permutation sorting the names", ok && seen.size() == SBS96);
  }
  {
    const Codes x = encode_dna("TACGTA");
    s.require("C>T at TpCpG and G>A at CpGpT are one class, T[C>T]G and A[C>T]G apart",
              sbs96_at(x, 2, T) == sbs96_index(A, C, T, G) && sbs96_at(x, 3, A) == sbs96_index(A, C, T, G) && sbs96_name(static_cast<std::size_t>(sbs96_at(x, 2, T))) == "A[C>T]G");
    s.require("no class at either end or for N", sbs96_at(x, 0, C) == -1 && sbs96_at(x, 5, C) == -1 && sbs96_index(N, C, T, A) == -1);
  }

  // ---- DBS-78 -----------------------------------------------------------------------------------------
  std::vector<int> pre(DBS78, 0);
  bool both = true, closed = true;
  std::size_t directed = 0;
  for (Code r0 = 0; r0 < 4; r0++)
    for (Code r1 = 0; r1 < 4; r1++)
      for (Code a0 = 0; a0 < 4; a0++)
        for (Code a1 = 0; a1 < 4; a1++) {
          const int k = dbs78_index(r0, r1, a0, a1);
          if (r0 == a0 || r1 == a1) { both = both && k == -1; continue; }
          directed++;
          if (k < 0) { closed = false; continue; }
          pre[static_cast<std::size_t>(k)]++;
          closed = closed && dbs78_index(complement(r1), complement(r0), complement(a1), complement(a0)) == k;
        }
  s.require("DBS-78: 144 directed doublets, each in a class, the same on both strands; none for a single change", directed == 144 && closed && both);
  {
    std::size_t ones = 0, twos = 0;
    bool sizes = true;
    for (std::size_t k = 0; k < DBS78; k++) {
      const auto n = DBS78_NAMES[k];
      const bool palindrome = complement(dna_code(n[1])) == dna_code(n[0]);
      const bool alt_palindrome = complement(dna_code(n[4])) == dna_code(n[3]);
      // Orbits of size 1: a palindromic ref with a palindromic alternate; all others size 2.
      const int expect = palindrome && alt_palindrome ? 1 : 2;
      sizes = sizes && pre[k] == expect;
      (pre[k] == 1 ? ones : twos)++;
    }
    s.require("DBS-78: orbits of size 1 exactly where reference and alternate are palindromes (12 of them: 4 references, 3 alternates each), the rest of size 2 (66)", sizes && ones == 12 && twos == 66,
              std::to_string(ones) + " and " + std::to_string(twos));
  }
  {
    std::ifstream f(std::string(GENOMICS_FIXTURES) + "/somatic-classes.txt");
    std::vector<std::string> dbs;
    std::string line;
    while (std::getline(f, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.rfind("DBS78\t", 0) == 0) dbs.push_back(line.substr(6));
    }
    std::vector<std::string> mine;
    for (std::size_t k = 0; k < DBS78; k++) mine.push_back(dbs78_name(k));
    if (dbs.empty()) s.skip("DBS-78 against SigProfilerMatrixGenerator: the fixture somatic-classes.txt is missing");
    else s.require("DBS-78: the 78 names and their order are SigProfilerMatrixGenerator 1.4.0's", dbs == mine);
    std::vector<std::string> sbs;
    std::ifstream g(std::string(GENOMICS_FIXTURES) + "/somatic-classes.txt");
    while (std::getline(g, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.rfind("SBS96\t", 0) == 0) sbs.push_back(line.substr(6));
    }
    std::vector<std::string> lex;
    for (std::size_t i = 0; i < SBS96; i++) lex.push_back(sbs96_name(sbs96_lexicographic(i)));
    if (sbs.empty()) s.skip("SBS-96 against SigProfilerMatrixGenerator: the fixture is missing");
    else s.require("SBS-96: the 96 names in SigProfilerMatrixGenerator 1.4.0's matrix order are the lexicographic order here", sbs == lex);
  }
  return s.done();
}
