// The grammar of the HGVS Nomenclature (version 21.1, https://hgvs-nomenclature.org; den Dunnen et
// al. 2016, Human Mutation 37, 564–569; Hart et al. 2024, Genome Medicine 16, 149), written as BNF
// over bytes and compiled by morphe's BNF reader (morphe/formal/bnf.hpp, morphe decision 0009).
//
// The standard states its syntax in an "EBNF-like style" (recommendations/grammar.md) that it does
// not yet make exact: syntax.yaml's frameshift form, `aa_position "fs" ["Ter" position]`, leaves out
// the new residue its own example has (p.Arg97ProfsTer23). This grammar is written from the
// recommendations' definitions and every one of their examples (the checks parse them all,
// tests/fixtures/hgvs-examples.tsv), and is unambiguous on each (the checks count the derivations).
// It is the language L(G) of the theory note (theory/nomenclature.md §2):
//
// - A description is a reference, a coordinate type and a list of variants (general.md): `g` linear
//   genomic, `m` mitochondrial, `o` circular, `c` coding DNA, `n` non-coding DNA, `r` RNA, `p`
//   protein. The positions each type admits differ (numbering.md): g., m. and o. are plain numbers
//   (and pter, qter, cen of the HGVS/ISCN extension, DNA/complex.md); c. and r. take the 5′ UTR's
//   '-', the 3′ UTR's '*' and intronic offsets (c.88+1, c.89-2); n. takes offsets only. DNA is
//   written in capitals and RNA in lower case, both in IUPAC's letters (background/standards.md).
//   The grammar therefore has one family of rules for each nucleotide coordinate type, made from one
//   template (`family`), and one for proteins.
// - Variants are joined into alleles: in cis `[a;b]`, in trans `[a];[b]`, of unknown phase
//   `a(;)b`, transcripts of one allele `[a,b]` (RNA and protein), alternatives `a^b`; a variant or a
//   list in parentheses is predicted (`p.(Arg97ProfsTer23)`, `r.[(578c>u;1339a>g)]`); `::` joins
//   the junctions of a ring (complex.md) and the two partners of an adjoined transcript
//   (RNA/adjoined_transcript.md); `=/` and `=//` are mosaic and chimeric (DNA/other.md).
// - Edits (DNA/*.md): substitution `>` of one base, `del`, `dup`, `ins`, `delins`, `inv`, `=`,
//   repeated sequences `CAG[23]` and copy numbers `[3]`, `sup`, methylation `|gom`; inserted
//   material as a sequence, a count of unknown bases `N[15]`, a range of the same reference
//   (inverted or not), or a range of another `[NC_000022.10:g.35788169_35788352]`, several in brackets.
//   Proteins (protein/*.md): substitutions in three- or one-letter code, `Ter` and `*`, `del`, `dup`,
//   `ins`, `delins`, `fs` with the new stop, `ext` at either end, `Met1?`, repeats.
// - Uncertainty (uncertain.md): `(a_b)` for a position somewhere between a and b, `?` for an unknown
//   one, `(a_b)_(c_d)` for a range whose ends are uncertain.
//
// The grammar is deliberately a little wider than the standard where a check of meaning gives a
// better message than a syntax error: a substitution's location may be a range, and a range's two
// ends may be equal or reversed; hgvs.hpp's check reports these. Rules beginning with '_' are hidden
// (spliced into their parents); every other rule is a node of the named tree, which hgvs.hpp reads.

#pragma once

#include <morphe/formal/bnf.hpp>

#include <string>
#include <string_view>

namespace genomics::nomenclature {

namespace detail {

// The shared rules: the start, references, numbers and counts, separators, the protein family.
inline constexpr std::string_view grammar_head = R"BNF(
# ---- descriptions -------------------------------------------------------------------------------
hgvs            ::= description ( "^" description )* | genome_alleles ;
genome_alleles  ::= genome_allele ( ";" genome_allele )* ;
genome_allele   ::= "[" description ( ";" description )* "]" ;
description     ::= ( reference ":" )? _variants | adjoined ;
_variants       ::= g_variants | c_variants | n_variants | r_variants | p_variants ;

# A reference: an accession with its version, and a selector, as NG_012232.1(NM_004006.2).
reference       ::= accession selector? ;
accession       ::= identifier ( "." version )? ;
identifier      ::= [A-Za-z] [A-Za-z0-9_\-]+ ;
version         ::= [0-9]+ ;
selector        ::= "(" accession ")" ;

# Adjoined transcripts (RNA/adjoined_transcript.md): NM_002354.2:r.-358_555::NM_000251.2:r.212_*279.
adjoined        ::= adjoined_part "::" ( linker "::" )? adjoined_part ;
adjoined_part   ::= reference ":" "r." r_upos "_" r_upos ;
linker          ::= [acgubdhkmnrsvwy]+ ;

# ---- shared pieces --------------------------------------------------------------------------------
number          ::= [0-9]+ ;
offset          ::= [+\-] offset_amount ;
offset_amount   ::= number | "?" ;
utr             ::= [\-*] ;
count           ::= "[" copies "]" ;
copies          ::= number | "?" | "(" number "_" copies_upper ")" ;
copies_upper    ::= number | "?" ;
other_allele    ::= ";" count ;
separator       ::= ";" | "(;)" | "," ;
mixture         ::= "/" | "//" ;
inverted        ::= "inv" ;
supernumerary   ::= "sup" ;
methylation     ::= "|" methylation_state ;
methylation_state ::= "gom" | "lom" | "met=" ;
whole           ::= "=" | "?" | "0" | "0?" | "(=)" | "(?)" ;
r_whole         ::= whole | "spl" | "spl?" | "(spl)" | "(spl?)" ;
remote          ::= reference ":" _remote_location inverted? ;
_remote_location ::= g_type "." g_loc | "c." c_loc | "n." n_loc | "r." r_loc ;

# ---- positions of each coordinate type (numbering.md) --------------------------------------------
g_type          ::= [gmo] ;
g_pos           ::= number | "?" | "pter" | "qter" | "cen" ;
c_pos           ::= utr? number offset? | "?" ;
n_pos           ::= number offset? | "?" ;
r_pos           ::= utr? number offset? | "?" ;

# ---- proteins (protein/*.md) -----------------------------------------------------------------------
p_variants      ::= "p" "." p_list ;
p_list          ::= p_alternative ( "^" p_alternative )* | whole ;
p_alternative   ::= p_group ( "(;)" p_group )* ;
p_group         ::= p_allele ( ";" p_allele )* | p_chain ;
p_allele        ::= "[" p_body "]" ;
p_body          ::= p_predicted_list | p_chain ( separator p_chain )* | whole ;
p_predicted_list ::= "(" p_chain separator p_chain ( separator p_chain )* ")" ;
p_chain         ::= p_item ( "::" p_item )* ;
p_item          ::= p_edit | p_predicted | p_mixed ;
p_predicted     ::= "(" p_edit ( "^" p_edit )* ")" ;
p_mixed         ::= "=" mixture p_edit ;
p_edit          ::= p_loc p_change | p_loc "=" mixture p_change ;
p_change        ::= p_substitution | p_deletion | p_duplication | p_insertion | p_delins
                  | p_frameshift | p_extension | p_n_extension | p_identity | p_unknown | p_repeat ;
p_loc           ::= p_upos ( "_" p_upos )? ;
p_upos          ::= p_pos | p_uncertain ;
p_uncertain     ::= "(" p_pos "_" p_pos ")" | "(" p_pos ")" ;
p_pos           ::= residue number ;
p_substitution  ::= residue ( "^" residue )* ;
p_deletion      ::= "del" ;
p_duplication   ::= "dup" ;
p_insertion     ::= "ins" p_inserted ;
p_delins        ::= "delins" p_inserted ;
p_inserted      ::= p_sequence count? | p_stop_length | p_inserted_range ;
p_stop_length   ::= stop number ;
p_inserted_range ::= p_pos "_" p_pos ;
p_frameshift    ::= p_fs_new? "fs" p_new_stop? ;
p_fs_new        ::= residue | "(" residue ( "^" residue )+ ")" ;
p_new_stop      ::= stop p_stop_at ;
p_stop_at       ::= number | "?" ;
p_extension     ::= residue "ext" p_new_stop ;
p_n_extension   ::= "ext" "-" number ;
p_identity      ::= "=" ;
p_unknown       ::= "?" ;
p_repeat        ::= count other_allele* ;
p_sequence      ::= _aa3_sequence | _aa1_sequence ;
_aa3_sequence   ::= aa3+ stop3? | stop3 ;
_aa1_sequence   ::= aa1+ "*"? ;
residue         ::= aa3 | aa1 | stop ;
stop            ::= "Ter" | "*" ;
stop3           ::= "Ter" | "*" ;
aa3             ::= "Ala" | "Arg" | "Asn" | "Asp" | "Cys" | "Gln" | "Glu" | "Gly" | "His" | "Ile"
                  | "Leu" | "Lys" | "Met" | "Phe" | "Pro" | "Ser" | "Thr" | "Trp" | "Tyr" | "Val"
                  | "Sec" | "Pyl" | "Asx" | "Glx" | "Xle" | "Xaa" ;
aa1             ::= [ACDEFGHIKLMNPQRSTVWYUOBZJX] ;
)BNF";

// One family of nucleotide rules; '@' is the family's name (g, c, n, r) and {TYPE}, {BASE},
// {WHOLE} its coordinate letters, its alphabet and its statements about a whole sequence.
inline constexpr std::string_view family = R"BNF(
# ---- the @ family ----------------------------------------------------------------------------------
@_variants      ::= {TYPE} "." @_list ;
@_list          ::= @_alternative ( "^" @_alternative )* | {WHOLE} ;
@_alternative   ::= @_group ( "(;)" @_group )* ;
@_group         ::= @_allele ( ";" @_allele )* | @_chain ;
@_allele        ::= "[" @_body "]" supernumerary? ;
@_body          ::= @_predicted_list | @_chain ( separator @_chain )* | {WHOLE} ;
@_predicted_list ::= "(" @_chain separator @_chain ( separator @_chain )* ")" ;
@_chain         ::= @_item ( "::" @_item )* ;
@_item          ::= @_edit | @_predicted | @_mixed ;
@_predicted     ::= "(" ( @_edit ( "^" @_edit )* | @_mixed ) ")" ;
@_mixed         ::= "=" mixture @_edit ;
@_edit          ::= @_loc @_change | @_loc "=" mixture @_change ;
@_change        ::= @_substitution | @_deletion | @_duplication | @_insertion | @_delins
                  | @_inversion | @_identity | @_repeat | supernumerary | methylation ;
@_loc           ::= @_upos ( "_" @_upos )? ;
@_upos          ::= @_pos | @_uncertain ;
@_uncertain     ::= "(" @_pos "_" @_pos ")" | "(" @_pos ")" ;
@_substitution  ::= @_was ">" @_now ;
@_was           ::= {BASE} ;
@_now           ::= {BASE} ;
@_deletion      ::= "del" @_stated? ;
@_duplication   ::= "dup" @_stated? ;
@_stated        ::= @_sequence count? ;
@_insertion     ::= "ins" @_inserted ;
@_delins        ::= "delins" @_inserted ;
@_inversion     ::= "inv" ;
@_identity      ::= @_was? "=" ;
@_repeat        ::= ( @_unit+ | count ) other_allele* ;
@_unit          ::= @_sequence count ;
@_inserted      ::= @_piece | "[" @_bracketed ( ";" @_bracketed )* "]" ;
@_bracketed     ::= @_piece | remote ;
@_piece         ::= @_sequence count? | @_loc inverted? ;
@_sequence      ::= {BASE}+ ;
)BNF";

inline void replace_all(std::string& s, std::string_view from, std::string_view to) {
  for (std::size_t i = s.find(from); i != std::string::npos; i = s.find(from, i + to.size())) s.replace(i, from.size(), to);
}

inline std::string family_rules(std::string_view name, std::string_view type, std::string_view base, std::string_view whole) {
  std::string s(family);
  replace_all(s, "{TYPE}", type);
  replace_all(s, "{BASE}", base);
  replace_all(s, "{WHOLE}", whole);
  replace_all(s, "@", name);
  return s;
}

}  // namespace detail

// The whole grammar's text (about 150 rules), as morphe::formal::compile_bnf reads it.
inline const std::string& hgvs_grammar_text() {
  static const std::string text = [] {
    constexpr std::string_view dna = "[ACGTBDHKMNRSVWY]", rna = "[acgubdhkmnrsvwy]";
    std::string t(detail::grammar_head);
    t += detail::family_rules("g", "g_type", dna, "whole");
    t += detail::family_rules("c", "\"c\"", dna, "whole");
    t += detail::family_rules("n", "\"n\"", dna, "whole");
    t += detail::family_rules("r", "\"r\"", rna, "r_whole");
    return t;
  }();
  return text;
}

// The compiled grammar, made once (compile_bnf is pure; the object is immutable after).
inline const morphe::formal::Bnf& hgvs_grammar() {
  static const morphe::formal::Bnf bnf = morphe::formal::compile_bnf(hgvs_grammar_text());
  return bnf;
}

}  // namespace genomics::nomenclature
