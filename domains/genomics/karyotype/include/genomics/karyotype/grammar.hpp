// The grammar of ISCN, the International System for Human Cytogenomic Nomenclature (ISCN 2024,
// Hastings, Moore and Chia, eds., Karger 2024; ISCN 2020, McGowan-Jordan, Hastings and Moore,
// eds., Karger 2020), written once as BNF over bytes and read by morphe's BNF reader
// (morphe/formal/bnf.hpp: Earley's parser, named trees, the furthest viable prefix and what may
// follow it).
//
// The language covered:
// - Karyotypes in the short system: the chromosome count (a number, a range 45~48 or 44-45, a
//   ploidy <3n>), the sex chromosomes (XX, XY, X, XXY, …, U where the sex is undisclosed, 2024;
//   c when constitutional, ? when uncertain), then the abnormalities separated by commas:
//   numerical (+21, −7, +X, +21c), structural (add, del, der, dic, dup, fis, fra, hsr, i, idic,
//   ider, ins, inv, neo, qdp, r, rec, rob, t, tas, trc, trp; with the prefixes psu, inv, dir, dic,
//   trc and der; with ? before them; a multiplier x2; c and the inheritance dn, mat, pat, inh,
//   dmat, dpat, dinh after them), markers (+mar, +2mar, +mar1, +r), double minutes (3~9dmin),
//   inc, and alternatives separated by "or".
// - Breakpoints by arm, region, band and sub-band (p13, q21.1, p36.33), pter, qter, cen, the
//   centromeric p10 and q10, uncertain bands (q2?, q?), intervals of bands (q21~q24, q11.2-13.1).
// - Cell lines with their counts of cells ([20], [cp10] for a composite karyotype), separated by
//   / (clones, mosaics) or // (chimeras, donor and recipient), mos and chi before them; idem, sl,
//   sdl and sdl1 standing for an earlier clone.
// - The detailed system: each derivative chromosome from end to end, segments a→b (or a->b)
//   joined by :: (break and reunion), : a break, ? material of unknown origin, a ring opened and
//   closed by ::; with nucleotide positions for sequencing results (Ordulu et al. 2014).
// - Other techniques after a period: arr (microarrays: 1p36.33p36.32(849,467_2,912,345)x1, x2~3,
//   [0.5], hmz, whole chromosomes (X,1-22)x2, the extended form with the flanking probes), seq
//   (sequencing: structural abnormalities in either system, regions as for arr, junctions as
//   1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650), (2,3,8)cx), and ish, nuc ish and rsa, which are
//   recognised and kept as text but not read further.
//
// The grammar is the infrastructure's reading of the standard: the rules come from ISCN 2024's
// summary of revisions (Hastings et al. 2025, Cytogenet Genome Res 165, 1–8), from the comparison
// of ISCN 2020 with 2016 (Miller and Madan 2021, E.C.A. Newsletter 47, 2–11), from the examples of
// ISCN 2009 and 2013 quoted by the Karyotype Ontology (Warrender, Moorman and Lord 2019,
// Bioinformatics 35, 5264–5270), from Ordulu et al. (2014, Am. J. Hum. Genet. 94, 695–709) for
// seq, and from Conlin et al. (2025, Clin. Chem. 71, 45–55) for arr; the theory note
// (theory/karyotypes.md, §2) says which rule comes from which. Multiplication signs × and x,
// arrows → and ->, minus signs − and -, and tildes ∼ and ~ are both accepted.
//
// Complexity: Earley's parser, O(n²) on these unambiguous rules for a text of n bytes (a few
// hundred for a karyotype); the grammar is compiled once.

#pragma once

#include <morphe/formal/bnf.hpp>

#include <string_view>

namespace genomics::karyotype {

inline constexpr std::string_view iscn_grammar = R"BNF(
# ---- a description: a karyotype, then other techniques after periods --------------------------
description     ::= _s karyotype ( "." _technique )* _s
                  | _s _technique ( "." _technique )* _s ;
karyotype       ::= ( mosaic _sp )? cell_line ( line_sep cell_line )* ;
mosaic          ::= "mos" | "chi" ;
line_sep        ::= _s ( "//" | "/" ) _s ;
cell_line       ::= count ( _c _head )? ( _c _items )? _s cells? ;
count           ::= number ( _tilde number | "-" number )? ploidy? | "?" ;
ploidy          ::= "<" number "n>" ;
_head           ::= sex | stemline ;
sex             ::= sex_chromosomes sex_flag? | "U" ;
sex_chromosomes ::= [XY]+ ;
sex_flag        ::= "?c" | "c" | "?" ;
stemline        ::= "idem" | "sl" | "sdl" number? ;
_items          ::= item ( _c item )* ;
item            ::= _abnormality alternative* ;
alternative     ::= _sp "or" _sp _abnormality ;
_abnormality    ::= numeric | structural | marker | ring_marker | dmin | inc ;
cells           ::= "[" ( composite number | number ) "]" ;
composite       ::= "cp" ;

# ---- numerical abnormalities ---------------------------------------------------------------------
numeric         ::= sign chromosome flag* ;
sign            ::= "+" | "-" | "\xE2\x88\x92" ;
chromosome      ::= [1-9] [0-9]? | "X" | "Y" ;

# ---- structural abnormalities ---------------------------------------------------------------------
structural      ::= sign? uncertain? prefix* head component* multiplier? flag* ;
uncertain       ::= "?" ;
prefix          ::= ( "psu" | "inv" | "dir" | "dic" | "trc" | "der" ) _sp ;
head            ::= abbr chroms _detail? ;
component       ::= abbr chroms _detail? | rec_part ;
rec_part        ::= ( "dup" | "del" ) "(" chromosome arm ")" ;
abbr            ::= "add" | "del" | "der" | "dic" | "dup" | "fis" | "fra" | "hsr" | "idic" | "ider"
                  | "ins" | "inv" | "i" | "neo" | "qdp" | "rec" | "rob" | "r" | "tas" | "trc"
                  | "trp" | "t" ;
chroms          ::= "(" chrom_ref ( ";" chrom_ref )* ")" ;
chrom_ref       ::= ( chromosome | unknown ) homologue? ;
unknown         ::= "?" ;
homologue       ::= "*" ;
_detail         ::= "(" group ( ";" group )* ")" | "(" long_body ")" ;
group           ::= breakpoint+ | unknown ;
breakpoint      ::= terminal | arm band_digits uncertain? ( range_sep band_end )? | arm uncertain ;
terminal        ::= "pter" | "qter" | "cen" ;
arm             ::= "p" | "q" ;
band_digits     ::= [0-9]+ ( "." [0-9]+ )? ;
range_sep       ::= _tilde | "-" ;
band_end        ::= arm? band_digits ;
multiplier      ::= _times [0-9]+ ;
flag            ::= _s ( "c" | "dn" | "mat" | "pat" | "inh" | "dmat" | "dpat" | "dinh" ) ;

# ---- the detailed system ---------------------------------------------------------------------------
long_body       ::= long_chromosome ( _s ";" _s long_chromosome )* ;
long_chromosome ::= lf_end? ( lf_unknown lf_join )* lf_segment ( lf_join _lf_part )* lf_end? ;
lf_end          ::= "::" | ":" ;
lf_join         ::= _s "::" _s ;
_lf_part        ::= lf_segment | lf_unknown | lf_bases | lf_region ;
lf_unknown      ::= "?" ;
lf_bases        ::= [ACGTN]+ ;
lf_segment      ::= lf_point _arrow lf_point ;
lf_point        ::= chromosome? ( terminal | arm band_digits ( range_sep band_end )? ) lf_coord? ;
lf_region       ::= chromosome? arm band_digits arm band_digits lf_coord multiplier ;
lf_coord        ::= "(" coord ")" | "(" coord _tilde coord? ")" | "(" coord _dash coord ")" ;
coord           ::= [0-9]+ ( "," [0-9] [0-9] [0-9] )* ;

# ---- markers, rings of unknown origin, double minutes ---------------------------------------------
marker          ::= sign? marker_count? "mar" marker_number? multiplier? flag* ;
marker_count    ::= [0-9]+ ;
marker_number   ::= [0-9]+ ;
ring_marker     ::= sign? marker_count? "r" multiplier? flag* ;
dmin            ::= ( number ( _tilde number )? )? "dmin" ;
inc             ::= "inc" ;

# ---- other techniques ---------------------------------------------------------------------------------
_technique      ::= arr | seq | ish | rsa ;
arr             ::= ( "arr" | "sseq" ) build? _s arr_entry ( _c arr_entry )* ;
seq             ::= "seq" build? _s _seq_entry ( _c _seq_entry )* ;
_seq_entry      ::= structural | arr_entry | junction ;
build           ::= "[" [A-Za-z0-9./_]+ "]" ;
arr_entry       ::= ( arr_region | arr_whole | arr_complex ) arr_flag* ;
arr_region      ::= chromosome arm band_digits ( arm band_digits )? "(" ( arr_span | arr_extended ) ")" copies? ;
arr_span        ::= coord _range coord ;
arr_extended    ::= coord copies "," coord _range coord copies "," coord copies ;
arr_whole       ::= "(" whole_ref ( "," whole_ref )* ")" copies ;
whole_ref       ::= chromosome ( "-" chromosome )? | "XY" ;
arr_complex     ::= "(" chromosome ( "," chromosome )* ")" "cx" ;
copies          ::= _s _times _s [0-9]+ ( _tilde [0-9]+ )? proportion? ;
proportion      ::= "[" [0-9]+ ( "." [0-9]+ )? "]" ;
arr_flag        ::= _s ( "hmz" | "htz" | "cth" | "cha" | "c" | "dn" | "mat" | "pat" | "inh" | "dmat" | "dpat" | "dinh" | "mos" ) ;
junction        ::= junction_end "::" ( lf_bases "::" )? junction_end ;
junction_end    ::= chromosome arm band_digits strand? lf_coord ;
strand          ::= "(+)" | "(-)" | "(\xE2\x88\x92)" ;
ish             ::= ( "nuc" _sp )? "ish" _s opaque ;
rsa             ::= "rsa" build? _s opaque ;
opaque          ::= _opaque_char+ ;
_opaque_char    ::= [^.] | "." [0-9] ;

# ---- tokens -----------------------------------------------------------------------------------------------
number          ::= [0-9]+ ;
_c              ::= _s "," _s ;
_s              ::= " "* ;
_sp             ::= " "+ ;
_tilde          ::= "~" | "\xE2\x88\xBC" ;
_dash           ::= "-" | "\xE2\x80\x93" ;
_range          ::= "_" | "-" | "\xE2\x80\x93" ;
_times          ::= "x" | "\xC3\x97" ;
_arrow          ::= _s ( "->" | "\xE2\x86\x92" ) _s ;
)BNF";

// The grammar compiled once (an immutable value, built on first use).
inline const morphe::formal::Bnf& iscn_bnf() {
  static const morphe::formal::Bnf g = morphe::formal::compile_bnf(iscn_grammar);
  return g;
}

}  // namespace genomics::karyotype
