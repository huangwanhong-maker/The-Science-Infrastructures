# sequence: Alphabets, Formats and k-mers

*Status: draft (2026-10-24).* The base of the genomics domain: how a genome and its products are written as codes, read from and written to their formats, translated by the genetic codes, cut into k-mers and sampled, and measured along their length. Every other module of the domain computes on what this one gives.

## 1. Concepts

- **Codes, not letters.** A sequence is a vector of bytes: DNA is A = 0, C = 1, G = 2, T = 3, N = 4, so that the complement of x < 4 is 3 − x, two bits hold a base, and the order of codes is the order of the letters; protein is the order of the NCBI matrices, `ARNDCQEGHILKMFPSTWYVBZX*`. Letters exist only at the boundaries, in files and requests.
- **Reading never fails on content.** A reader returns what it could read, with a list of problems giving the line and the reason. A GenBank location is a small context-free language, and is read as one.
- **Translation is a table and its exceptions.** The twenty-seven NCBI translation tables are transcribed from NCBI's own file; a coding sequence's protein follows its `/codon_start`, `/transl_table` and `/transl_except` (selenocysteine at TGA, stops completed by polyadenylation), and is held to RefSeq's own `/translation` on every CDS of fourteen genomes.
- **k-mers name the pieces of a double strand.** A k-mer and its reverse complement have one canonical name; hash64, a bijection of [0, 4^k), orders them at random without collisions; minimizers and syncmers choose a sparse subset that two overlapping sequences share.

## 2. Content

| Header | Content |
|---|---|
| [`alphabet.hpp`](include/genomics/sequence/alphabet.hpp) | DNA, IUPAC and protein codes; complements and reverse complements of codes and of letters |
| [`packed.hpp`](include/genomics/sequence/packed.hpp) | `PackedDna`: two bits a base, runs of N beside the words; the reverse complement, k-mers and counts a word at a time |
| [`genetic_code.hpp`](include/genomics/sequence/genetic_code.hpp) | The 27 NCBI translation tables (gc.prt 4.6); codons, ambiguity resolved when every reading agrees; translation; ORFs on both strands |
| [`fasta.hpp`](include/genomics/sequence/fasta.hpp) | FASTA and FASTQ read and written; the Phred scale; `Problem` |
| [`genbank.hpp`](include/genomics/sequence/genbank.hpp) | GenBank records: LOCUS, DEFINITION, lineage, features with their qualifiers, ORIGIN; locations by recursive descent (join, order, complement, partial ends, sites between bases, remote parts); extraction and `translate_cds` |
| [`kmer.hpp`](include/genomics/sequence/kmer.hpp) | k-mers ≤ 32, their reverse complements and canonical forms; minimap2's hash64; ntHash for any k; (w, k)-minimizers; closed and open syncmers |
| [`composition.hpp`](include/genomics/sequence/composition.hpp) | Base counts, GC in windows, GC skew and its cumulative sum (origin and terminus), CpG observed/expected, k-mer entropy, codon counts, RSCU, Wright's effective number of codons |

## 3. Dependencies

None in the infrastructure; the checks use calx (`verification`, `random`).

## 4. Checks

`check_genomics_sequence`, `check_genomics_formats`, `check_genomics_kmers`:

| Check | What it verifies |
|---|---|
| alphabets | Letters and codes round-trip; U reads as T; the complement is an involution; the fifteen IUPAC sets and their complements; protein codes in the matrices' order |
| packed DNA | Against codes, base by base, count by count, k-mer by k-mer and reverse complement by reverse complement, with runs of N of 1, 31, 32 and 33 bases; a megabase in a quarter of a byte a base |
| genetic codes | The 27 tables well formed; the standard code's stops and starts; the vertebrate mitochondrial AGA, AGG, TGA and ATA; the bacterial starts; ambiguous codons; a peptide back-translated and translated again; frames |
| ORFs | On random sequence, on both strands: each starts at ATG, ends at its first stop in frame, and begins at the first start after the previous stop |
| composition | GC; CpG o/e near 1 in i.i.d. sequence and 0 without CpG; k-mer entropy; the cumulative skew's extremes at a planted origin and terminus; RSCU 1 and N_c = 61 for even usage; N_c = 20 for one codon per amino acid |
| FASTA, FASTQ | Comments, wrapping, descriptions, CRLF; four-line and multi-line records; the faults found on their lines; round trips; Phred |
| locations | Spans, partial ends, sites between bases, complement of a join, mixed strands, order, remote parts, and the errors of each |
| GenBank | A circular genome written by the check, with genes across the origin, on the complement strand and from GTG, read whole and translated as its `/translation`; and, with the data cache, **every CDS of the fourteen RefSeq genomes of `ncbi-genomes` (4,522 proteins, the three selenoproteins of *E. coli* among them) translated exactly as RefSeq's own `/translation`** |
| k-mers | Encoding, reverse complements and canonical forms for k = 1 to 32; every k-mer free of N visited once; hash64 a bijection for k ≤ 8; ntHash rolled against ntHash computed, k = 1 to 100, and its canonical value strand-symmetric |
| minimizers, syncmers | Against their definitions by brute force (ties, palindromes, runs of N); mirrored on the reverse strand; densities 2/(w + 1), 2/(k − s + 1) and 1/(k − s + 1) on a random megabase |

## 5. References

- Cock, P. J. A., Fields, C. J., Goto, N., Heuer, M. L. and Rice, P. M. (2010), The Sanger FASTQ file format, *Nucleic Acids Research* 38, 1767–1771.
- Cornish-Bowden, A. (1985), Nomenclature for incompletely specified bases in nucleic acid sequences, *Nucleic Acids Research* 13, 3021–3030.
- Edgar, R. (2021), Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences, *PeerJ* 9, e10805.
- Gardiner-Garden, M. and Frommer, M. (1987), CpG islands in vertebrate genomes, *Journal of Molecular Biology* 196, 261–282.
- Grigoriev, A. (1998), Analyzing genomes with cumulative skew diagrams, *Nucleic Acids Research* 26, 2286–2290.
- INSDC, *The DDBJ/ENA/GenBank Feature Table Definition*, §3.4 (locations) and §7 (qualifiers); NCBI, *GenBank Release Notes*, §3.4.
- Li, H. (2018), Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics* 34, 3094–3100.
- Mohamadi, H., Chu, J., Vandervalk, B. P. and Birol, I. (2016), ntHash: recursive nucleotide hashing, *Bioinformatics* 32, 3492–3494.
- NCBI, *The Genetic Codes* (gc.prt, version 4.6), https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt, retrieved 2026-10-24.
- Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M. and Yorke, J. A. (2004), Reducing storage requirements for biological sequence comparison, *Bioinformatics* 20, 3363–3369.
- Schleimer, S., Wilkerson, D. S. and Aiken, A. (2003), Winnowing: local algorithms for document fingerprinting, *SIGMOD*, 76–85.
- Sharp, P. M., Tuohy, T. M. F. and Mosurski, K. R. (1986), Codon usage in yeast, *Nucleic Acids Research* 14, 5125–5143.
- Wright, F. (1990), The 'effective number of codons' used in a gene, *Gene* 87, 23–29.

## 6. Open Questions

1. **hash64 at k ≤ 15.** Its last step vanishes modulo 2^{2k}, and the density of minimizers rises 2–3 per cent above 2/(w + 1) (`kmer.hpp`, Remark). A stronger invertible mix for short k would cost a step; whether to adopt one, at the price of differing from minimap2's minimizers, is open.
2. **Large genomes.** `PackedDna` holds a genome in memory; chromosomes of a few hundred megabases fit, a whole human genome in about 750 MB. Memory-mapped packed files would serve larger ones.
