# classify: Taxonomic Classification of Reads

*Status: draft (2026-10-24).* Reads assigned to the taxa they come from by exact k-mers and the lowest common ancestor, Kraken's method (Wood and Salzberg 2014): the taxonomy built from the lineages of GenBank records, the database that maps every canonical k-mer of the reference genomes to the LCA of the genomes containing it, and the heaviest root-to-leaf path of a read's hits. It serves use, as the first look at what a sample holds, and learning, as an algorithm on a semilattice. The theory is [theory/signals-classification-and-comparison.md](../theory/signals-classification-and-comparison.md), §2.

## 1. Concepts

- **The taxonomy is the trie of the records' lineages.** A GenBank record names its organism and the taxa above it; the trie of those lists under one root is the tree, the organisms its leaves. Nothing else is typed in. A flat file's lineage carries no ranks, so a taxon's depth is its place along the lineages, not a Linnaean rank.
- **The lowest common ancestor is the meet of a semilattice.** It is associative, commutative and idempotent, so the taxon a k-mer receives, the LCA of the genomes that contain it, is the same in any order of building.
- **A read goes to the heaviest root-to-leaf path** of the tree of its hit taxa and their ancestors, a path weighed by the hits along it; ties go to the LCA of their leaves. As the hits are positive, the heaviest paths end at hit taxa, so the tree need not be built.
- **Exact k-mers classify as specifically as the database allows, and no more.** Reads of a genome in the database go to its leaf, or to an ancestor where relatives share their k-mers: of 8400 reads of the fourteen cached genomes, 0.9996 land on their own lineage, 0.963 at their leaf. With its genome left out, a read goes to the nearest relatives present: E. coli's, phiX174's and SARS-CoV-2's stay unclassified, lambda's go to E. coli K-12 (with which it shares 31-mers) or nowhere, and a primate mitochondrion's stay among the primates, 0.88 of them within the clade of the deepest ancestor still in the database (mostly on a relative's own leaf, which Kraken's rule cannot avoid), the rest at a higher ancestor.

## 2. Content

| Header | Content |
|---|---|
| [`taxonomy.hpp`](include/genomics/classify/taxonomy.hpp) | `TaxonId`, `ROOT`, `NO_TAXON`; `Taxon{name, parent, depth}`; `Taxonomy` with `child(parent, name)`, `add(lineage, organism)`, `add(GenBank)`, `find(name)`, `lca(a, b)`, `above(a, b)`, `path(t)`, `lineage(t)`, `size()` |
| [`kraken.hpp`](include/genomics/classify/kraken.hpp) | `Database(taxonomy, k)` with `add(genome, taxon)`, `lookup(canonical)`, `prefetch`, `size`, `genomes`, `memory_bytes`, `for_each` (on kmers' `KmerMap`); `ClassifyOptions{min_hits, confidence}`; `Classification{taxon, kmers, hits, clade_hits}`; `heaviest_path(taxonomy, hits)`, `clade_hits(taxonomy, label, hits)`; `classify(db, read, options)`, `classify_all(db, reads, pool, options)` |

## 3. Dependencies

`genomics::kmers` (the open-addressing table and its mix), `genomics::sequence` (canonical k-mers, the GenBank reader); calx `scheduler` (the pool). The checks use `genomics::simulate` and calx `random`.

## 4. Checks

| Check | What it verifies |
|---|---|
| `taxonomy` (11) | On 30 random taxonomies, the LCA of 9000 pairs against the deepest taxon of both root paths; its laws (commutative, associative, idempotent, the root absorbing); `above` against the root path; root paths through parents. A lineage name wrapped over two lines of a GenBank record read as one name. The fourteen cached genomes' lineages joined into a tree of 73 taxa: human and Neanderthal two leaves of Homo, chimpanzee and bonobo of Pan, the orangutans of Pongo; the LCAs of human and chimpanzee (Hominidae), of human and gibbon or macaque (Catarrhini), of macaque and baboon (Cercopithecinae), of E. coli and lambda (the root), of the phages (Viruses) |
| `kraken` (13) | The database against brute force on 12 sets of 10 genomes sharing segments on both strands (k = 7 to 11): every canonical k-mer on the LCA of the genomes containing it, absent ones on none, the same built in the opposite order. The heaviest path against the classification tree built explicitly, on 3000 random sets of hits; the confidence climb by hand. With the cache: a database of the fourteen genomes at k = 31 (4 772 316 k-mers, 0.5 s); 600 simulated reads of each, 0.9996 on their lineage (at least 0.99 asserted), 0.997 at their leaf for the genomes without a close relative (at least 0.98); the same on one thread and four. Each genome left out in turn: E. coli, phiX174 and SARS-CoV-2 unclassified (at least 0.99); lambda unclassified (0.892, at least 0.85) or on E. coli and nowhere else; a primate never off the primates, and 0.880 of its classified reads (at least 0.85) within the clade of its deepest ancestor still in the database |

## 5. References

- Wood, D. E. and Salzberg, S. L. (2014), Kraken: ultrafast metagenomic sequence classification using exact alignments, *Genome Biology* 15, R46.
- Aho, A. V., Hopcroft, J. E. and Ullman, J. D. (1973), On finding lowest common ancestors in trees, *Proceedings of the Fifth ACM Symposium on Theory of Computing*, for the problem; the climb by depths used here is its simplest solution.
- Wood, D. E., Lu, J. and Langmead, B. (2019), Improved metagenomic analysis with Kraken 2, *Genome Biology* 20, 257.
- Lu, J., Breitwieser, F. P., Thielen, P. and Salzberg, S. L. (2017), Bracken: estimating species abundance in metagenomics data, *PeerJ Computer Science* 3, e104.

## 6. Open Questions

1. **Trees.** The taxonomy is a rooted tree with LCAs; `msa` keeps a guide tree of its own and evolution's `phylogeny` the full ones. Rooted trees, their LCAs (by Euler tours and range minima, for large taxonomies) and their semilattice belong in a foundation (the domain README's open question 3).
2. **The NCBI taxonomy.** A database of thousands of genomes wants NCBI's taxonomy with its ranks (`nodes.dmp`, `names.dmp`), fetched by the data library; the trie of lineages is exact for the records at hand but knows no ranks.
3. **Minimizers and memory.** Kraken stores k-mers sorted by minimizer for locality, and Kraken 2 (Wood, Lu and Langmead 2019) stores minimizers alone in a compact hash table; for a database of many bacteria the table here (12 bytes a slot) should follow.
4. **Abundances.** Reads at internal taxa can be redistributed to species by their genomes' unique k-mers, as Bracken does (Lu et al. 2017); it belongs with the EM of `expression`.
