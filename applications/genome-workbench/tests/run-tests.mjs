// Checks of what is the Genome Workbench's own: how it writes lengths, positions and ticks; and,
// when the engine is built, its answers held to the protocol of the README (§3): a genome's
// summary and tracks (phiX174's length, GC and genes; the origin of Escherichia coli from the
// cumulative skew within a few hundred bases of oriC), the search of a pattern on both strands
// against a count by brute force, a feature's protein as RefSeq translates it, the gene finder
// against RefSeq's annotation, every schedule of an alignment held to Gotoh's optimum, the batch
// kernels to one another, the consensus of a partial order to its truth, and a job whose answer
// is the same however it is cut into shards and spread over workers; and the group grammars: the
// draft's orbit-only scan accepting a word outside the language, Burnside's count of rule orbits,
// the canonical search of a monotone grammar, the codes' stabilisers, Σ χ(1)² = |G|, hairpins
// against their definition, and Chargaff's second rule on Escherichia coli. The science is
// checked with the genomics and evolution domains.
// Run from the repository root:  node applications/genome-workbench/tests/run-tests.mjs

import { spawnSync } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';

import { suite, ROOT, CACHE } from '../../../development/tools/checks.mjs';
import * as fmt from '../src/format.js';
import { readHash, writeHash } from '../src/state.js';

const { assert, check, skip, done } = suite();
const APP = resolve(ROOT, 'applications/genome-workbench');

// ---- writing --------------------------------------------------------------------------------------
{
  assert('lengths: 5,386 bp, 48.5 kb, 4.64 Mb', fmt.length(5386) === '5,386 bp' && fmt.length(48502) === '48.5 kb' && fmt.length(4641652) === '4.64 Mb', `${fmt.length(5386)} ${fmt.length(48502)} ${fmt.length(4641652)}`);
  assert('positions are written 1-based, intervals as b..e', fmt.position(0) === '1' && fmt.interval(99, 250) === '100..250');
  assert('ticks step by 1, 2 or 5 times a power of ten', [fmt.tickStep(4641652, 12), fmt.tickStep(16569, 8), fmt.tickStep(5386, 14)].every((s) => /^(1|2|5)0*$/.test(String(s))), `${fmt.tickStep(4641652, 12)} ${fmt.tickStep(16569, 8)}`);
  assert('a megabase genome is labelled in Mb', fmt.tickLabel(3800000, 200000, 4641652) === '3.8 Mb' && fmt.tickLabel(2000, 1000) === '2 kb', fmt.tickLabel(3800000, 200000, 4641652));
  const s = { view: 'align', genome: 'ncbi:NC_000913.3', params: { begin: 100, end: 900, a: 'ACGT' } };
  const back = readHash(writeHash(s));
  assert('the address holds the view, the genome and its parameters', back.view === 'align' && back.genome === s.genome && back.params.begin === '100' && back.params.a === 'ACGT', writeHash(s));
  const tb = readHash('#toolbox&tool=align-msa'), old = readHash('#view=trees&gene=ND5');
  assert('the view is the address\'s first word, and the older view= is still read', writeHash({ view: 'toolbox', params: { tool: 'align-msa' } }) === '#toolbox&tool=align-msa' && tb.view === 'toolbox' && tb.params.tool === 'align-msa' && !('toolbox' in tb.params) && old.view === 'trees' && old.params.gene === 'ND5');
}

// ---- the toolbox: tools are declarations -------------------------------------------------------
{
  await import('../src/tools/index.js');
  const registry = await import('../src/core/registry.js');
  const { match } = await import('../src/core/palette.js');
  const groups = new Set(registry.TOOL_GROUPS.map(([g]) => g));
  const TYPES = new Set(['sequence', 'sequences', 'file', 'genome', 'number', 'range', 'choice', 'select', 'switch', 'text']);
  const ids = registry.tools.map((t) => t.id);
  const bad = registry.tools.filter((t) => !groups.has(t.group) || typeof t.run !== 'function' || typeof t.render !== 'function' || !(t.inputs || []).every((i) => i.group || TYPES.has(i.type)));
  assert('every tool is declared whole: a known group, inputs of known types, run and render', ids.length >= 25 && new Set(ids).size === ids.length && bad.length === 0, bad.map((t) => t.id).join(', '));
  assert('every kind of file a tool takes is one the library recognises', registry.tools.flatMap((t) => (t.inputs || []).flatMap((i) => i.kinds || [])).every((k) => ['fasta', 'fastq', 'genbank', 'vcf', 'newick', 'sam'].includes(k)));
  assert('a file of reads finds the tools that take it', ['reads-qc', 'reads-map', 'reads-call', 'assembly-reads'].every((id) => registry.toolsFor('fastq').some((t) => t.id === id)));
  const m = match('msa', 'Multiple alignment msa muscle'), f = match('mlt', 'Multiple alignment');
  assert('the palette finds a word directly, and letters scattered in order, but not noise', m && m.at.join() === '19,20,21' && f && match('tree', 'A small genome to assemble, 30 kb') === null);
}

// ---- the engine, when it is built --------------------------------------------------------------
{
  const status = spawnSync('python', [join(APP, 'engine.py'), 'status'], { cwd: ROOT, encoding: 'utf8' });
  let state = null;
  try { state = JSON.parse(status.stdout); } catch { /* not built */ }
  if (!state?.current) {
    skip('the engine is not built from the present sources: its answers are not checked');
  } else {
    // Requests through the cluster, as the server makes them.
    const genomes = resolve(CACHE, 'ncbi-genomes');
    const script = `
import json, sys
sys.path.insert(0, ${JSON.stringify(APP)})
import cluster, engine
kind = (engine.read_build() or {}).get("toolchain") or "native"
def g(acc): return {"file": engine.path_for(kind, engine.Path(${JSON.stringify(genomes)}) / (acc + ".gb"))}
c = cluster.Cluster(local=2)
out = {}
sim = {"simulated": {"length": 50000, "gc": 0.5, "seed": 4}}
try:
    if (engine.Path(${JSON.stringify(genomes)}) / "NC_001422.1.gb").exists():
        out["phix"] = c.request({"op": "genome.summary", "genome": g("NC_001422.1")})
        out["phix_tracks"] = c.request({"op": "genome.tracks", "genome": g("NC_001422.1"), "bins": 50})
        out["phix_search"] = c.request({"op": "genome.search", "genome": g("NC_001422.1"), "pattern": "GGCTA"})
        cds = [f for f in out["phix_tracks"]["features"] if f["key"] == "CDS"]
        out["phix_feature"] = c.request({"op": "genome.feature", "genome": g("NC_001422.1"), "index": cds[0]["i"]})
        out["phix_letters"] = c.request({"op": "genome.tracks", "genome": g("NC_001422.1"), "begin": 0, "end": 5386, "bins": 10, "letters_max": 6000})["letters"]
        # The toolbox's operations, on phiX174's file and its letters, and on reads drawn from it.
        phix = g("NC_001422.1")["file"]
        letters = out["phix_letters"]
        t = {}
        t["stats"] = c.request({"op": "tool.seq.stats", "fasta_file": phix})
        rc = c.request({"op": "tool.seq.revcomp", "text": ">x\\n" + letters, "what": "revcomp"})["records"][0]["seq"]
        t["revcomp"] = [rc == letters, c.request({"op": "tool.seq.revcomp", "text": ">x\\n" + rc, "what": "revcomp"})["records"][0]["seq"] == letters]
        t["translate"] = c.request({"op": "tool.seq.translate", "text": ">cds\\n" + out["phix_feature"]["sequence"], "table": 11, "frames": 1})["records"][0]["protein"]
        t["find"] = c.request({"op": "tool.seq.find", "fasta_file": phix, "patterns": ["GGCTA", "GAATTC"], "both": True})
        t["orfs"] = c.request({"op": "tool.seq.orfs", "fasta_file": phix, "table": 11, "min_codons": 60})
        t["kmers"] = c.request({"op": "tool.seq.kmers", "fasta_file": phix, "k": 11, "top": 5})
        t["evalue"] = [c.request({"op": "tool.stats.evalue", "score": s, "matrix": "BLOSUM62", "gap_open": 11, "gap_extend": 1, "m": 300, "n": 1e8}) for s in (40, 60, 80)]
        t["tree"] = c.request({"op": "tool.tree.read", "text": "((a:1,b:2):1,c:4,(d:1,e:1):2);", "root": "midpoint"})
        reads = c.request({"op": "align.batch_input", "genome": g("NC_001422.1"), "targets": 2000, "length": 150, "seed": 3})["targets"]
        import os, shutil, tempfile
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, "phix_reads.fasta")
            with open(path, "w", encoding="utf-8", newline="\\n") as f:
                f.write("".join(f">r{i}\\n{s}\\n" for i, s in enumerate(reads)))
            rf = engine.path_for(kind, engine.Path(path))
            t["qc"] = c.request({"op": "tool.reads.qc", "reads_file": rf})
            t["map"] = c.request({"op": "tool.reads.map", "genome": g("NC_001422.1"), "reads_file": rf, "kind": "short"})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        out["tools"] = t
    if (engine.Path(${JSON.stringify(genomes)}) / "NC_000913.3.gb").exists():
        out["ecoli"] = c.request({"op": "genome.tracks", "genome": g("NC_000913.3"), "bins": 100, "orfs": False})
        out["ecoli_genes"] = c.request({"op": "genome.annotate", "genome": g("NC_000913.3"), "what": "genes"})
    pair = c.request({"op": "align.random", "length": 1500, "divergence": 0.1, "seed": 5})
    out["schedules"] = {k + ":" + m: c.request({"op": "align.pair", "a": pair["a"], "b": pair["b"], "algorithm": k, "mode": m})["reference"]
                        for k, m in [("gotoh", "global"), ("hirschberg", "global"), ("banded", "global"), ("wfa", "global"), ("biwfa", "global"), ("striped", "local"), ("gotoh", "local")]}
    out["compare"] = c.request({"op": "align.compare", "a": pair["a"], "b": pair["b"], "mode": "global"})
    batch = c.request({"op": "align.batch_input", "genome": sim, "targets": 3000, "seed": 2})
    out["batch"] = {k: c.request({"op": "align.batch", "query": batch["query"], "targets": batch["targets"], "kernel": k})["scores"] for k in ("gotoh", "striped", "inter")}
    out["poa"] = c.request({"op": "align.poa", "length": 150, "reads": 10, "error": 0.08, "seed": 2})
    out["job5"] = c.job("reads.sample", {"genome": sim, "reads": 6000, "shards": 5, "seed": 9})
    out["job13"] = c.job("reads.sample", {"genome": sim, "reads": 6000, "shards": 13, "seed": 9})
    # Group grammars: the draft's scan, Burnside's count, a monotone search, the codes, characters.
    gr = {}
    hp = [["S", "A P T | T P A | C P G | G P C"], ["P", "A Q T | T Q A | C Q G | G Q C"], ["Q", "A Q T | T Q A | C Q G | G Q C | A L T | T L A | C L G | G L C"], ["L", "N N N | N N N N | N N N N N"], ["N", "A | C | G | T"]]
    split = lambda rules: [[l, alt.strip()] for l, r in rules for alt in r.split("|")]
    gr["scan"] = c.request({"op": "grammar.parse", "alphabet": "ab", "rules": [["S", "a a"], ["S", "b b"]], "group": [["ba", False]], "word": "ab"})
    gr["brackets"] = c.request({"op": "grammar.parse", "alphabet": "()", "rules": split([["S", "( ) | ( S ) | S S"]]), "group": [[")(", True]], "word": "(()())()"})
    gr["hairpin"] = c.request({"op": "grammar.parse", "alphabet": "ACGT", "rules": split(hp), "group": [["TGCA", False], ["ACGT", True]], "word": "GCATTCGTGC"})
    gr["monotone"] = c.request({"op": "grammar.parse", "alphabet": "abcd", "rules": split([["S", "a S C D | b S C D | a C D | b C D"], ["D C", "C D"], ["a C", "a c"], ["b C", "b c"], ["c C", "c c"], ["c D", "c d"], ["d D", "d d"]]), "group": [["bacd", False]], "word": "abccdd"})
    gr["code1"] = c.request({"op": "grammar.code", "table": 1})
    gr["code2"] = c.request({"op": "grammar.code", "table": 2})
    gr["chars"] = {name: c.request({"op": "grammar.characters", "group": name}) for name in ["cyclic 4", "klein", "dihedral 4", "quaternion", "symmetric 4", "alternating 5"]}
    gr["nucleotides"] = c.request({"op": "grammar.nucleotides", "base_point": "T"})
    if (engine.Path(${JSON.stringify(genomes)}) / "NC_001422.1.gb").exists():
        gr["hairpins"] = c.request({"op": "grammar.hairpins", "fasta_file": g("NC_001422.1")["file"], "begin": 0, "end": 300, "min_stem": 4, "min_loop": 3, "max_loop": 8, "max_span": 80})
    if (engine.Path(${JSON.stringify(genomes)}) / "NC_000913.3.gb").exists():
        gr["symmetry"] = c.request({"op": "grammar.symmetry", "genome": g("NC_000913.3"), "k": 4})
    out["grammar"] = gr
finally:
    c.stop()
print(json.dumps(out))
`;
    const r = spawnSync('python', ['-c', script], { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 26 });
    let out = null;
    try { out = JSON.parse(r.stdout); } catch { assert('the engine answers through the cluster', false, (r.stderr || r.stdout).slice(-400)); }
    if (out) {
      if (out.phix) {
        const c = out.phix.contigs[0];
        assert('phiX174: 5,386 bases, circular, GC 44.8%, eleven CDSs, genetic code 11', out.phix.length === 5386 && c.circular && Math.abs(c.gc - 0.44764) < 1e-4 && c.features.CDS === 11 && out.phix.table === 11, JSON.stringify(c.features));
        const t = out.phix_tracks;
        assert('tracks: fifty bins of GC, skew, CpG and entropy, and the cumulative skew at their 51 edges', t.gc.length === 50 && t.skew.length === 50 && t.cpg.length === 50 && t.cumulative.length === 51);
        const letters = out.phix_letters;
        // GGCTA on either strand: GGCTA, or its reverse complement TAGCC, on the forward strand.
        let sites = 0;
        for (let i = 0; i + 5 <= letters.length; i++) sites += ['GGCTA', 'TAGCC'].includes(letters.slice(i, i + 5));
        assert('search: GGCTA on both strands of phiX174, as counted by brute force on the letters', out.phix_search.total === sites && sites > 0, `${out.phix_search.total} and ${sites}`);
        assert('a CDS of phiX174 translated as RefSeq translates it', out.phix_feature.agrees === true && out.phix_feature.translation.startsWith('M'), out.phix_feature.name);

        // The toolbox's operations, held to closed forms and to counts made here on the letters.
        const tools = out.tools;
        const gc = [...letters].filter((x) => x === 'G' || x === 'C').length / letters.length;
        const st = tools.stats.records[0];
        assert('tool.seq.stats: phiX174 read from its GenBank file, its length and GC those counted on its letters, N50 its length', st.length === 5386 && Math.abs(st.gc - gc) < 1e-5 && tools.stats.total.n50 === 5386, `${st.length} ${st.gc} ${gc}`);
        assert('tool.seq.revcomp: the reverse complement differs, and the reverse complement of it is the sequence again', tools.revcomp[0] === false && tools.revcomp[1] === true);
        const tr = out.phix_feature.translation, pr = tools.translate;
        assert('tool.seq.translate: frame 1 of a CDS of phiX174 under code 11 is RefSeq\'s /translation, then its stop', pr.slice(1).replace(/\*$/, '') === tr.slice(1) && pr.endsWith('*') && !pr.slice(0, -1).includes('*'), `${pr.slice(0, 20)}… ${pr.length} and ${tr.length}`);
        const byMotif = Object.fromEntries(tools.find.by_motif.map((x) => [x.motif, x]));
        assert('tool.seq.find: GGCTA on both strands of the file by one automaton, as many as brute force counts on the letters; GAATTC, a palindrome, once a site', byMotif.GGCTA.count === sites && byMotif.GAATTC.reverse === 0, `${byMotif.GGCTA.count} and ${sites}`);
        const orfsOk = tools.orfs.orfs.every((o) => (o.end - o.begin) === 3 * o.codons && o.codons >= 60 && o.protein.startsWith('M') && !o.protein.includes('*') && o.protein.length === o.codons - 1);
        assert('tool.seq.orfs: every ORF of phiX174 whole codons from ATG to its first stop, its protein one residue a codon but the stop', tools.orfs.orfs.length > 5 && orfsOk, `${tools.orfs.orfs.length} ORFs`);
        assert('tool.seq.kmers: the k-mers of phiX174 in all are L − k + 1 = 5,376 for k = 11, its distinct ones no more, the commonest first', tools.kmers.total === 5386 - 11 + 1 && tools.kmers.distinct <= tools.kmers.total && tools.kmers.top[0][1] >= tools.kmers.top[4][1], `${tools.kmers.total} ${tools.kmers.distinct}`);
        const [e40, e60, e80] = tools.evalue;
        assert('tool.stats.evalue: BLOSUM62 with gaps 11 + L·1 takes BLAST\'s gapped λ 0.267, K 0.041; the E-value falls and the bits rise with the score, E = m n 2^(−bits)',
          e60.gapped && e60.lambda === 0.267 && e60.K === 0.041 && e40.evalue > e60.evalue && e60.evalue > e80.evalue && e40.bits < e60.bits && Math.abs(Math.log(e60.evalue) - Math.log(300 * 1e8 * 2 ** -e60.bits)) < 1e-9, `${e40.evalue} ${e60.evalue} ${e80.evalue}`);
        assert('tool.tree.read: a tree rooted at its midpoint keeps its five leaves and its length 12', tools.tree.leaves === 5 && Math.abs(tools.tree.length - 12) < 1e-9 && tools.tree.tree.nodes.filter((x) => x.parent === -1).length === 1, `${tools.tree.leaves} ${tools.tree.length}`);
        const depth = (2000 * 150) / 5386;
        assert('tool.reads.qc: 2,000 reads of 150 bases read from FASTA, their GC near the genome\'s', tools.qc.reads === 2000 && Math.abs(tools.qc.mean_length - 150) < 0.5 && Math.abs(tools.qc.gc - gc) < 0.02, `${tools.qc.reads} ${tools.qc.gc}`);
        assert('tool.reads.map: the reads mapped back to phiX174, over 95% of them, at a mean depth within 5% of 2,000 · 150 / 5,386', tools.map.mapped / tools.map.reads > 0.95 && Math.abs(tools.map.mean_depth / depth - 1) < 0.05 && tools.map.coverage.length === 1000, `${tools.map.mapped} ${tools.map.mean_depth}`);
      } else skip('phiX174 is not in the data cache: its answers are not checked');
      if (out.ecoli) {
        // oriC of E. coli K-12 MG1655 lies at 3,925,744..3,925,975 (RefSeq's own annotation of NC_000913.3, rep_origin).
        check('Escherichia coli: the minimum of the cumulative GC skew lies within 500 bases of oriC', out.ecoli.origin - 3925860, 500, 'bases');
        // Two replichores of nearly equal length put the terminus about opposite the origin.
        const L = out.ecoli.length, opposite = (out.ecoli.origin + L / 2) % L;
        check('Escherichia coli: its maximum, the terminus, lies within a tenth of the genome of the point opposite the origin', (out.ecoli.terminus - opposite) / L, 0.1);
        // A self-trained gene finder of two periodic Markov chains, scored by the stop codons of
        // RefSeq's CDSs: GeneMarkS reports some 95% on bacteria; this simpler one, above 85%.
        const genes = out.ecoli_genes;
        assert('Escherichia coli: the gene finder, trained on the genome alone, finds over 85% of RefSeq\'s CDSs, and over 95% of its predictions are among them', genes.sensitivity > 0.85 && genes.specificity > 0.95, `${genes.sensitivity} ${genes.specificity}`);
      } else skip('Escherichia coli is not in the data cache: its origin and genes are not checked');
      const sch = out.schedules;
      assert('align: Gotoh, Myers–Miller, the band, the WFA, BiWFA and the stripes reach Gotoh\'s optimum on 1.5 kb relatives', Object.values(sch).every((x) => x.same === true), JSON.stringify(sch));
      const runs = out.compare.runs;
      assert('align: every schedule of the same order agrees with the reference, morphe\'s generic sum among them', runs.filter((x) => x.agrees !== null).length >= 6 && runs.every((x) => x.agrees !== false), runs.map((x) => `${x.name} ${x.agrees}`).join(', '));
      const bt = out.batch;
      assert('align: 3,000 reads scored by scalar rows, stripes and one target to a lane give the same integers', JSON.stringify(bt.gotoh) === JSON.stringify(bt.striped) && JSON.stringify(bt.gotoh) === JSON.stringify(bt.inter) && bt.gotoh.length === 3000);
      const poa = out.poa;
      const mean = poa.read_distances.reduce((x, y) => x + y, 0) / poa.read_distances.length;
      assert('a partial order: ten reads of 150 bases with 8% errors (about 12 edits each) give a consensus within 2 edits of the truth', poa.distance <= 2 && mean > 6, `${poa.distance} from reads ${poa.read_distances}`);
      const a = out.job5, b = out.job13;
      assert('a job: the same answer in 5 shards and in 13, whichever worker computed which', a.errors === b.errors && JSON.stringify(a.quality) === JSON.stringify(b.quality) && a.minimizers === b.minimizers && a.reads === 6000, `${a.errors} ${b.errors}`);
      assert('a job: its trace has every shard, each on a worker, within the job\'s time', a._job.trace.length === 5 && a._job.trace.every((t) => t.end >= t.start && t.end <= a._job.seconds + 1e-6));
      check('a job: the simulated error rate rises along the read (last tenth over first tenth, 6.7 for 0.001 → 0.01 over 150 bases)', b.error_rate.slice(-15).reduce((s, v) => s + v, 0) / b.error_rate.slice(0, 15).reduce((s, v) => s + v, 0) / 6.7 - 1, 0.4);

      // Group grammars (domains/genomics/theory/group-grammars.md), through the engine's grammar.* operations.
      const gr = out.grammar;
      const sc = gr.scan;
      assert('grammar.parse: S → aa | bb under a ↔ b is invariant with one rule orbit; the orbit-only scan accepts ab, which Earley\'s parser and the equivariant parser reject (the draft\'s scan, Proposition 2.1)',
        sc.invariant && sc.orbits.count === 1 && sc.equivariant.scan.accepts && !sc.earley.accepts && !sc.equivariant.consistent.accepts && sc.equivariant.survey.outside > 0,
        `scan ${sc.equivariant.scan.accepts}, Earley ${sc.earley.accepts}, equivariant ${sc.equivariant.consistent.accepts}; ${sc.equivariant.survey.outside} short words outside the language accepted by the scan`);
      const orbitCounts = ['scan', 'brackets', 'hairpin'].map((k) => [gr[k].orbits.count, gr[k].orbits.burnside]);
      assert('grammar.parse: Burnside\'s count equals the rule orbits: 1 for S → aa | bb, 3 for brackets under the mirror (above |R|/|G| = 1.5), 13 for the hairpin under the strand group',
        orbitCounts.every(([o, b]) => o === b) && orbitCounts.map(([o]) => o).join() === '1,3,13' && gr.brackets.orbits.ratio === 1.5, JSON.stringify(orbitCounts));
      const hpq = gr.hairpin.equivariant;
      assert('grammar.parse: on the hairpin, the equivariant parser under {1, κ} accepts the folded word, its chart standing for Earley\'s exactly (Proposition 2.3); they agree on every short word',
        gr.hairpin.earley.accepts && hpq.consistent.accepts && hpq.order === 2 && hpq.chart_agrees && hpq.survey.equivariant_agrees && sc.equivariant.chart_agrees && sc.equivariant.survey.equivariant_agrees,
        `${hpq.consistent.plain_equivalent} and ${gr.hairpin.earley.items} items`);
      const ms = gr.monotone.monotone_search;
      assert('grammar.parse: the monotone w cⁿ dⁿ under a ↔ b derives abccdd; its canonical search visits fewer classes than forms, and at least half as many (Proposition 2.4)',
        !gr.monotone.context_free && gr.monotone.monotone && ms.classes.derives && ms.forms.derives && ms.classes.visited < ms.forms.visited && 2 * ms.classes.visited >= ms.forms.visited, `${ms.classes.visited} classes, ${ms.forms.visited} forms`);
      const c1 = gr.code1, c2 = gr.code2;
      assert('grammar.code: the vertebrate mitochondrial code\'s stabiliser is {1, τ₃}; the standard code\'s is {1}, and τ₃ fails there at ATA/ATG and TGA/TGG alone (Proposition 6.1)',
        c2.stabiliser.map((s) => s.bits).join() === `0,${c2.tau3}` && c2.stabiliser[1].name === 'τ₃' && c1.stabiliser.map((s) => s.bits).join() === '0' && c1.element.keeps === 60 && c1.element.failures.join() === 'ATA/ATG,TGA/TGG',
        `${c2.stabiliser.map((s) => s.name).join(', ')}; ${c1.element.failures.join(' ')}`);
      assert('grammar.code: the power spectrum is the Walsh–Hadamard transform of the profile, and invariance is vanishing on the odd characters, in both codes (Proposition 6.2)', c1.wiener_khinchin && c2.wiener_khinchin && c1.duality && c2.duality);
      const ch = gr.chars;
      const degrees = (x) => x.characters.map((c) => c.degree).sort((p, q) => p - q).join();
      assert('grammar.characters: Σ χ(1)² = |G| and the rows orthonormal for C₄, V₄, D₄, Q₈, S₄ and A₅; each Wedderburn block of dimension χ(1)², measured',
        Object.values(ch).every((x) => x.sum_of_squares === x.order && x.orthonormal && x.wedderburn.measured && x.wedderburn.blocks.every((bl) => bl.dimension === bl.degree ** 2)),
        Object.values(ch).map((x) => `${x.name} ${x.sum_of_squares}/${x.order}`).join(', '));
      assert('grammar.characters: C₄ and C₂ × C₂, D₄ and Q₈ have the same degrees (the algebra does not determine the group); the indicators of D₄\'s and Q₈\'s two-dimensional characters are +1 and −1; S₄\'s degrees are the hook lengths\'',
        degrees(ch['cyclic 4']) === degrees(ch.klein) && degrees(ch['dihedral 4']) === degrees(ch.quaternion) && degrees(ch.quaternion) === '1,1,1,1,2'
          && ch['dihedral 4'].characters.find((c) => c.degree === 2).indicator === 1 && ch.quaternion.characters.find((c) => c.degree === 2).indicator === -1 && ch['symmetric 4'].symmetric.degrees_agree && degrees(ch['alternating 5']) === '1,3,3,4,5',
        `${degrees(ch['dihedral 4'])} and ${degrees(ch.quaternion)}; A5 ${degrees(ch['alternating 5'])}`);
      const nu = gr.nucleotides;
      const byName = Object.fromEntries(nu.classifications.map((c) => [c.name, c]));
      assert('grammar.nucleotides: the cosets of V₄\'s subgroups are the bases, the six two-letter codes and N, and B, D, H, V none; at the base point T the purine character changes sign, the weak–strong one does not',
        nu.coset_count === 11 && !nu.three_letter_codes_are_cosets && nu.product_rule && byName['purine–pyrimidine'].on_bases.join() === '-1,1,-1,1' && byName['weak–strong'].on_bases.join() === '1,-1,-1,1' && byName['purine–pyrimidine'].letters.join('/') === 'Y/R',
        nu.classifications.map((c) => `${c.letters.join('/')} ${c.on_bases.join(' ')}`).join('; '));
      if (gr.hairpins && out.phix_letters) {
        // The definition, again, here: every span of at most 80 of phiX174's first 300 bases that is a
        // stem of four Watson–Crick pairs or more closing a loop of 3 to 8.
        const wc = new Set(['AT', 'TA', 'CG', 'GC']);
        const hairpin = (w) => { for (let st = 4; 2 * st + 3 <= w.length; st++) { const loop = w.length - 2 * st; if (loop > 8) continue; let ok = true; for (let i = 0; i < st && ok; i++) ok = wc.has(w[i] + w[w.length - 1 - i]); if (ok) return true; } return false; };
        const text = out.phix_letters.slice(0, 300);
        let count = 0;
        for (let len = 1; len <= 80; len++) for (let i = 0; i + len <= 300; i++) count += hairpin(text.slice(i, i + len));
        const h = gr.hairpins;
        assert('grammar.hairpins: on phiX174\'s first 300 bases the inside pass finds every hairpin span, as many as counted here by the definition; graduated breaking over R/Y, W/S and M/K keeps them all',
          h.agrees_with_definition && h.spans === count && h.by_definition === count && h.coarse.exact && h.coarse.by.length === 3, `${h.spans} spans, ${count} by the definition, ${h.maximal_count} maximal`);
      } else skip('phiX174 is not in the data cache: its hairpins are not checked');
      if (gr.symmetry) {
        const sy = gr.symmetry;
        const rcTV = sy.elements.find((e) => e.short === 'κρ').total_variation, rest = Math.min(...sy.elements.filter((e) => !['1', 'κρ'].includes(e.short)).map((e) => e.total_variation));
        assert('grammar.symmetry: Escherichia coli\'s 4-mer statistics keep the reverse complement alone among the eight symmetries (Chargaff\'s second rule); the eight channels sum to ‖f‖²',
          rcTV < 0.05 && rest > 4 * rcTV && Math.abs(sy.channel_sum - sy.norm) / sy.norm < 1e-12 && sy.channels.length === 8 && sy.odd_rc < 1e-3, `κρ ${rcTV.toFixed(4)}, others from ${rest.toFixed(3)}`);
      } else skip('Escherichia coli is not in the data cache: the symmetry of its k-mers is not checked');
    }
  }
}

done();
