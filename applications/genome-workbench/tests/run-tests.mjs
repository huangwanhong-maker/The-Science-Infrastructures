// Checks of what is the Genome Workbench's own: how it writes lengths, positions and ticks; and,
// when the engine is built, its answers held to the protocol of the README (§3): a genome's
// summary and tracks (phiX174's length, GC and genes; the origin of Escherichia coli from the
// cumulative skew within a few hundred bases of oriC), the search of a pattern on both strands
// against a count by brute force, a feature's protein as RefSeq translates it, the gene finder
// against RefSeq's annotation, every schedule of an alignment held to Gotoh's optimum, the batch
// kernels to one another, the consensus of a partial order to its truth, and a job whose answer
// is the same however it is cut into shards and spread over workers. The science is checked
// with the genomics and evolution domains.
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
    }
  }
}

done();
