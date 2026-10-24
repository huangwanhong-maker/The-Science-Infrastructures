// Checks of what is the Genome Atlas's own: how it writes lengths, positions and ticks; and,
// when the engine is built, its answers held to the protocol of the README (§3): a genome's
// summary and tracks (phiX174's length, GC and genes; the origin of Escherichia coli from the
// cumulative skew within a few hundred bases of oriC), the search of a pattern on both strands
// against a count by brute force, a feature's protein as RefSeq translates it, the gene finder
// against RefSeq's annotation, every schedule of an alignment held to Gotoh's optimum, the batch
// kernels to one another, the consensus of a partial order to its truth, and a job whose answer
// is the same however it is cut into shards and spread over workers. The science is checked
// with the genomics and evolution domains.
// Run from the repository root:  node applications/genome-atlas/tests/run-tests.mjs

import { spawnSync } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';

import { suite, ROOT, CACHE } from '../../../development/tools/checks.mjs';
import * as fmt from '../src/format.js';
import { readHash, writeHash } from '../src/state.js';

const { assert, check, skip, done } = suite();
const APP = resolve(ROOT, 'applications/genome-atlas');

// ---- writing --------------------------------------------------------------------------------------
{
  assert('lengths: 5,386 bp, 48.5 kb, 4.64 Mb', fmt.length(5386) === '5,386 bp' && fmt.length(48502) === '48.5 kb' && fmt.length(4641652) === '4.64 Mb', `${fmt.length(5386)} ${fmt.length(48502)} ${fmt.length(4641652)}`);
  assert('positions are written 1-based, intervals as b..e', fmt.position(0) === '1' && fmt.interval(99, 250) === '100..250');
  assert('ticks step by 1, 2 or 5 times a power of ten', [fmt.tickStep(4641652, 12), fmt.tickStep(16569, 8), fmt.tickStep(5386, 14)].every((s) => /^(1|2|5)0*$/.test(String(s))), `${fmt.tickStep(4641652, 12)} ${fmt.tickStep(16569, 8)}`);
  assert('a megabase genome is labelled in Mb', fmt.tickLabel(3800000, 200000, 4641652) === '3.8 Mb' && fmt.tickLabel(2000, 1000) === '2 kb', fmt.tickLabel(3800000, 200000, 4641652));
  const s = { view: 'align', genome: 'ncbi:NC_000913.3', params: { begin: 100, end: 900, a: 'ACGT' } };
  const back = readHash(writeHash(s));
  assert('the address holds the view, the genome and its parameters', back.view === 'align' && back.genome === s.genome && back.params.begin === '100' && back.params.a === 'ACGT', writeHash(s));
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
