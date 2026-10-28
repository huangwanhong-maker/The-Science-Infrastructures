// The Genome Workbench: what it offers, declared (core/registry.js), and the shell started
// (shell.js). A new workspace is a view declared here; a new instrument is a tool declared in
// src/tools/ and imported below; a new action of the command palette is a command.

import { view, command } from './core/registry.js';
import { start } from './shell.js';
import * as api from './api.js';
import { toast } from './ui/kit.js';
import * as genome from './views/genome.js';
import * as align from './views/align.js';
import * as indexView from './views/index.js';
import * as reads from './views/reads.js';
import * as assembly from './views/assembly.js';
import * as expression from './views/expression.js';
import * as trees from './views/trees.js';
import * as populations from './views/populations.js';
import * as grammar from './views/grammar.js';
import * as toolbox from './views/toolbox.js';
import * as data from './views/data.js';
import * as run from './views/run.js';
import './tools/index.js';

view({ id: 'genome', title: 'Genome', tagline: 'its map, its features, and what the sequence alone says', icon: 'genome', group: 'explore', module: genome, full: true, cls: 'view-genome' });
view({ id: 'align', title: 'Align', tagline: 'two sequences, every algorithm, every semiring', icon: 'align', group: 'analyse', module: align });
view({ id: 'index', title: 'Index', tagline: 'suffix arrays, the Burrows–Wheeler transform, the FM-index', icon: 'index', group: 'analyse', module: indexView });
view({ id: 'reads', title: 'Reads', tagline: 'a diploid read, mapped and called against its truth', icon: 'reads', group: 'analyse', module: reads });
view({ id: 'assembly', title: 'Assembly', tagline: 'a genome rebuilt from its reads', icon: 'assembly', group: 'analyse', module: assembly });
view({ id: 'expression', title: 'Expression', tagline: 'RNA-seq quantified, and two conditions compared', icon: 'expression', group: 'analyse', module: expression });
view({ id: 'trees', title: 'Trees', tagline: 'genes and genomes on the tree of their descent', icon: 'trees', group: 'analyse', module: trees });
view({ id: 'populations', title: 'Populations', tagline: 'genealogies, drift and structure', icon: 'populations', group: 'analyse', module: populations });
view({ id: 'grammar', title: 'Grammar', tagline: 'genomes under a group: the bases, parsing up to symmetry, strands, the code', icon: 'grammar', group: 'analyse', module: grammar });
view({ id: 'toolbox', title: 'Toolbox', tagline: 'every instrument, on your own sequences and files', icon: 'toolbox', group: 'make', module: toolbox, full: true });
view({ id: 'data', title: 'Your data', tagline: 'genomes and files, yours and the library\'s', icon: 'data', group: 'make', module: data, needsEngine: false });
view({ id: 'run', title: 'Engine', tagline: 'the machine, the workers, the jobs', icon: 'run', group: 'system', module: run, needsEngine: false });

command({ id: 'theme', title: 'Switch between day and night', icon: 'moon', keywords: ['theme', 'dark', 'light'], run: () => document.getElementById('theme-btn')?.click() });
command({ id: 'add-file', title: 'Add files to your data…', icon: 'upload', keywords: ['upload', 'open', 'import', 'fasta', 'fastq', 'vcf'], run: (ctx) => ctx.go('data', { add: 1 }) });
command({ id: 'restart', title: 'Restart the engine\'s workers', icon: 'refresh', keywords: ['engine', 'workers'], run: async () => { await api.restart(); toast('The workers will start again with the next request.', 'good'); } });
command({ id: 'rail', title: 'Narrow or widen the rail', icon: 'sidebar', run: () => document.querySelector('.rail-foot .icon-btn')?.click() });

start().catch((e) => toast(`The workbench could not start: ${e.message}`, 'error', 20000));
