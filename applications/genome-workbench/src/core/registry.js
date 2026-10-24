// What the workbench offers, each part declared once and found from here by the rail, the
// toolbox, the command palette and the data library:
//
//   view({id, title, tagline, about, group, icon, module})     a workspace in the rail
//   tool({id, title, group, icon, summary, about, refs, inputs, run, render, accepts})
//                                                              an instrument of the toolbox
//   command({id, title, icon, keywords, run})                  an action of the palette
//
// A new algorithm joins by declaring a tool (src/tools/); a new workspace by declaring a view. The
// tools' inputs are a schema the toolbox turns into a form (ui/form.js); their results are drawn by
// the renderers of ui/results.js or by the tool's own render.

export const RAIL_GROUPS = [
  ['explore', 'Explore'],
  ['analyse', 'Analyse'],
  ['make', 'Instruments'],
  ['system', 'System'],
];

export const TOOL_GROUPS = [
  ['sequence', 'Sequences', 'sequence'],
  ['alignment', 'Alignment', 'align'],
  ['annotation', 'Annotation', 'annotation'],
  ['index', 'Indexes and search', 'index'],
  ['reads', 'Reads and variants', 'reads'],
  ['assembly', 'Assembly', 'assembly'],
  ['expression', 'Expression', 'expression'],
  ['phylogeny', 'Trees and selection', 'trees'],
  ['population', 'Populations', 'populations'],
  ['statistics', 'Statistics', 'statistics'],
];

export const views = [];
export const tools = [];
export const commands = [];

export function view(def) { views.push(def); return def; }
export function tool(def) {
  if (tools.some((t) => t.id === def.id)) throw new Error(`a tool ${def.id} is declared twice`);
  tools.push(def);
  return def;
}
export function command(def) { commands.push(def); return def; }

export const viewById = (id) => views.find((v) => v.id === id);
export const toolById = (id) => tools.find((t) => t.id === id);
// The tools that take a file of this kind (fasta, fastq, vcf …) as one of their inputs.
export const toolsFor = (kind) => tools.filter((t) => (t.inputs || []).some((i) => (i.kinds || []).includes(kind) || (i.type === 'genome' && ['fasta', 'genbank'].includes(kind))));
