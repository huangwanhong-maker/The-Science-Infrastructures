// What this workbench draws on, and what it does not: each tool and dataset on the list it was asked to integrate, what was
// taken from it, and, of those not taken, why. The list stays beside the work, so that a reader can judge what stands behind a picture.

const TAKEN = [
  ['Overture Maps (2026-10-21.0)', 'Buildings, ways, water and open land of any place, read by range requests from its world archives as the camera arrives.', 'ODbL 1.0; footprints from OpenStreetMap, Google Open Buildings, Microsoft and others'],
  ['GHS-BUILT-H R2023A', 'The estimated height of a building that records none: the average of the buildings on its 100 m cell in 2018. Drawn paler, and said to be an estimate.', 'CC BY 4.0, European Commission JRC'],
  ['Space2Stats', 'The H3 level 6 grid of people, lights, built land, urbanisation and hazards, for any region.', 'World Bank; each field keeps its source\'s terms'],
  ['GHS-SMOD by country, R2023A', 'The Degree of Urbanisation on the ground (from the global grid), and in Space2Stats\' cells.', 'CC BY 4.0, European Commission JRC'],
  ['World Settlement Footprint Tracker', 'Its data, the WSF Evolution (1985–2015), on the ground with a scale of years. The Tracker itself offers no open grids.', 'CC BY 4.0, DLR'],
  ['blackmarblepy', 'What it reads: Black Marble, on the globe through GIBS, and its sums through Space2Stats. The package itself needs a NASA token, and is not used.', 'NASA open data'],
  ['ThinkHazard! 3.0', 'Its public reports, each hazard\'s level for an administrative unit.', 'GFDRR, by its open API; its code (GPL) is not copied'],
  ['Risk Data Library Standard 1.0', 'The risk layers\' records are described in its fields; the schema is fetched, not kept here.', 'CC BY-SA 4.0'],
  ['Freetown Flood Risk Assessment', 'Depths at 1 m by return period, the exposed portfolio and the event and exceedance results, on the map.', 'World Bank, CC BY 4.0'],
  ['Poverty and Inequality Platform', 'Every economy\'s latest survey, over the globe.', 'CC BY 4.0, World Bank'],
  ['data360-mcp', 'What stands behind it, Data360\'s API, in the tool “Any indicator on the globe”. The MCP server is for language models, and is not run here.', 'World Bank Data360'],
];
const METHODS = [
  ['giga-spatial', 'Binning any quantity on H3 cells: written independently, on a port of H3 to physis held to h3-py.', 'AGPL-3.0; its code is not copied'],
  ['GOSTnets', 'Shortest paths and accessibility over a way network: written independently, in the access analyses.', 'its method, not its code'],
  ['eca-resilience', 'The space usage index, on Melbourne\'s pedestrian counts.', 'its method, not its code'],
  ['GeoE3', 'A composite of factors scaled to a common range: written independently, in nomos.', 'GPL-3.0; its code is not copied'],
  ['RAPIDA', 'Exposure as hazard over people and buildings: written independently.', 'its method, not its code'],
  ['CCDR-tools', 'Expected annual loss from an exceedance curve: written independently.', 'GPL-3.0; its code is not copied'],
];
const LEFT = [
  ['montevideo-brt', 'A study of one bus corridor, on partner data that is not open.'],
  ['GRAS (procurement)', 'Screens procurement documents; it holds nothing spatial.'],
  ['Metadata Editor', 'A cataloguer\'s application; the workbench needs no second catalogue beside the data library\'s manifests.'],
  ['Proof-Carrying Numbers', 'A protocol for checking numbers a language model reports; no language model stands between these numbers and the page.'],
  ['DevDataBench', 'A benchmark of language models on development data.'],
  ['REaLTabFormer, and the guide to synthetic data', 'They make synthetic microdata; this workbench draws only what was recorded, or says it is a model.'],
  ['MANAGE-WB', 'Written in GAMS, a proprietary language, and for one model of a national economy.'],
  ['Mobile Phone Data Initiatives manual', 'A manual for working with operators\' records, which are not open; nothing here could be drawn from it.'],
  ['ai4data', 'Language-model tools for discovering datasets; the workbench\'s datasets are named.'],
];

export async function render(root, ctx) {
  const { dom: { h } } = ctx;
  const table = (rows, heads) => h('table', { class: 'sources' }, h('thead', {}, h('tr', {}, heads.map((x) => h('th', {}, x)))), h('tbody', {}, rows.map((r) => h('tr', {}, r.map((x, i) => h(i === 0 ? 'th' : 'td', {}, x))))));
  root.append(h('h1', {}, 'Sources and what stands behind them'),
    h('p', { class: 'muted' }, 'What this workbench takes from the tools and data it was asked to draw on, and what it leaves. Whatever is computed here is computed by the foundations and domains of the infrastructure, in C++, and checked there.'),
    h('h2', {}, 'Data on the map'), table(TAKEN, ['Source', 'What is taken', 'Terms']),
    h('h2', {}, 'Methods written here'), table(METHODS, ['Tool', 'What was learnt from it', 'Note']),
    h('h2', {}, 'Left aside, and why'), table(LEFT, ['Tool', 'Why']));
}
