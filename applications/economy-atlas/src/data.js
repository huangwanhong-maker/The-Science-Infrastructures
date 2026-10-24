// One's own data (ADR 0046): flows between places, or currencies held and exchanged, given as CSV or
// JSON, kept by the server as products of the library's own schemas (datasets.py), and drawn by every
// view and given to every plug-in as the library's are.

import * as api from './api.js';
import { el, clear, put, fill } from './dom.js';
import { shared, setDatasets, emit, datasetsOf } from './state.js';

const FORMATS = {
  flows: {
    title: 'Flows between places',
    help: 'One file of flows, and, for places that are not economies of the library, one of places. Ids that are ISO alpha-3 codes (USA, DEU) or the names of economies need no place.',
    example: 'flows.csv\nfrom,to,value,year\nFRA,DEU,104.2,2024\nDEU,FRA,121.7,2024\n\nplaces.csv (optional)\nid,name,latitude,longitude,group\nlyon,Lyon,45.76,4.84,France',
  },
  currencies: {
    title: 'Currencies held and exchanged',
    help: 'Monthly rates against a base currency, and each currency\'s policy rate in per cent a year. The base\'s rows carry its own policy rate, with a rate of 1 or none.',
    example: 'rates.csv\nmonth,currency,rate,policy\n2024-01,USD,,5.33\n2024-01,JPY,146.5,-0.10\n2024-01,GBP,0.787,5.25',
  },
};

const dialogState = { kind: 'flows', files: [] };

async function refresh() {
  const { datasets } = await api.datasets();
  setDatasets(datasets);
}

function render() {
  const dialog = document.getElementById('data-dialog');
  const F = FORMATS[dialogState.kind];
  const input = el('input', { type: 'file', multiple: true, accept: '.csv,.tsv,.txt,.json', hidden: true,
    onchange: (e) => { dialogState.files = [...e.target.files]; render(); } });
  const drop = el('div', { class: 'drop', tabindex: 0, role: 'button', 'aria-label': 'Choose files',
    onclick: () => input.click(), onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') input.click(); },
    ondragover: (e) => { e.preventDefault(); drop.classList.add('over'); }, ondragleave: () => drop.classList.remove('over'),
    ondrop: (e) => { e.preventDefault(); dialogState.files = [...e.dataTransfer.files]; render(); } },
  el('b', { text: 'Drop CSV or JSON files here' }), el('span', { class: 'note', text: 'or choose them' }), input);
  const title = el('input', { type: 'text', id: 'data-name', placeholder: 'a title for these data', value: dialogState.title || '',
    oninput: (e) => { dialogState.title = e.target.value; } });
  const base = el('input', { type: 'text', id: 'data-base', placeholder: 'e.g. USD (found by itself if its rate is 1)', value: dialogState.base || '',
    oninput: (e) => { dialogState.base = e.target.value; } });
  const status = el('span', { class: 'note', id: 'data-status' });
  const own = datasetsOf(dialogState.kind).filter((d) => !d.library);
  fill(dialog, 
    el('div', { class: 'dialog-head' }, el('h2', { id: 'data-title', text: 'Your data' }),
      el('button', { class: 'quiet', onclick: () => dialog.close(), text: 'close' })),
    el('div', { class: 'dialog-body' },
      el('div', { class: 'kinds' }, Object.entries(FORMATS).map(([k, f]) => el('button', { 'aria-pressed': String(k === dialogState.kind),
        onclick: () => { dialogState.kind = k; dialogState.files = []; render(); } }, f.title))),
      el('p', { class: 'note', text: F.help }),
      drop,
      dialogState.files.length ? el('div', { class: 'files' }, dialogState.files.map((f) => el('span', { text: `${f.name} · ${(f.size / 1024).toFixed(1)} kB` }))) : null,
      el('label', { class: 'field' }, el('span', { text: 'Title' }), title),
      dialogState.kind === 'currencies' ? el('label', { class: 'field' }, el('span', { text: 'Base currency' }), base) : null,
      el('details', { class: 'format' }, el('summary', { text: 'The format, with an example' }), el('pre', { text: F.example })),
      own.length ? el('div', {}, el('div', { class: 'heading' }, el('span', { text: 'kept here' })),
        own.map((d) => el('div', { class: 'dataset' }, el('b', { text: d.title }),
          el('small', { text: describe(d) }),
          el('div', {}, el('button', { class: 'quiet', text: 'forget', onclick: async () => { await api.removeDataset(d.id); await refresh(); render(); } }))))) : null),
    el('div', { class: 'dialog-foot' }, status,
      el('button', { class: 'primary', text: 'Add these data', disabled: !dialogState.files.length || null, onclick: () => add(status) })),
  );
}

async function add(status) {
  status.textContent = 'reading…';
  status.classList.remove('error');
  try {
    const files = await Promise.all(dialogState.files.map(async (f) => ({ name: f.name, text: await f.text() })));
    const about = await api.addDataset(dialogState.kind, dialogState.title || '', files, dialogState.base || '');
    await refresh();
    shared[dialogState.kind === 'flows' ? 'flows' : 'currencies'] = about.id;
    emit('dataset-chosen', dialogState.kind);
    dialogState.files = [];
    dialogState.title = '';
    document.getElementById('data-dialog').close();
  } catch (error) {
    status.textContent = error.message;
    status.classList.add('error');
  }
}

export function describe(d) {
  if (d.kind === 'flows') {
    const years = d.years || [];
    const span = years.length > 1 ? `${years[0]}–${years.at(-1)}` : years[0] ? String(years[0]) : 'no year';
    return `${d.places} places · ${(d.flows ?? 0).toLocaleString('en-GB')} flows · ${span}${d.unit ? ` · ${d.unit}` : ''}`;
  }
  return `${(d.currencies || []).length} currencies on ${d.base} · ${(d.months || []).join(' to ')}`;
}

export function openDataDialog(kind = 'flows') {
  dialogState.kind = kind;
  render();
  document.getElementById('data-dialog').showModal();
}

// A choice of dataset for a view, with the means to add one's own.
export function datasetChooser(kind, onChange) {
  const key = kind === 'flows' ? 'flows' : 'currencies';
  const select = el('select', { 'aria-label': 'Dataset', onchange: (e) => { shared[key] = e.target.value; onChange(); } },
    datasetsOf(kind).map((d) => el('option', { value: d.id, text: d.library ? d.title : `${d.title} (yours)`, selected: d.id === shared[key] })));
  const chosen = datasetsOf(kind).find((d) => d.id === shared[key]);
  return el('div', {},
    el('label', { class: 'field' }, el('span', { text: kind === 'flows' ? 'Flows' : 'Record' }), select),
    chosen ? el('div', { class: 'note', text: describe(chosen) }) : null,
    el('div', { style: { marginTop: '8px' } }, el('button', { class: 'quiet', text: '+ your data', onclick: () => openDataDialog(kind) })));
}
