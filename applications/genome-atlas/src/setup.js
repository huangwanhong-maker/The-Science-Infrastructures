// The recipe from the page (recipe.json through the server's /setup): which steps are done, a
// button that runs them, and the log while they run. Shown where a view needs the engine or the
// genomes and they are not there yet.

import { el, clear } from './dom.js';
import * as api from './api.js';

export function setupPanel(onDone) {
  const box = el('div', { class: 'setup' });
  let timer = null;
  const render = (v) => {
    clear(box);
    box.append(el('h2', { class: 'view-title' }, 'Setup'));
    box.append(el('p', { class: 'lede' }, 'The atlas computes with its engine, the genomics and evolution domains in C++ compiled for this processor, and shows genomes fetched from NCBI. Both are made once by the recipe: about two minutes and 12 MB.'));
    const list = el('ol');
    for (const s of v.steps || []) list.append(el('li', {}, el('span', { class: `badge ${s.done ? 'good' : ''}` }, s.done ? 'done' : 'to do'), ' ', s.title));
    box.append(list);
    const button = el('button', { class: 'button primary', disabled: v.running || !v.writable, onclick: async () => {
      try { render(await api.runSetup(v.steps.every((s) => s.done))); poll(); } catch (e) { box.append(el('p', { class: 'error' }, e.message)); }
    } }, v.running ? 'Running…' : v.steps?.every((s) => s.done) ? 'Run again' : 'Build the engine and fetch the genomes');
    box.append(el('div', { class: 'row' }, button, v.writable ? null : el('span', { class: 'faint small' }, v.refusal || '')));
    if (v.log) box.append(el('pre', {}, v.log));
    if (v.state === 'failed') box.append(el('p', { class: 'error' }, `The recipe failed at: ${(v.failed || []).join(', ')}.`));
  };
  const poll = () => {
    clearTimeout(timer);
    timer = setTimeout(async () => {
      const v = await api.setup().catch(() => null);
      if (!v) return;
      render(v);
      if (v.running) poll();
      else if (v.state === 'done' && onDone) onDone();
    }, 1200);
  };
  api.setup().then((v) => { render(v); if (v.running) poll(); }).catch((e) => box.append(el('p', { class: 'error' }, e.message)));
  return box;
}
