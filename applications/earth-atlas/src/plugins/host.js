// The host of plug-ins: loads a script, resolves what it needs from the
// library, runs it with the atlas's context, and hands its results on: field
// products to the library and the catalogue (through the application), sketches
// to the sketch layer, tables and notes to the panel. The application gives the
// host its hooks; the host never touches the page itself.

import { Context } from './context.js';
import * as api from './index.js';
import { modelRank, levelOf } from '../layers/catalogue.js';

const defaults = (parameters) => Object.fromEntries(Object.entries(parameters || {}).map(([k, p]) => [k, p.default]));

export class PluginHost {
  // hooks: { publish(product, show), sketches(list), changed(), progress(), time() -> ms, pinned() -> [lon, lat] | null }
  constructor(library, hooks) {
    this.library = library;
    this.hooks = hooks;
    this.available = [];      // the built-in plug-ins, from plugins/index.json
    this.loaded = new Map();  // id -> entry
    this.error = null;        // the last failure to load
  }

  // The built-in plug-ins offered by the page.
  async index(url = 'plugins/index.json') {
    try { this.available = (await (await fetch(url)).json()).plugins || []; } catch { this.available = []; }
    this.hooks.changed();
    return this.available;
  }

  // Loads a script by URL (one of plugins/, or any other) and asks it for its plug-in.
  async load(url) {
    this.error = null;
    try {
      const module = await import(new URL(url, document.baseURI).href);
      if (typeof module.create !== 'function') throw new Error(`${url} does not export create(api)`);
      const plugin = await module.create(api);
      if (!plugin || !plugin.id || typeof plugin.run !== 'function') throw new Error(`${url}: create(api) did not return a plug-in`);
      const needs = this.resolve(plugin), chosen = {};
      for (const need of needs) {
        const taken = new Set(Object.values(chosen));
        const pick = need.distinct ? need.candidates.find((c) => !taken.has(c.key)) || need.candidates[0] : need.candidates[0];
        chosen[need.name] = pick ? pick.key : null;
      }
      const entry = { id: plugin.id, url, plugin, params: defaults(plugin.parameters), needs, chosen, running: false, progress: 0, text: '', results: [], error: null };
      this.loaded.set(plugin.id, entry);
      return entry;
    } catch (error) {
      console.error(error);
      this.error = `${url}: ${error.message || error}`;
      return null;
    } finally {
      this.hooks.changed();
    }
  }

  unload(id) {
    if (!this.loaded.delete(id)) return;
    this.hooks.sketches(this.sketches());
    this.hooks.changed();
  }

  // The products that can serve each need: one candidate per product and level
  // carrying every standard name asked, the plug-in's own output excepted.
  resolve(plugin) {
    return plugin.needs().map((need) => {
      const standards = need.standard_names || [need.standard_name];
      const candidates = [];
      for (const [name, meta] of this.library.fieldMeta) {
        if (name === plugin.productName || !meta.variables) continue;
        if (need.times && !(meta.times && meta.times.length >= need.times)) continue;
        // Levels as the catalogue reads them, stated or inferred from the name ('2 m' for GFS's tmp2m).
        const plain = meta.variables.filter((v) => !/realization/.test(v.cell_methods || ''));
        const levels = need.level ? [need.level] : [...new Set(plain.map((v) => levelOf(v)))];
        for (const level of levels) {
          const variables = standards.map((std) => plain.find((v) => v.standard_name === std && levelOf(v) === (level || '')));
          if (!variables.every(Boolean)) continue;
          const model = meta.model || meta.title;
          candidates.push({ key: `${name}|${level || ''}`, product: name, level: level || null, model: level && !need.level ? `${model}, ${level}` : model,
            rank: modelRank(model), names: variables.map((v) => v.name), variables, meta });
        }
      }
      // The models in the panel's order, and within a model the levels from the ground up.
      candidates.sort((a, b) => a.rank - b.rank || a.model.localeCompare(b.model, undefined, { numeric: true }));
      return { ...need, candidates };
    });
  }

  choose(id, needName, key) {
    const entry = this.loaded.get(id);
    if (entry) { entry.chosen[needName] = key; this.hooks.changed(); }
  }

  setParameter(id, name, value) {
    const entry = this.loaded.get(id);
    if (entry) entry.params[name] = value;
  }

  async run(id) {
    const entry = this.loaded.get(id);
    if (!entry || entry.running) return;
    entry.running = true; entry.error = null; entry.progress = 0; entry.text = '';
    this.hooks.changed();
    try {
      const inputs = {};
      for (const need of entry.needs) {
        const c = need.candidates.find((x) => x.key === entry.chosen[need.name]);
        if (!c) throw new Error(`nothing in the library serves "${need.label || need.name}"`);
        inputs[need.name] = { product: await this.library.field(c.product), names: c.names, variables: c.variables, meta: c.meta, model: c.model, level: c.level };
      }
      const context = new Context(this.library, {
        time: this.hooks.time(), pinned: this.hooks.pinned(),
        progress: (fraction, text) => { entry.progress = fraction; entry.text = text || ''; this.hooks.progress(); },
      });
      entry.results = (await entry.plugin.run(inputs, { ...entry.params }, context)) || [];
      for (const r of entry.results) if (r.kind === 'field') await this.hooks.publish(r.product, r.show || {});
      this.hooks.sketches(this.sketches());
    } catch (error) {
      console.error(error);
      entry.error = `The plug-in could not run: ${error.message || error}`;
    }
    entry.running = false;
    this.hooks.changed();
  }

  // Every loaded plug-in's sketches, for the sketch layer.
  sketches() { return [...this.loaded.values()].flatMap((e) => e.results.filter((r) => r.kind === 'sketch')); }

  // What the panel shows.
  view() {
    return {
      available: this.available,
      error: this.error,
      loaded: [...this.loaded.values()].map((e) => ({
        id: e.id, url: e.url, title: e.plugin.title, description: e.plugin.description,
        parameters: e.plugin.parameters || {}, params: e.params,
        needs: e.needs.map((n) => ({ name: n.name, label: n.label || n.name, chosen: e.chosen[n.name], candidates: n.candidates.map((c) => ({ key: c.key, model: c.model })) })),
        running: e.running, progress: e.progress, text: e.text, error: e.error,
        results: e.results.filter((r) => r.kind === 'table' || r.kind === 'note'),
        sketches: e.results.filter((r) => r.kind === 'sketch').map((s) => s.title),
      })),
    };
  }
}
