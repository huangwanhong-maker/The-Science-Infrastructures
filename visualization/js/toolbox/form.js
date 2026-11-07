// A tool's form, built from its declaration (catalogue.js): the columns its roles take, and its parameters. It holds no state of
// its own: `form` is { table, paste, where, roles, params }, and a control that changes what can be chosen asks for a redraw.
// Styled by each application's stylesheet (.field, .chips, .chip, .toggles).

import { el } from './dom.js';

export function freshForm(t) {
  const p = t.preset || {};
  const params = {};
  for (const q of t.params || []) params[q.name] = q.default;
  return {
    table: t.noTable ? null : p.paste ? '__paste' : p.table || 'economies',
    paste: p.paste || '',
    where: { ...(p.where || {}) },
    roles: JSON.parse(JSON.stringify(p.roles || {})),
    params: { ...params, ...(p.params || {}) },
  };
}

// A table typed or pasted: its columns, numbers where every filled cell is one. Separated by tabs, semicolons or commas.
export function parsePaste(text) {
  const lines = text.trim().split(/\r?\n/).filter((l) => l.trim());
  if (lines.length < 2) return { columns: [], data: {} };
  const sep = lines[0].includes('\t') ? '\t' : lines[0].includes(';') ? ';' : ',';
  const names = lines[0].split(sep).map((s) => s.trim());
  const data = Object.fromEntries(names.map((n) => [n, []]));
  for (const line of lines.slice(1)) {
    const cells = line.split(sep);
    names.forEach((n, i) => { const c = (cells[i] ?? '').trim(); data[n].push(c === '' ? null : c); });
  }
  const columns = names.map((n) => {
    const numeric = data[n].every((v) => v === null || Number.isFinite(Number(v)));
    if (numeric) data[n] = data[n].map((v) => (v === null ? null : Number(v)));
    return { name: n, kind: numeric ? 'number' : 'text' };
  });
  return { columns, data };
}

export function field(label, control, help) {
  return el('label', { class: 'field', title: help || '' }, el('span', { text: label }), control);
}

export function roleControl(role, form, columns, redraw) {
  const fits = columns.filter((c) => role.kind !== 'number' || c.kind === 'number');
  const value = form.roles[role.role];
  if (role.many) {
    const chosen = new Set(Array.isArray(value) ? value : value ? [value] : []);
    return el('div', { class: 'field' }, el('span', { text: `${role.label}${role.optional ? ' (optional)' : ''}` }),
      el('div', { class: 'chips' }, fits.map((c) => el('button', { class: `chip${chosen.has(c.name) ? ' on' : ''}`, type: 'button', title: c.label || '',
        onclick: () => { chosen.has(c.name) ? chosen.delete(c.name) : chosen.add(c.name); form.roles[role.role] = [...chosen]; redraw(); } }, c.name))));
  }
  return field(`${role.label}${role.optional ? ' (optional)' : ''}`, el('select', { onchange: (e) => { form.roles[role.role] = e.target.value || undefined; } },
    el('option', { value: '', text: role.optional ? '—' : 'choose a column' }),
    fits.map((c) => el('option', { value: c.name, text: c.label ? `${c.name} · ${c.label}` : c.name, selected: c.name === value }))));
}

export function paramControl(q, form) {
  const v = form.params[q.name];
  const set = (x) => { form.params[q.name] = x; };
  if (q.kind === 'choice') {
    return field(q.label, el('select', { onchange: (e) => set(e.target.value) },
      q.options.map(([k, label]) => el('option', { value: k, text: label, selected: k === v }))), q.help);
  }
  // One of the columns a role has taken.
  if (q.kind === 'chosen') {
    const among = [form.roles[q.role]].flat().filter(Boolean);
    return field(q.label, el('select', { onchange: (e) => set(e.target.value) },
      among.map((k) => el('option', { value: k, text: k, selected: k === v }))), q.help);
  }
  if (q.kind === 'boolean') {
    return el('label', { class: 'toggles' }, el('label', {}, el('input', { type: 'checkbox', checked: Boolean(v), onchange: (e) => set(e.target.checked) }), q.label));
  }
  if (q.kind === 'area') {
    return field(q.label, el('textarea', { rows: 5, class: 'mono', oninput: (e) => set(e.target.value) }, v ?? ''), q.help);
  }
  if (q.kind === 'text') return field(q.label, el('input', { type: 'text', value: v ?? '', oninput: (e) => set(e.target.value) }), q.help);
  return field(q.label, el('input', { type: 'number', value: v ?? '', step: q.step ?? (q.kind === 'integer' ? 1 : 'any'), min: q.min ?? undefined, max: q.max ?? undefined,
    oninput: (e) => set(e.target.value === '' ? null : q.kind === 'integer' ? parseInt(e.target.value, 10) : Number(e.target.value)) }), q.help);
}

// What a request takes of a form: the roles filled, and whether something required is not. Answers { missing } or { roles }.
export function rolesOf(t, form) {
  for (const r of t.roles || []) {
    const v = form.roles[r.role];
    if (!r.optional && (!v || (Array.isArray(v) && !v.length))) return { missing: r.label.toLowerCase() };
  }
  return { roles: Object.fromEntries(Object.entries(form.roles).filter(([, v]) => v && (!Array.isArray(v) || v.length))) };
}
