// The hazards named for a place: ThinkHazard!'s levels for an administrative unit, by a name anywhere on Earth. It classifies a unit, not a
// street: a level says how likely damaging events are somewhere in the unit within a decade or so.

const LEVEL = { HIG: ['high', 3], MED: ['medium', 2], LOW: ['low', 1], VLO: ['very low', 0] };

export async function render(root, ctx) {
  const { area, dom: { h } } = ctx;
  const input = h('input', { type: 'search', class: 'wide', value: area.short || area.name || '', 'aria-label': 'an administrative unit', placeholder: 'A district, a province, a country' });
  const list = h('div', { class: 'choices' });
  const out = h('div', {});
  root.append(h('h1', {}, 'The hazards named for this place'),
    h('p', { class: 'muted' }, 'Levels for an administrative unit, from ThinkHazard! (GFDRR): how likely it is that a damaging event of each kind strikes somewhere in the unit. A first look, to be followed by a study of the place.'),
    h('form', { class: 'search-row', onsubmit: (e) => { e.preventDefault(); search(); } }, input, h('button', { class: 'primary' }, 'Find')), list, out);
  async function search() {
    const q = input.value.trim();
    if (!q) return;
    list.replaceChildren(h('p', { class: 'busy' }, 'asking…'));
    const r = await fetch(`/world/thinkhazard?name=${encodeURIComponent(q)}`);
    const a = await r.json();
    if (!r.ok) { list.replaceChildren(h('p', { class: 'error-box' }, a.message)); return; }
    if (!a.units?.length) { list.replaceChildren(h('p', { class: 'muted' }, 'No unit of that name: its units are districts, provinces and countries, and a city often lies within one of another name.')); return; }
    list.replaceChildren(...a.units.slice(0, 12).map((u) => h('button', { class: 'choice', onclick: () => report(u) }, h('b', {}, u.admin2 || u.admin1 || u.admin0), h('span', {}, [u.admin1, u.admin0].filter(Boolean).join(', ')))));
    if (a.units.length === 1) report(a.units[0]);
  }
  async function report(u) {
    out.replaceChildren(h('p', { class: 'busy' }, 'asking…'));
    const code = u.code ?? u.admin2_code ?? u.admin1_code ?? u.admin0_code;
    const r = await fetch(`/world/thinkhazard?code=${encodeURIComponent(code)}`);
    const a = await r.json();
    if (!r.ok) { out.replaceChildren(h('p', { class: 'error-box' }, a.message)); return; }
    const rows = [...a.report].map((x) => ({ name: x.hazardtype?.hazardtype, level: x.hazardlevel?.mnemonic })).filter((x) => x.name).sort((p, q) => (LEVEL[q.level]?.[1] ?? -1) - (LEVEL[p.level]?.[1] ?? -1));
    out.replaceChildren(h('h2', {}, [u.admin2, u.admin1, u.admin0].filter(Boolean).join(', ')),
      h('div', { class: 'hazards' }, rows.map((x) => h('div', { class: `hazard level-${LEVEL[x.level]?.[1] ?? 'none'}` }, h('span', {}, x.name), h('i', {}, h('b', {})), h('em', {}, LEVEL[x.level]?.[0] || 'no data')))),
      h('p', { class: 'note' }, 'ThinkHazard!, GFDRR. Its levels come from global datasets at coarse scales, and are not a substitute for a study of the place.'));
  }
  search();
}
