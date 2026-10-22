// A domain's function called through the service, from any JavaScript (a page
// elsewhere, or Node): the Sun's apparent place, the subsolar point now, and the
// list of what can be called.
//
//   node api/examples/js/compute.mjs

import { Client } from '../../js/science.mjs';

const api = new Client();
const jd = Date.now() / 86400000 + 2440587.5;
console.log('astronomy.sunPosition', await api.compute('astronomy', 'sunPosition', jd));
console.log('astronomy.subsolarPoint', await api.compute('astronomy', 'subsolarPoint', Date.now()));
for (const unit of await api.units()) {
  console.log(`  ${unit.unit.padEnd(14)} ${unit.kind.padEnd(10)} ${unit.functions.length} functions: ${unit.functions.slice(0, 4).map((f) => f.name).join(', ')}, ...`);
}
