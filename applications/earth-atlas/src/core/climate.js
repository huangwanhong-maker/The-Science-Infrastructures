// The climate of a place from its normals: twelve monthly mean temperatures
// T (degrees C) and precipitation totals P (mm), January first.
//
// Provisional stand-in for core/domains/climate (zones; laboratory ADR 0017).
//
// The Köppen–Geiger classification follows Beck et al. (2018), Table 2, who
// follow Peel et al. (2007), which differs from Köppen (1936) in three ways:
// temperate (C) and cold (D) are divided at 0 degrees C rather than -3; the
// desert and steppe thresholds depend on whether 70 % of the precipitation falls
// in summer or in winter; and where a C or D climate meets both the dry-summer
// (s) and the dry-winter (w) criteria, it is s if more precipitation falls in
// winter, w otherwise. The arid group B takes precedence over all others.
//
// Summer is the warmer of the two half-years April-September and October-March,
// winter the other. The classification is thus defined by temperature, not by
// hemisphere or calendar: shifting a climate by six months leaves its class
// unchanged whenever its two half-years differ in temperature (the classes are
// invariant under the action of Z/2 on the months by a half-year shift), as
// tests/run-tests.mjs checks.
//
// Beck, H. E., et al. (2018), Present and future Köppen-Geiger climate
//   classification maps at 1-km resolution, Scientific Data 5, 180214.
// Beck, H. E., et al. (2023), High-resolution (1 km) Köppen-Geiger maps for
//   1901-2099 based on constrained CMIP6 projections, Scientific Data 10, 724.
// Peel, M. C., Finlayson, B. L. and McMahon, T. A. (2007), Updated world map of
//   the Köppen-Geiger climate classification, HESS 11, 1633-1644.

// Codes in the order of Beck et al.'s legend: the value of a class is its index + 1.
export const KOPPEN_CODES = ['Af', 'Am', 'Aw', 'BWh', 'BWk', 'BSh', 'BSk', 'Csa', 'Csb', 'Csc', 'Cwa', 'Cwb', 'Cwc',
  'Cfa', 'Cfb', 'Cfc', 'Dsa', 'Dsb', 'Dsc', 'Dsd', 'Dwa', 'Dwb', 'Dwc', 'Dwd', 'Dfa', 'Dfb', 'Dfc', 'Dfd', 'ET', 'EF'];
export const koppenValue = (code) => KOPPEN_CODES.indexOf(code) + 1;

// April to September, as indices from January.
const AMJJAS = [3, 4, 5, 6, 7, 8];
const ONDJFM = [9, 10, 11, 0, 1, 2];
const sum = (a, idx) => idx.reduce((s, i) => s + a[i], 0);

// The quantities the classification uses (Beck et al. 2018, Table 2, note a).
export function koppenQuantities(T, P) {
  const tA = sum(T, AMJJAS), tO = sum(T, ONDJFM);
  // Summer is the warmer half-year; an exact tie (no seasons) counts April-September as summer.
  const [summer, winter] = tA >= tO ? [AMJJAS, ONDJFM] : [ONDJFM, AMJJAS];
  const MAT = (tA + tO) / 12;
  const MAP = sum(P, AMJJAS) + sum(P, ONDJFM);
  const Psummer = sum(P, summer), Pwinter = sum(P, winter);
  const pick = (idx, f) => f(...idx.map((i) => P[i]));
  let Pthreshold = 2 * MAT + 14;
  if (Pwinter > 0.7 * MAP) Pthreshold = 2 * MAT;
  else if (Psummer > 0.7 * MAP) Pthreshold = 2 * MAT + 28;
  return {
    MAT, MAP, Pthreshold, Psummer, Pwinter,
    Thot: Math.max(...T), Tcold: Math.min(...T),
    Tmon10: T.filter((t) => t > 10).length,
    Pdry: Math.min(...P),
    Psdry: pick(summer, Math.min), Pswet: pick(summer, Math.max),
    Pwdry: pick(winter, Math.min), Pwwet: pick(winter, Math.max),
    summer: summer === AMJJAS ? 'April–September' : 'October–March',
  };
}

// The Köppen–Geiger class of a climate, as its code ('Cfb'), or null where a value is missing.
export function koppen(T, P) {
  for (let i = 0; i < 12; i++) if (!(Number.isFinite(T[i]) && Number.isFinite(P[i]))) return null;
  const q = koppenQuantities(T, P);
  if (q.MAP < 10 * q.Pthreshold) {
    return `B${q.MAP < 5 * q.Pthreshold ? 'W' : 'S'}${q.MAT >= 18 ? 'h' : 'k'}`;
  }
  if (q.Tcold >= 18) {
    if (q.Pdry >= 60) return 'Af';
    return q.Pdry >= 100 - q.MAP / 25 ? 'Am' : 'Aw';
  }
  if (q.Thot <= 10) return q.Thot > 0 ? 'ET' : 'EF';
  const group = q.Tcold > 0 ? 'C' : 'D';
  const s = q.Psdry < 40 && q.Psdry < q.Pwwet / 3;
  const w = q.Pwdry < q.Pswet / 10;
  const second = s && w ? (q.Pwinter > q.Psummer ? 's' : 'w') : s ? 's' : w ? 'w' : 'f';
  let third;
  if (q.Thot >= 22) third = 'a';
  else if (q.Tmon10 >= 4) third = 'b';
  else if (group === 'D' && q.Tcold < -38) third = 'd';
  else third = 'c';
  return `${group}${second}${third}`;
}

// The mean temperature of a month from the means of its daily maxima and minima.
export const meanTemperature = (tmax, tmin) => tmax.map((x, i) => (x + tmin[i]) / 2);

// Dry months in the sense of Gaussen: P <= 2 T, the rule drawn by the
// climate diagrams of Walter and Lieth (1960), where precipitation is plotted
// at 20 mm to the 10 degrees of temperature.
export const dryMonths = (T, P) => T.map((t, i) => P[i] <= 2 * t);

// Classes of a whole grid of normals: Uint8Array of class values (0 where missing).
// tmax, tmin, ppt: arrays of 12 Float32Arrays (one per month) on the same grid.
export function koppenGrid(tmax, tmin, ppt) {
  const n = tmax[0].length;
  const out = new Uint8Array(n);
  const T = new Array(12), P = new Array(12);
  for (let k = 0; k < n; k++) {
    let ok = true;
    for (let m = 0; m < 12 && ok; m++) {
      T[m] = (tmax[m][k] + tmin[m][k]) / 2;
      P[m] = ppt[m][k];
      ok = T[m] === T[m] && P[m] === P[m];
    }
    if (ok) out[k] = koppenValue(koppen(T, P));
  }
  return out;
}
