// The workbench's clock: a civil date and hour in the shown city's time zone, the UTC instant they are (the
// zone's offset from the browser's time-zone database, so daylight saving is the zone's own), the hour of
// the civil week (Monday 00:00 = 0, as the layers of hours read it), and the Sun's altitude and azimuth there
// (the astronomy domain's subsolar point). Layers listen to it.

import { subsolarPoint, sunAltAz } from '/domains/astronomy/js/index.js';

// The zone's offset from UTC, hours, at an instant (ms).
export function offsetHours(ms, timezone) {
  const parts = Object.fromEntries(new Intl.DateTimeFormat('en-GB', { timeZone: timezone, hourCycle: 'h23', year: 'numeric', month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit', second: '2-digit' })
    .formatToParts(new Date(ms)).map((p) => [p.type, p.value]));
  const asUtc = Date.UTC(+parts.year, +parts.month - 1, +parts.day, +parts.hour, +parts.minute, +parts.second);
  return Math.round((asUtc - Math.floor(ms / 1000) * 1000) / 60000) / 60;
}

// The UTC instant of a civil date (YYYY-MM-DD) and hour (fractional) in a zone.
export function instantOf(date, hour, timezone) {
  const [y, m, d] = date.split('-').map(Number);
  const naive = Date.UTC(y, m - 1, d) + hour * 3600000;
  let t = naive - offsetHours(naive, timezone) * 3600000;
  t = naive - offsetHours(t, timezone) * 3600000;
  return t;
}

export const DAYS = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'];

export function createClock(area, { date, hour = 14 } = {}) {
  const today = new Date();
  const state = {
    area,
    date: date || today.toISOString().slice(0, 10),
    hour,
    listeners: [],
    get weekday() { const [y, m, d] = this.date.split('-').map(Number); return (new Date(Date.UTC(y, m - 1, d)).getUTCDay() + 6) % 7; },
    get weekHour() { return 24 * this.weekday + this.hour; },
    get instant() { return instantOf(this.date, this.hour, this.area.timezone); },
    get offset() { return offsetHours(this.instant, this.area.timezone); },
    get month() { return +this.date.slice(5, 7); },
    get sun() {
      const a = sunAltAz(this.area.lon, this.area.lat, subsolarPoint(this.instant));
      return { altitude: a.altitude, azimuth: a.azimuth };
    },
    label() {
      const h = Math.floor(this.hour), mi = Math.round((this.hour - h) * 60);
      const off = this.offset;
      return `${DAYS[this.weekday].slice(0, 3)} ${this.date} ${String(h).padStart(2, '0')}:${String(mi).padStart(2, '0')} (UTC${off >= 0 ? '+' : '−'}${Math.abs(off)})`;
    },
    set(changes) {
      Object.assign(this, changes);
      for (const f of this.listeners) f(this);
    },
    on(f) { this.listeners.push(f); },
  };
  return state;
}
