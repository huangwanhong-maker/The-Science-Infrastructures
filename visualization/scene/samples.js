// A field read for drawing: the answers of an engine's field operations decoded into typed arrays (the
// protocol of Cosmopoiesis's fields.hpp, which any engine serving physis's fields may speak).
//
// A sampling's answer carries, for its `size` points, `components` numbers each (1 for a scalar or a
// reduction, 3 for a vector, interleaved), NaN where there is no value; each point's quality (physis's
// order: valid, outside, missing, invalid) and, by index, its reason (when refused) or note (when valid
// but qualified), with the counts of each. Encoded "f32", the numbers are base64 of little-endian
// float32, the qualities of bytes and the indices of uint16; encoded "json", they are arrays, null where
// there is no value. A point without a value is never drawn as a zero (formalism §3.2).

export const QUALITIES = ['valid', 'outside', 'missing', 'invalid'];

// Base64 (RFC 4648 §4) to bytes, in the browser and in Node alike.
export function bytesOf(text) {
  const s = atob(text), out = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}
// Typed arrays of base64 (little-endian, as every browser's and Node's typed arrays are here).
export const float32Of = (text) => new Float32Array(bytesOf(text).buffer);
export const uint16Of = (text) => new Uint16Array(bytesOf(text).buffer);

// The samples of a field.sample answer.
export function samplesOf(answer) {
  if (answer.error) throw new Error(answer.error);
  const size = answer.size, components = answer.components;
  const values = answer.values_f32 !== undefined ? float32Of(answer.values_f32) : Float32Array.from(answer.values, (v) => (v === null ? NaN : v));
  const quality = answer.quality_u8 !== undefined ? bytesOf(answer.quality_u8) : Uint8Array.from(answer.quality);
  const why = answer.why_u16 !== undefined ? uint16Of(answer.why_u16) : Uint16Array.from(answer.why);
  if (values.length !== size * components || quality.length !== size || why.length !== size)
    throw new Error(`field.sample: ${values.length} numbers, ${quality.length} qualities and ${why.length} indices for ${size} points of ${components}`);
  return {
    field: answer.field, sampling: answer.sampling, time: answer.time, reduction: answer.reduction,
    size, components, values, quality, why,
    reasons: answer.reasons || [], notes: answer.notes || [], counts: answer.counts, range: answer.range || null,
  };
}

export const isValid = (s, i) => s.quality[i] === 0;
// Why a point has no value, or null; the note on a valid point, or null.
export const reasonOf = (s, i) => (s.quality[i] !== 0 && s.why[i] ? s.reasons[s.why[i] - 1].text : null);
export const noteOf = (s, i) => (s.quality[i] === 0 && s.why[i] ? s.notes[s.why[i] - 1].text : null);

// The scalar drawn at each point: the value of a scalar, the length of a vector; NaN where none.
export function scalarOf(s) {
  if (s.components === 1) return s.values;
  const out = new Float32Array(s.size);
  for (let i = 0; i < s.size; i++) {
    const x = s.values[3 * i], y = s.values[3 * i + 1], z = s.values[3 * i + 2];
    out[i] = Math.hypot(x, y, z);
  }
  return out;
}

// The curves of a field.curves answer: each one's points (x, y, z interleaved, on the field's frame's
// axes), |v| along it, its seed, and how each end ended.
export function curvesOf(answer) {
  if (answer.error) throw new Error(answer.error);
  return answer.curves.map((c) => {
    if (!c.seeded) return { seeded: false, quality: c.quality, reason: c.reason, points: new Float32Array(0), magnitude: new Float32Array(0) };
    const points = c.points_f32 !== undefined ? float32Of(c.points_f32) : Float32Array.from(c.points.flat());
    const magnitude = c.magnitude_f32 !== undefined ? float32Of(c.magnitude_f32) : Float32Array.from(c.magnitude, (v) => (v === null ? NaN : v));
    return { seeded: true, points, magnitude, seedIndex: c.seed_index, length_m: c.length_m, backward: c.backward, forward: c.forward };
  });
}
