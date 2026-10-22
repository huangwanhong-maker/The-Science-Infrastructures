// A small globe for choosing where on the Earth to stand: the Earth by day and
// its lights by night, lit by the Sun at the instant shown, in orthographic
// projection. Drag to turn it; click to choose a place.
//
// Orthographic projection centred on (lon0, lat0) (Snyder 1987, eqs. 20-3, 20-4,
// 20-14 to 20-16):
//   x = cos(lat) sin(lon - lon0)
//   y = cos(lat0) sin(lat) - sin(lat0) cos(lat) cos(lon - lon0)
// visible where cos c = sin(lat0) sin(lat) + cos(lat0) cos(lat) cos(lon - lon0) >= 0;
// inversely, with rho = |(x, y)| and c = asin(rho),
//   lat = asin(cos c sin(lat0) + y sin c cos(lat0) / rho)
//   lon = lon0 + atan2(x sin c, rho cos c cos(lat0) - y sin c sin(lat0)).
//
// Snyder, J. P. (1987), Map Projections: A Working Manual, USGS Professional
// Paper 1395, pp. 145-153.

const DEG = Math.PI / 180;
const SOURCE_WIDTH = 1024;

export function orthographic(lon, lat, lon0, lat0) {
  const dl = lon - lon0;
  const cosc = Math.sin(lat0) * Math.sin(lat) + Math.cos(lat0) * Math.cos(lat) * Math.cos(dl);
  return {
    x: Math.cos(lat) * Math.sin(dl),
    y: Math.cos(lat0) * Math.sin(lat) - Math.sin(lat0) * Math.cos(lat) * Math.cos(dl),
    visible: cosc >= 0,
  };
}

export function orthographicInverse(x, y, lon0, lat0) {
  const rho = Math.hypot(x, y);
  if (rho > 1) return null;
  if (rho < 1e-12) return { lon: lon0, lat: lat0 };
  const c = Math.asin(rho), sc = Math.sin(c), cc = Math.cos(c);
  const lat = Math.asin(cc * Math.sin(lat0) + (y * sc * Math.cos(lat0)) / rho);
  const lon = lon0 + Math.atan2(x * sc, rho * cc * Math.cos(lat0) - y * sc * Math.sin(lat0));
  return { lon: ((lon + 3 * Math.PI) % (2 * Math.PI)) - Math.PI, lat };
}

// An equirectangular image as pixels at SOURCE_WIDTH x SOURCE_WIDTH / 2.
function pixelsOf(image) {
  if (!image) return null;
  const canvas = document.createElement('canvas');
  canvas.width = SOURCE_WIDTH;
  canvas.height = SOURCE_WIDTH / 2;
  const g = canvas.getContext('2d', { willReadFrequently: true });
  g.drawImage(image, 0, 0, canvas.width, canvas.height);
  return g.getImageData(0, 0, canvas.width, canvas.height).data;
}

const smooth = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };

export class GlobePicker {
  // onPick(latRad, lonRad) is called when the user clicks a place on the globe.
  constructor(onPick) {
    this.onPick = onPick;
    this.lon0 = 0;
    this.lat0 = 0;
    this.day = null;
    this.night = null;
    this.sun = null;                 // { lon, lat } of the subsolar point, radians
    this.place = null;               // { lon, lat }, radians
    this.canvas = null;
    this.pending = false;
  }

  setImages(dayImage, nightImage) {
    this.day = pixelsOf(dayImage);
    this.night = pixelsOf(nightImage);
    this.request();
  }

  setSun(sun) { this.sun = sun; this.request(); }

  setPlace(place, centre = false) {
    this.place = place;
    if (centre) { this.lon0 = place.lon; this.lat0 = Math.max(-70 * DEG, Math.min(70 * DEG, place.lat)); }
    this.request();
  }

  // Draws into a canvas (the panel is rebuilt now and then, so the canvas may change).
  attach(canvas) {
    if (this.canvas === canvas) return;
    this.canvas = canvas;
    let down = null, moved = false;
    canvas.addEventListener('pointerdown', (e) => {
      down = { x: e.clientX, y: e.clientY, lon0: this.lon0, lat0: this.lat0 };
      moved = false;
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!down) return;
      const R = this.radius() / (window.devicePixelRatio || 1);
      const dx = e.clientX - down.x, dy = e.clientY - down.y;
      if (Math.hypot(dx, dy) > 4) moved = true;
      if (!moved) return;
      this.lon0 = down.lon0 - dx / R;
      this.lat0 = Math.max(-80 * DEG, Math.min(80 * DEG, down.lat0 + dy / R));
      this.request();
    });
    canvas.addEventListener('pointerup', (e) => {
      if (down && !moved) {
        const rect = canvas.getBoundingClientRect();
        const s = canvas.width / rect.width;
        const R = this.radius();
        const x = ((e.clientX - rect.left) * s - canvas.width / 2) / R;
        const y = -((e.clientY - rect.top) * s - canvas.height / 2) / R;
        const p = orthographicInverse(x, y, this.lon0, this.lat0);
        if (p) this.onPick(p.lat, p.lon);
      }
      down = null;
    });
    this.request();
  }

  radius() { return Math.min(this.canvas.width, this.canvas.height) / 2 - 6 * (window.devicePixelRatio || 1); }

  request() {
    if (this.pending || !this.canvas) return;
    this.pending = true;
    requestAnimationFrame(() => { this.pending = false; this.draw(); });
  }

  draw() {
    const canvas = this.canvas;
    if (!canvas || !canvas.isConnected) return;
    const dpr = window.devicePixelRatio || 1;
    const css = canvas.getBoundingClientRect();
    const w = Math.round(css.width * dpr), h = Math.round(css.height * dpr);
    if (!w || !h) return;
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    const g = canvas.getContext('2d');
    const R = this.radius(), cx = w / 2, cy = h / 2;
    const image = g.createImageData(w, h);
    const out = image.data;
    const { lon0, lat0 } = this;
    const sl0 = Math.sin(lat0), cl0 = Math.cos(lat0);
    const sun = this.sun;
    const sx = sun ? Math.cos(sun.lat) * Math.cos(sun.lon) : 0, sy = sun ? Math.cos(sun.lat) * Math.sin(sun.lon) : 0, sz = sun ? Math.sin(sun.lat) : 0;
    const W = SOURCE_WIDTH, H = SOURCE_WIDTH / 2;
    for (let j = 0; j < h; j++) {
      const y = -(j + 0.5 - cy) / R;
      for (let i = 0; i < w; i++) {
        const x = (i + 0.5 - cx) / R;
        const rho2 = x * x + y * y;
        if (rho2 > 1) continue;
        const rho = Math.sqrt(rho2), cc = Math.sqrt(1 - rho2);
        let lat, lon;
        if (rho < 1e-9) { lat = lat0; lon = lon0; } else {
          lat = Math.asin(cc * sl0 + y * cl0);
          lon = lon0 + Math.atan2(x, cc * cl0 - y * sl0);
        }
        const u = ((((lon / (2 * Math.PI) + 0.5) % 1) + 1) % 1) * W;
        const v = Math.min(H - 1, Math.max(0, (0.5 - lat / Math.PI) * H));
        const k = 4 * (Math.floor(v) * W + Math.floor(u));
        let r = 34, gg = 52, b = 92;
        if (this.day) { r = this.day[k]; gg = this.day[k + 1]; b = this.day[k + 2]; }
        // Sunlight: the cosine of the Sun's zenith distance, with a soft twilight.
        let light = 1, lights = 0;
        if (sun) {
          const cosz = Math.cos(lat) * Math.cos(lon) * sx + Math.cos(lat) * Math.sin(lon) * sy + Math.sin(lat) * sz;
          light = smooth(-0.1, 0.12, cosz);
          if (this.night) lights = (1 - smooth(-0.15, 0.02, cosz)) * (this.night[k] * 0.3 + this.night[k + 1] * 0.55 + this.night[k + 2] * 0.15) / 255;
        }
        const limb = 0.72 + 0.28 * cc;
        const dark = 0.16 + 0.84 * light;
        const o = 4 * (j * w + i);
        out[o] = Math.min(255, (r * dark + 36 * (1 - light)) * limb + 255 * lights ** 1.3 * 1.4);
        out[o + 1] = Math.min(255, (gg * dark + 26 * (1 - light)) * limb + 210 * lights ** 1.3 * 1.4);
        out[o + 2] = Math.min(255, (b * dark + 58 * (1 - light)) * limb + 140 * lights ** 1.3 * 1.4);
        // Soft antialiased edge.
        out[o + 3] = Math.round(255 * Math.min(1, (1 - rho) * R * 1.2));
      }
    }
    g.clearRect(0, 0, w, h);
    g.putImageData(image, 0, 0);
    // The graticule every 30 degrees.
    g.save();
    g.lineWidth = 0.6 * dpr;
    g.strokeStyle = 'rgba(246, 236, 244, 0.16)';
    const line = (points) => {
      g.beginPath();
      let pen = false;
      for (const [lo, la] of points) {
        const p = orthographic(lo, la, lon0, lat0);
        if (!p.visible) { pen = false; continue; }
        const px = cx + p.x * R, py = cy - p.y * R;
        if (pen) g.lineTo(px, py); else g.moveTo(px, py);
        pen = true;
      }
      g.stroke();
    };
    for (let lo = -180; lo < 180; lo += 30) line(Array.from({ length: 91 }, (_, k) => [lo * DEG, (-90 + 2 * k) * DEG]));
    for (let la = -60; la <= 60; la += 30) line(Array.from({ length: 181 }, (_, k) => [(-180 + 2 * k) * DEG, la * DEG]));
    // The rim, in the atlas's orchid.
    g.beginPath();
    g.arc(cx, cy, R, 0, 2 * Math.PI);
    g.lineWidth = 1.2 * dpr;
    g.strokeStyle = 'rgba(202, 141, 190, 0.7)';
    g.stroke();
    // The point under the Sun.
    if (sun) {
      const p = orthographic(sun.lon, sun.lat, lon0, lat0);
      if (p.visible) {
        g.beginPath();
        g.arc(cx + p.x * R, cy - p.y * R, 4 * dpr, 0, 2 * Math.PI);
        g.fillStyle = '#f6c47e';
        g.fill();
      }
    }
    // The place chosen.
    if (this.place) {
      const p = orthographic(this.place.lon, this.place.lat, lon0, lat0);
      if (p.visible) {
        const px = cx + p.x * R, py = cy - p.y * R;
        g.beginPath();
        g.arc(px, py, 7 * dpr, 0, 2 * Math.PI);
        g.lineWidth = 2 * dpr;
        g.strokeStyle = '#CA8DBE';
        g.stroke();
        g.beginPath();
        g.arc(px, py, 2.4 * dpr, 0, 2 * Math.PI);
        g.fillStyle = '#fdf4fa';
        g.fill();
      }
    }
    g.restore();
  }
}
