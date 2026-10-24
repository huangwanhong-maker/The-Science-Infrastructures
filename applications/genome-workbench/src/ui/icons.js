// Line icons, drawn for the workbench on a 24-unit grid with one stroke (styles.css, svg.icon): the
// views, the groups of tools, and the actions. icon(name, size) returns an SVG element.

const SVG = 'http://www.w3.org/2000/svg';

const PATHS = {
  // views
  genome: '<circle cx="12" cy="12" r="7.5"/><path d="M12 4.5a7.5 7.5 0 0 1 7.3 5.8"/><path d="M18.9 15a7.5 7.5 0 0 1-9.4 4" stroke-width="2.6"/><circle cx="12" cy="12" r="3"/>',
  align: '<path d="M4 7h4l2 3h4l2-3h4"/><path d="M4 17h5l1.5-2.5h3L15 17h5"/><path d="M8 7v10M16 7v10" stroke-dasharray="1.6 2"/>',
  index: '<path d="M5 4.5h10l4 4V19.5H5z"/><path d="M15 4.5v4h4"/><path d="M8.5 12h7M8.5 15h4"/><circle cx="16.5" cy="16.5" r="2.2"/><path d="m18.1 18.1 1.8 1.8"/>',
  reads: '<path d="M3.5 7h7M13 7h7.5M5.5 11h9M16.5 11h4M3.5 15h5M10.5 15h8M6.5 19h11"/>',
  assembly: '<circle cx="5" cy="12" r="2"/><circle cx="12" cy="6" r="2"/><circle cx="12" cy="18" r="2"/><circle cx="19" cy="12" r="2"/><path d="m6.6 10.8 3.8-3.4M13.6 7.2l3.8 3.6M6.6 13.2l3.8 3.4M13.6 16.8l3.8-3.6"/>',
  expression: '<path d="M4 19.5V9M9 19.5V5M14 19.5v-8M19 19.5V7"/><path d="M3 19.5h18"/>',
  trees: '<path d="M4 5h4v14H4M8 12h5M13 8.5h3v7h-3M16 8.5h4M16 15.5h4M8 5h12M8 19h12"/>',
  populations: '<circle cx="8" cy="9" r="2.6"/><circle cx="16" cy="9" r="2.6"/><path d="M3.5 18.5c.7-3 2.4-4.5 4.5-4.5s3.8 1.5 4.5 4.5M11.5 18.5c.7-3 2.4-4.5 4.5-4.5s3.8 1.5 4.5 4.5"/>',
  run: '<rect x="6" y="6" width="12" height="12" rx="2"/><rect x="9.5" y="9.5" width="5" height="5" rx="1"/><path d="M9 3v3M15 3v3M9 18v3M15 18v3M3 9h3M3 15h3M18 9h3M18 15h3"/>',
  toolbox: '<rect x="3.5" y="8" width="17" height="11" rx="2"/><path d="M9 8V6.2a1.7 1.7 0 0 1 1.7-1.7h2.6A1.7 1.7 0 0 1 15 6.2V8"/><path d="M3.5 12.5h17M10.5 11.2v2.6M13.5 11.2v2.6"/>',
  data: '<ellipse cx="12" cy="6" rx="7.5" ry="2.8"/><path d="M4.5 6v12c0 1.5 3.4 2.8 7.5 2.8s7.5-1.3 7.5-2.8V6"/><path d="M4.5 12c0 1.5 3.4 2.8 7.5 2.8s7.5-1.3 7.5-2.8"/>',
  // tool groups
  sequence: '<path d="M7 3c0 4.5 10 4.5 10 9S7 16.5 7 21M17 3c0 4.5-10 4.5-10 9s10 4.5 10 9"/><path d="M8.5 6h7M9.5 9.5h5M9.5 14.5h5M8.5 18h7"/>',
  annotation: '<path d="M4 15h16"/><rect x="5" y="11.5" width="5" height="3.5" rx="1"/><rect x="12" y="15" width="6" height="3.5" rx="1"/><path d="M7.5 11.5V7l2 1.5L11 6"/>',
  sketch: '<circle cx="9" cy="12" r="5.5"/><circle cx="15" cy="12" r="5.5"/>',
  statistics: '<path d="M3.5 19.5c2.5 0 3.3-13 8.5-13s6 13 8.5 13"/><path d="M12 6.5v13" stroke-dasharray="1.6 2"/>',
  variant: '<path d="M4 8h16M4 16h16"/><path d="M11 6.5 13 9.5M11 9.5l2-3" /><circle cx="12" cy="16" r="1.6"/>',
  file: '<path d="M6 3.5h8l4 4v13H6z"/><path d="M14 3.5v4h4"/>',
  upload: '<path d="M12 15.5V4.5M7.5 9 12 4.5 16.5 9"/><path d="M4.5 15.5v3a1.5 1.5 0 0 0 1.5 1.5h12a1.5 1.5 0 0 0 1.5-1.5v-3"/>',
  download: '<path d="M12 4.5v11M7.5 11 12 15.5 16.5 11"/><path d="M4.5 15.5v3a1.5 1.5 0 0 0 1.5 1.5h12a1.5 1.5 0 0 0 1.5-1.5v-3"/>',
  // actions
  search: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 5 5"/>',
  command: '<path d="M9 6.5A2.5 2.5 0 1 0 6.5 9H9V6.5zM15 6.5A2.5 2.5 0 1 1 17.5 9H15V6.5zM9 17.5A2.5 2.5 0 1 1 6.5 15H9v2.5zM15 17.5a2.5 2.5 0 1 0 2.5-2.5H15v2.5zM9 9h6v6H9z"/>',
  play: '<path d="M7.5 5.5v13l11-6.5z"/>',
  copy: '<rect x="8.5" y="8.5" width="11" height="11" rx="2"/><path d="M15.5 8.5V6a1.5 1.5 0 0 0-1.5-1.5H6A1.5 1.5 0 0 0 4.5 6v8A1.5 1.5 0 0 0 6 15.5h2.5"/>',
  info: '<circle cx="12" cy="12" r="8.5"/><path d="M12 11v5.5M12 7.6v.4"/>',
  warn: '<path d="M12 4 21 19.5H3z"/><path d="M12 10v4.5M12 17v.3"/>',
  check: '<path d="m5 12.5 4.5 4.5L19 7.5"/>',
  x: '<path d="M6 6l12 12M18 6 6 18"/>',
  plus: '<path d="M12 5v14M5 12h14"/>',
  trash: '<path d="M4.5 7h15M9.5 7V4.8h5V7M6.5 7l1 12.5h9l1-12.5"/>',
  refresh: '<path d="M19 12a7 7 0 1 1-2.1-5"/><path d="M19 4.5V9h-4.5"/>',
  sidebar: '<rect x="3.5" y="4.5" width="17" height="15" rx="2"/><path d="M9.5 4.5v15"/>',
  sun: '<circle cx="12" cy="12" r="4"/><path d="M12 2.5v2M12 19.5v2M4.6 4.6l1.4 1.4M18 18l1.4 1.4M2.5 12h2M19.5 12h2M4.6 19.4 6 18M18 6l1.4-1.4"/>',
  moon: '<path d="M19.5 14.5A8 8 0 0 1 9.5 4.5a8 8 0 1 0 10 10z"/>',
  chevron: '<path d="m9 6 6 6-6 6"/>',
  left: '<path d="m15 6-6 6 6 6"/>',
  up: '<path d="m6 15 6-6 6 6"/>',
  'zoom-in': '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 5 5M8 10.5h5M10.5 8v5"/>',
  'zoom-out': '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 5 5M8 10.5h5"/>',
  fit: '<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/><path d="M8 12h8" stroke-dasharray="1.6 2"/>',
  target: '<circle cx="12" cy="12" r="6.5"/><circle cx="12" cy="12" r="1.8"/><path d="M12 2.5v3.5M12 18v3.5M2.5 12H6M18 12h3.5"/>',
  down: '<path d="m6 9 6 6 6-6"/>',
  external: '<path d="M14 4.5h5.5V10M19.5 4.5 11 13"/><path d="M18 14v4.5a1.5 1.5 0 0 1-1.5 1.5h-11A1.5 1.5 0 0 1 4 18.5v-11A1.5 1.5 0 0 1 5.5 6H10"/>',
  activity: '<path d="M3 12h4l2.5-6 5 12 2.5-6h4"/>',
  history: '<path d="M4.5 12a7.5 7.5 0 1 0 2.2-5.3"/><path d="M4.5 4.5V9H9"/><path d="M12 8v4.5l3 2"/>',
  spark: '<path d="M12 3.5 13.8 10.2 20.5 12 13.8 13.8 12 20.5 10.2 13.8 3.5 12 10.2 10.2z"/>',
  layers: '<path d="M12 4 20.5 8.5 12 13 3.5 8.5z"/><path d="m3.5 12.5 8.5 4.5 8.5-4.5M3.5 16.5 12 21l8.5-4.5"/>',
  gpu: '<rect x="3.5" y="7" width="17" height="10" rx="1.5"/><circle cx="9" cy="12" r="2.5"/><path d="M14 10.5h4M14 13.5h4M6 17v2.5M10 17v2.5M14 17v2.5"/>',
  keyboard: '<rect x="3" y="6.5" width="18" height="11" rx="2"/><path d="M6.5 10h1M10 10h1M13.5 10h1M17 10h.5M7 14h10"/>',
};

export function icon(name, size = '') {
  const svg = document.createElementNS(SVG, 'svg');
  svg.setAttribute('viewBox', '0 0 24 24');
  svg.setAttribute('class', `icon${size ? ` ${size}` : ''}`);
  svg.setAttribute('aria-hidden', 'true');
  svg.innerHTML = PATHS[name] || PATHS.spark;
  return svg;
}

export const ICONS = Object.keys(PATHS);
