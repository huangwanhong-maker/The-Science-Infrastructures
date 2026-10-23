// The lead fields of the forward models, away from the page's main thread.
//
// Message in:  { id, kind: 'eeg' | 'meg', model (EEG: {radius, radii, conductivities, terms}),
//                centre, sensors (mm), orientations (MEG), sources (mm), normals }
// Message out: { id, L (Float64Array sensors x sources) } or { id, error }

import { shellModel, eegLeadField, megLeadField } from '../../../../domains/neuroscience/js/index.js';

self.onmessage = (event) => {
  const { id, kind, model, centre, sensors, orientations, sources, normals } = event.data;
  try {
    const L = kind === 'meg'
      ? megLeadField(centre, sensors, orientations, sources, normals)
      : eegLeadField(shellModel(model), centre, sensors, sources, normals);
    self.postMessage({ id, L }, [L.buffer]);
  } catch (e) {
    self.postMessage({ id, error: e.message });
  }
};
