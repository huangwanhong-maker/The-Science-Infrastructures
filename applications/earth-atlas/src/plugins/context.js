// What a plug-in sees of the atlas while it runs: a narrow view, so that
// scripts depend on little, and so that the same view can be made from files
// when a plug-in is checked without a page.
export class Context {
  constructor(library, { time = Date.now(), pinned = null, progress = () => {} } = {}) {
    this.library = library;
    this.time = time;          // the instant shown, ms UTC
    this.pinned = pinned;      // the chosen place, [lon, lat], or null
    this.progress = progress;  // (fraction, text) -> the panel's progress bar
  }

  // The metadata of every field product in the library.
  fields() { return [...this.library.fieldMeta.values()]; }

  // A field product by name.
  field(name) { return this.library.field(name); }

  // The geography product (Natural Earth): coasts, countries, places.
  geography() { return this.library.geography(); }

  // The library's index entries of a kind: 'fields', 'events', 'observations', 'tracks', 'regions', ...
  list(kind) { return this.library.list(kind); }

  // An event set by name (earthquakes, mechanisms, disasters): { fields, events: [[row]], window, source }.
  events(name) { return this.library.eventSet(name); }

  // An observation set by name (stations, buoys, sensors; data/README.md §4): { variables, fields, rows, window, source }.
  observations(name) { return this.library.observationSet(name); }
}
