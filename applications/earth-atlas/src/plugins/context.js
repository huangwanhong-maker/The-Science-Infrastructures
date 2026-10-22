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
}
