# 0068. The Economy Atlas as a Scholarly Instrument: Porcelain and Gold, and Figures for Publication

- **Status:** Accepted, provisionally (the colours, and the proposal of §3, are the user's to confirm); amended by [0069](0069-names-of-places-and-a-2d-gis.md), on which its map now stands
- **Date:** 2026-11-04
- **Amends:** [ADR 0046](0046-the-economy-atlas-data-and-plugins.md), [ADR 0056](0056-nomos-and-a-toolbox-of-quantitative-economics.md)
- **Decided by:** the user handed in a brief for the interface's redesign ("a sophisticated scientific research instrument with the visual refinement of an exceptional contemporary cartographic publication") and asked that only its really good points be taken: "Not need fully trust the design." Claude chose which, adapted some, left others, and wrote this record.

## 1. Context

After ADR 0056 the Economy Atlas computed much and looked unconsidered: gold washed over its surfaces, its controls were pills, pale italics were its default type, and nothing on its map said what projection it was or where its numbers came from. The brief named these faults well. It also asked for things that would be speculative now (projections switched, spatial queries, PDF) and for a second typeface, which would change a convention of CLAUDE.md (§7: every interface is set in Source Serif 4).

## 2. Decision

### 2.1 Three Kinds of Colour, Kept Apart

- **The interface**: porcelain surfaces (`#FAFAF7`, `#F3F3EF`), graphite ink (`#2B2F33`, secondary `#555D63`, metadata `#6B7378`, each at least 4.5:1), hairline borders, and a champagne gold (`#B08A3E`) only for what is chosen, in focus or emphasised. As text, gold is a darker `#7D5F1E` (5.6:1); the brief's `#B99149` is 2.9:1, too faint for text. Rose (`#E0607E`), the user's minor colour, stays: in the data, and in the interface for errors.
- **The base of the map**: a white sea, land without data a neutral grey, the coasts in the canvas's colour, the graticule nearly gone; quieter than anything drawn on them.
- **The data** (`palette.js`) read only the data's tokens: the validated colours of ADR 0046 and 0056, and the ramps' own gold, separate from the interface's.

Night is graphite (`#15181B`), with the same separation.

### 2.2 Type and Geometry

One family, Source Serif 4, by its optical sizes; lining tabular figures wherever numbers are read; section heads in small capitals on a hairline; italics only for citations, captions and notes about method. Corners of 3 px, controls 28 px high, hairlines rather than shadows; the view chosen underlined in gold; the zoom a compact column; no gradients, no blur.

### 2.3 The Map as the Subject

- The dock and the inspector **fold away** (kept in the browser), and the map takes the room.
- The map's **instruments** name the projection (Equal Earth), the longitude and latitude under the pointer, by the geography domain's inverse projection, and the zoom.
- Each **legend** ends with its source, year, unit and projection.
- Country names are set in small capitals; a chosen place is outlined in graphite, not gold.
- **The address keeps the view**: dataset, year, colouring, arcs and their number, the layers' switches and the place chosen, so that any view, and any figure, can be linked and drawn again.

### 2.4 Figures for Publication

The map's numbers were already saved as CSV. Its drawing is now saved as a figure composed for publication: a title and a line saying what is coloured and drawn; the map as shown, its arcs solid and headed rather than moving; the legend's ramp, ends and keys; and a line of source, projection, date and the link that draws it again. It is SVG with Source Serif 4 written into it, or PNG at twice its 1800 pixels. The figure is drawn by the map's own renderer, cloned with its styles, so the analytical view and the figure cannot disagree.

## 3. Alternatives Considered

- **A second typeface for controls and figures** (the brief's suggestion). Not taken: it would change a convention of CLAUDE.md, which needs the user's word. *Proposed*: Source Sans 3, Source Serif's companion, for controls, tables and the map's instruments, as Cosmopoiesis took IBM Plex (its decision 0006), if the user wishes it.
- **The brief's tokens as given.** Adapted: its gold for text fails contrast, and its cool grey surface (`#EFF0ED`) was warmed a little toward the palette the user chose.
- **Projections and coordinate systems switched, spatial queries, PDF export, saved configurations.** Not built: the address is the saved configuration; the readout and the legend's projection line leave room for the rest when the work needs them.

## 4. Consequences

- The Economy Atlas's visual style in CLAUDE.md §7 becomes porcelain and graphite with champagne gold, rose its minor colour.
- The gallery's thirteen pictures of the Economy Atlas were taken again in it.
- The application's checks follow the ramps' new neutral (`#F1F0EA`).
