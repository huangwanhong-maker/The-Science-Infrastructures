# 0017. Climate Classes, Climate Normals, and Two New Ways of Acquiring Data

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review.
- **Date:** 2026-10-22
- **Decided by:** the user asked for more complex visualization of geography, and for creativity. Claude chose the design below.

## 1. Context

The Earth Atlas showed weather: forecasts at an instant. Geography also asks what a place is like over the years, which is its climate. The classic answers are two. One is the Köppen–Geiger classification, the map of climates every atlas carries. The other is the climate diagram of Walter and Lieth, the normal year of one place in one figure.

Both need data the library did not have:
- maps of climate classes, which are categories rather than quantities;
- monthly normals, which describe a typical year rather than a sequence of instants.

Both sources are served in ways the library could not yet use:
- Beck et al.'s maps are members of a 131 MB zip archive.
- TerraClimate's normals sit on a THREDDS server whose NetCDF Subset Service answers with empty files; its OPeNDAP service works.

## 2. Decision

1. **Members of remote archives.** A file of a `files` dataset may name a `member` of a zip archive. The fetcher reads the archive's directory and that member by HTTP range requests (`data/acquisition/remotezip.py`). The twelve maps and their legend take about 2.2 MB of transfer instead of 131 MB. Checksums are recorded per member.
2. **A new acquisition kind, `opendap`.**
   - It reads strided subsets of gridded datasets from DAP2 servers (`data/acquisition/opendap.py`), using their attributes (DAS), structure (DDS) and binary values.
   - It stores each subset as a netCDF classic file, through a new writer in `data/formats/netcdf3.py`, so the product builders read it like any other.
   - Registered: `terraclimate-normals`, the 1991–2020 normals of daily maximum and minimum temperature and of precipitation (CC0), every 1/3°.
3. **Climatologies in field products.**
   - A manifest may declare a `climatology` period. Its times then become the middles of the months of the period's middle year (CF Conventions, §7.4), and the product records the period.
   - An application shows a climatology by the month of the time shown, in any year. Climatologies never extend the timeline.
4. **A new product kind, `classes`.**
   - `data/products/class_raster.py` turns GeoTIFF maps of categories into one byte per cell. It writes the legend (code, name, the authors' colour) and the share of land in each class, weighted by area.
   - The GeoTIFF georeferencing is read and checked, not assumed.
   - Registered: `koppen-geiger`, Beck et al. (2023), 0.1°:
     - 1901–1930, 1931–1960, 1961–1990 and 1991–2020;
     - 2041–2070 and 2071–2099, each under SSP1-2.6, SSP2-4.5, SSP3-7.0 and SSP5-8.5.
5. **Classes in the Earth Atlas.**
   - Class maps are a kind of surface, chosen by period, scenario and source.
   - Classes are never interpolated: the shader reads the cell that contains each pixel, from an integer texture.
   - Any two maps can be compared. Cells whose class changed keep their colour, and unchanged cells fade.
   - The conventional colours (Peel et al. 2007) keep their hues but are softened in OKLab, to sit with the atlas's palette.
6. **Our own classification.**
   - The Köppen–Geiger rules are implemented in `src/core/climate.js`, a stand-in for `core/domains/climate`.
   - The atlas classifies the TerraClimate normals itself, as a computed class map, and offers it beside the published one.
   - The rules are tested on constructed climates and on their symmetry: invariance under (S₆ × S₆) ≀ ℤ/2.
   - The result is verified against Beck et al.'s map: 89.1 % of land cells agree, 95.9 % by main group.
   - The theory is in [theory/climate/koppen-geiger-classification.md](../../theory/climate/koppen-geiger-classification.md).
7. **The climate diagram.** The place card draws the normal year of any place on land after Walter and Lieth (1960), at a fixed scale, so that diagrams compare.

## 3. Alternatives Considered

- **Only the published maps.** Simpler, but then the atlas would only display a classification. Classifying the normals ourselves makes the rules inspectable and testable, and the comparison shows where climatologies disagree.
- **Classifying in the product builder (Python).** Rejected. The classification is domain computation, which belongs to the core (ADR 0002), here its provisional JavaScript stand-in, tested with the rest of it.
- **WorldClim or CHELSA normals.** Finer, but with no subsetting service: whole global files of hundreds of megabytes. WorldClim's terms also restrict redistribution. TerraClimate is CC0 and served by OPeNDAP.
- **Class maps as PNG images.** About thirty times smaller. But the browser and the tests would then decode them differently. One byte per cell is read identically by both, and a map is only loaded when chosen.
- **Downloading the whole archive.** It works, but takes fifty times the transfer for what is used. Range requests are standard, and the fetcher falls back to an error that says so when a server refuses them.

## 4. Consequences

- **The data library** gains:
  - the kind `opendap`, and members of remote archives;
  - the product kind `classes`;
  - climatologies;
  - a netCDF writer.
- **The atlas reads climate.**
  - Classes from 1901 to 2099, and how they change under four scenarios.
  - Normals by month.
  - A climate diagram for any place on land.
- **Tests:** the Earth Atlas has 93 checks.
- **Licences:** the maps are CC BY 4.0 and credited; TerraClimate is CC0, as stated in the Google Earth Engine catalogue entry.
- **Downloads:**
  - about 2.2 MB for the twelve maps (12 × 6.5 MB of products);
  - 56 MB for the normals (a 42 MB product).
