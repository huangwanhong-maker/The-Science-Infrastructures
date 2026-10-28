# tiles: Terrain of the Earth and the Moon, Tile by Tile

*Status: draft (2026-10-28).* Terrain tiles in the geographic pyramid of the shared visualization layer's globes, built on demand from registered sources with the provenance of every sample ([ADR 0049](../../development/decisions/0049-the-spacetime-of-the-earth-and-the-moon.md)).

## 1. Purpose

A globe that a camera can bring from space to a mountainside needs heights at every scale, and no one source has them everywhere: a global grid of a few kilometres, regional models of tens of metres, nothing finer in most places, and on the Moon a laser altimeter's grids. This module cuts what the sources have into tiles of one pyramid, says of each sample where it came from and how good it is, and never pretends to more: a tile finer than its data says so, and detail made up below the data's resolution, when asked for, is labelled synthetic.

## 2. Concepts

- **The pyramid** (`pyramid.py`). Level z has 2^(z+1) × 2^z tiles in longitude and latitude, from 180° W and 90° N (NASA Trek's scheme, and the visualization layer's `globe/tiling.js`, held equal by the checks). A tile carries 65 × 65 samples at the corners of its 64 × 64 cells, so that neighbours share their edge, and a border of one sample beyond each edge for the normals there.
- **Sources** (`sources.py`). ETOPO 2022 at 5′ (the data library's field product) for the coarsest levels and the polar caps; Mapzen's Terrarium tiles on the AWS Open Data registry (SRTM, EU-DEM, NED, GMTED, ETOPO1 and others, each tile naming its sources) for the rest of the Earth, fetched on a pool of threads with kept connections; LOLA's gridded data records of the Moon at 4, 16, 64 and 128 pixels per degree, the finer two read by blocks with HTTP multi-range requests; EGM96 at 15′ for heights above the ellipsoid; Natural Earth's land (the product land-mask) to tell the sea from land below sea level.
- **The sea.** An Earth tile's heights are the surface one sees: over the sea its surface at the geoid, where the source is below sea level outside the land; the sea floor beneath is a second layer of the tile.
- **Provenance**, a byte per sample: 1 measured (a LOLA pixel with laser shots), 2 filled by the source's makers, 3 compiled (the source does not say), 4 finer than the data (interpolated here), 5 synthetic, 6 the sea's surface (the geoid).
- **Synthetic detail** (`synthetic.py`), only on request and only where a tile is finer than its data: octaves of value noise from half the data's resolution down to two cells, each of RMS amplitude g R^(1−H) λ^H and entering over two levels of the pyramid (none at two cells, the whole from eight), so that neighbouring levels differ by at most half an octave; the self-affine extrapolation of the region's gradient g at the coarse scale R (Mark and Aronson 1984 for the Earth, H = 0.8; Rosenburg et al. 2011 for the Moon, H = 0.9). A function of the point, so that tiles agree at their seams. It draws texture, not terrain.

## 3. The Encoding

```
'CPT1'          4 bytes
L               uint32, little-endian: the header's length
header          L bytes of JSON, padded with spaces to a multiple of 4
heights         m × m float32, rows north to south (m = 65 + 2)
provenance      m × m uint8
floor           m × m float32 after zero bytes to a multiple of 4, when the header's "layers" names it
```

The header gives the body, the tile, its bounds, the vertical datum, the sources and their resolutions, the least and greatest heights (and the floor's), the cell's size and the counts of each provenance code. Version 3 (the sea and its floor) is kept in `data/cache/tiles/terrain/<body>-v3/z/x/y.cpt`.

## 4. Content

| File | Content |
|---|---|
| `pyramid.py` | `TileKey`, `bounds`, `sample_lonlat`, `cell_degrees`, `cell_metres` |
| `sources.py` | `Geoid`, `LandMask`, `Etopo`, `Terrarium`, `Lola`; the provenance codes; the fetching pool |
| `terrain.py` | `encode`, `decode`; `TerrainService(cache)`: `tile(body, key)` built once and kept, `available()` |
| `synthetic.py` | `with_detail(service, body, tile)`, `detail`, `value_noise` |

## 5. Checks

`data/tests/check_tiles.py`: the encoding with and without its floor; the pyramid's closed forms and shared edges; the sources' rules (Terrarium's resolution for land and sea floor; the land mask at the Dead Sea, the Caspian, the Qattara Depression, Flevoland, Death Valley and four seas); a built tile's sea at the geoid and its floor beneath; the noise's unit RMS, determinism and seams, the octaves' weights, and their amplitude by the self-affine law.

## 6. Used By

Cosmopoiesis's server (`/tiles/terrain/<body>/<z>/<x>/<y>`), and through it the visualization layer's globes.
