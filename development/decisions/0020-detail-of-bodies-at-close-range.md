# 0020. Detail of the Moon, Mars and Mercury at Close Range

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Amends [ADR 0016](0016-maps-of-bodies-and-iau-rotation.md).
- **Date:** 2026-10-22
- **Decided by:** the user asked for finer textures in the 3D views of both atlases. Claude chose the design below.

## 1. Context

ADR 0016 gave each body a global map of 4096 pixels: about 2.7 km to a pixel on the Moon and 5 km on Mars. Focused on a body, the view can come much closer, and the maps blur.

NASA's Trek portals serve the mission mosaics as equirectangular tiles, and allow browsers to fetch them (CORS):
- the Moon from LRO's wide-angle camera, about 100 m;
- Mars from Viking, 232 m;
- Mercury from MESSENGER, 665 m.

## 2. Decision

1. **The services are declared by the products.**
   - In `solar-system-maps`, the Moon's, Mars's and Mercury's entries name their tile service (`live`: a WMTS template, its tile size and deepest level). `body_maps.py` passes it into the products.
   - As in ADR 0015, an application uses only a service that a product declares.
2. **Detail where the user looks.**
   - When the focused body fills a third of the view or more, rays from a grid on the screen find the region in view, in body-fixed longitude and latitude. They go through the same mapping as the global map, so the detail follows the body's rotation.
   - The level whose pixels match the screen's is chosen, at most 48 tiles (`src/core/tiles.js`, tested).
   - The tiles are composed into one texture and drawn over the global map by the body's shader, fading at the edges (`src/render/detail.js`).
3. **Closer approach.**
   - When focused, the camera may come to 1.03 radii from the body's centre.
   - Spheres are finer (192 × 96), so the limb stays round.
4. **Two fixes to the local server.**
   - Products are revalidated rather than kept for a day, so a rebuilt product is seen.
   - Its connection queue is deeper, as the Earth Atlas's already was.

## 3. Alternatives Considered

- **Finer global maps.** An 8k or 16k global map would cost hundreds of megabytes of graphics memory, for detail seen only up close.
- **A tiled sphere (a quadtree of patches).** It is the general solution, used by globe engines, but it means a new geometry and its seams. One detail region on the existing sphere is simpler and serves a view that looks at one body.
- **All bodies.** Venus's clouds, the gas giants and the Sun have no finer mosaics worth fetching. Trek also serves Vesta, Ceres and others, which the atlas does not show yet.

## 4. Consequences

- **Detail.** Close to the Moon, its craters show at about 1 km to a screen pixel, and finer as the view approaches. Mars and Mercury show their own mosaics.
- **Tests:** the Celestial Atlas has 61 checks.
- **Tone.** The mosaics differ a little in tone from the Solar System Scope maps. The edge fade softens the difference but does not remove it.
