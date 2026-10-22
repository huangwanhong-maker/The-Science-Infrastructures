# 0016. Maps of the Bodies, IAU Rotation, and Imagery Shared Across Applications

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Amended by [ADR 0020](0020-detail-of-bodies-at-close-range.md) (detail at close range).
- **Date:** 2026-10-22
- **Decided by:** the user asked for finer textures in the 3D views of both atlases. Claude chose the design below.

## 1. Context

The Solar System view drew each body as a sphere of one colour. A map on the sphere teaches more than a colour does, but only if the map is turned correctly. A map of Mars that does not show the meridian that actually faces the Earth at that instant is decoration, not astronomy. The orientation of each body is therefore part of the computation and must be verified like its position.

The Earth Atlas had already introduced imagery products (ADR 0015). The Blue Marble and Black Marble it uses are wanted in the Celestial Atlas too.

## 2. Decision

1. **The IAU rotation models.**
   - The Sun, the Moon and the planets are oriented by the models of the IAU Working Group on Cartographic Coordinates and Rotational Elements (Archinal et al. 2018): the right ascension α₀ and declination δ₀ of the north pole, and the angle W of the prime meridian, all in the ICRF.
   - The Moon's model includes its thirteen periodic terms; Mercury and Mars include theirs.
   - The rotation from the ICRF to body-fixed axes is B = R_z(W) R_x(90° − δ₀) R_z(90° + α₀).
   - The code is `src/core/rotation.js`, a provisional stand-in for `core/domains/astronomy`.
2. **Verification against JPL Horizons.**
   - For every body at five epochs from 2000 to 2036, the sub-Earth point (the longitude and latitude on the body of the direction to the Earth's centre) is compared with Horizons.
   - The fixture is fetched by `tests/fetch_orientation.py` and stored with its source and retrieval date.
   - Every body agrees within 0.02°; the worst is Mercury, at 0.011°.
3. **Maps as imagery products.**
   - A new dataset, `solar-system-maps`, lists maps of the Sun, the planets, the Moon, Saturn's rings and the Milky Way. Each file names its body.
   - `data/products/body_maps.py` makes each map an imagery product with one of three new roles: `body`, `rings` and `sky`.
   - The Milky Way map's orientation is checked, not assumed. The brightest pixels near twelve bright stars of known galactic position must lie where the stated orientation puts them, or the build fails.
4. **One mapping convention.**
   - A body map is in east longitude, 0 at the centre, north at the top.
   - The local frame of a textured sphere is related to body-fixed axes by (x, y, z)_body = (x, −z, y)_local, and to galactic axes the same way for the sky.
   - These constants live in `src/render/maps.js`, beside the rule that derives them.
5. **Imagery shared across applications.**
   - The Celestial Atlas reads the Blue Marble and Black Marble products that the Earth Atlas uses. Its Earth shows the lights of cities on its night side.
   - The product index, not either application, owns the images. An application uses whatever imagery the index lists and degrades gracefully when it lists none.

## 3. Alternatives Considered

- **Orientation from SPICE PCK kernels.** The PCK files hold the same IAU models. Reading them would add a parser and a download for no gain in accuracy, and the coefficients are short enough to cite in the code.
- **Maps from each mission's archive** (LRO for the Moon, MGS for Mars, Cassini for Saturn). Truer and finer, but in many projections, sizes and licences. Solar System Scope's set is uniform, CC BY 4.0, and documents its sources. A body's map can be replaced by a finer one later without any change to the code.
- **Rings from geometry alone.** The rings were previously drawn from their radii as a plain band. The radial map shows the Cassini division and the C ring's transparency, at no cost.

## 4. Consequences

- **Bodies show their true faces.** Each map is turned to the instant shown, to within the accuracy of the IAU models. Saturn's rings lie in its equatorial plane.
- **The Celestial Atlas has 57 checks**, ten more than before: one sub-Earth point per body.
- **Licences.** Solar System Scope textures are CC BY 4.0 and credited in the atlas. The Blue Marble and Black Marble are NASA's, acknowledged as GIBS requests.
- **Sizes.** Maps are resampled to 4096 pixels wide, and the sky to 8192. About 17 MB is loaded after the first frame, and the atlas works without them.
