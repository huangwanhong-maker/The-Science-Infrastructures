# 0021. The Sky from a Place: Choosing Where to Stand, and Links Between the Atlases

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review.
- **Date:** 2026-10-22
- **Decided by:** the user suggested that clicking somewhere on the Earth should show the sky there. Claude chose the design below.

## 1. Context

The Celestial Atlas's Sky Chart was drawn for a place chosen from a short list or typed as coordinates. The Earth Atlas knows places but not their sky, and the Celestial Atlas knows the sky but not the places. The two applications are separate, and neither may depend on the other's code (CLAUDE.md §3.3).

## 2. Decision

1. **Choosing a place by clicking the Earth, in two ways.**
   - In the Sky Chart, a small globe shows the Earth by day (Blue Marble) and its lights by night (Black Marble), lit by the Sun at the instant shown, in orthographic projection. Drag it to turn it; click it to stand there.
   - In the Solar System view, a click on the focused Earth chooses the place under the pointer, marks it, and offers the sky from there.
2. **One frame for the Earth.**
   - In the Solar System view, the Earth is turned by Greenwich mean sidereal time and precession: the frame in which the Sky Chart places its observer.
   - A click on the map is therefore the place whose sky is drawn. The IAU model of the Earth's rotation, used for the other bodies, differs from this frame by 0.05° to 0.26° of longitude between 2000 and 2100.
   - Clicked points are converted from geocentric to geodetic latitude on WGS 84 (`src/core/geodesy.js`).
3. **Places named from the shared geography product.**
   - The Celestial Atlas reads the data library's geography product (Natural Earth's populated places) when first needed. It names a chosen place after the nearest town within 120 km.
   - The place takes that town's time zone within 500 km, or nautical time at sea. The Sky Chart shows local civil time there.
   - Products are shared between applications, as imagery already is (ADR 0016).
4. **Links between the atlases, by their addresses.**
   - The Celestial Atlas reads a view, a place and an instant from its address (`#view=sky&lat=…&lon=…&t=…&name=…`).
   - The Earth Atlas's place card links to the sky at its place and time, and the Sky Chart links back to the weather there.
   - Each application only builds an address for the other, so neither depends on the other's code. On the local servers the address names the other's port; on a site holding the repository, its path.

## 3. Alternatives Considered

- **Drawing a sky chart inside the Earth Atlas.** This would duplicate the Celestial Atlas's ephemeris, star catalogue and chart. Sharing that code belongs in `visualization/`, a move still awaiting the user's agreement.
- **A 3D globe for the picker.** The Solar System view already has one. In the panel, a 2D orthographic globe is light, needs no second WebGL context, and is tested by its projection's round trip.
- **Reverse geocoding by a web service.** It would send the user's chosen places to a third party. The Natural Earth places are local and sufficient at this scale.

## 4. Consequences

- **Two ways to choose where to stand,** and a path from the weather at a place to its sky and back.
- **Tests:** the Celestial Atlas has 66 checks, adding:
  - the latitude conversions and their greatest difference (0.1924° at 45°);
  - the great-circle distance;
  - the globe's projection, both ways.
- **Data:** the geography product (1.6 MB) is read by the Celestial Atlas only when a place is first chosen by clicking, or when a linked place needs its time zone.
