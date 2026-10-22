---
title: The Subsolar Point, the Terminator and Twilight
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/sun.js
  - applications/geography/earth-atlas/src/render/surface.js
references: [AstronomicalAlmanac2019, Meeus1998]
---

# The Subsolar Point, the Terminator and Twilight

## 1. The Subsolar Point

**Definition 1.** The *subsolar point* at an instant is the place on the Earth where the Sun is at the zenith.

**Proposition 1.** Let α and δ be the Sun's apparent geocentric right ascension and declination, and θ the Greenwich apparent sidereal time, all in degrees. The subsolar point has latitude δ and east longitude

  λ_s = α − θ  (mod 360°).  (1)

*Proof.* The local hour angle of the Sun at east longitude λ is H = θ + λ − α. The Sun is on the meridian (H = 0) at λ = α − θ. Its altitude h satisfies

  sin h = sin φ sin δ + cos φ cos δ cos H,  (2)

which equals 1 exactly when φ = δ and H = 0. ∎

*Remark.* JPL Horizons also offers a "sub-solar point" for an observer placed at the Sun. That point is the place below the observer at the time the light left the Earth, so it is later by the light-time of about 8.3 min, during which the Earth turns about 2.08°. It is a different quantity, and the tests use (1) instead.

## 2. The Sun's Position, to 0.01°

The Astronomical Almanac (section C) gives, for n days from J2000.0:

- L = 280.460° + 0.9856474° n, the mean longitude, corrected for aberration;
- g = 357.528° + 0.9856003° n, the mean anomaly;
- λ = L + 1.915° sin g + 0.020° sin 2g, the ecliptic longitude;
- ε = 23.439° − 0.0000004° n, the obliquity of the ecliptic;
- α = atan2(cos ε sin λ, cos λ) and δ = asin(sin ε sin λ). (3)

These are accurate to 0.01° from 1950 to 2050.

The atlas takes Greenwich mean sidereal time from Meeus (1998), eq. 12.4, in place of θ:
- the equation of the equinoxes is at most 1.2 s, or 0.005°;
- UT1 − UTC is at most 0.9 s, or 0.004°.

*Verification.* Against Horizons at 20 instants from 1955 to 2045, the subsolar point from (1) and (3) is within 0.007° in longitude and 0.004° in latitude, about 0.8 km.

## 3. The Terminator and Twilight

By (2), the Sun's altitude at the point **p** of the unit sphere is h = asin(**p** · **s**), where **s** is the unit vector toward the subsolar point. Each class of daylight is a zone between two small circles about **s**:

| Altitude h | Class | Reason |
|---|---|---|
| h ≥ −0.833° | day | The upper limb is above the horizon: refraction lifts it by 34′, and the Sun's radius is 16′ |
| −6° ≤ h < −0.833° | civil twilight | Outdoor work without lights |
| −12° ≤ h < −6° | nautical twilight | The horizon at sea is still visible |
| −18° ≤ h < −12° | astronomical twilight | The sky is not yet fully dark |
| h < −18° | night | |

The atlas darkens the map in four steps at these altitudes, so that the twilight zones are seen as bands. Data drawn over the map are darkened half as much, so that they stay readable at night.

The tropics (latitudes ±ε) and the polar circles (latitudes ±(90° − ε)) are drawn on the graticule:
- the subsolar point stays between the tropics;
- beyond the polar circles the Sun can stay above or below the horizon for a whole day.

*Caveat.* The terminator is drawn for a geocentric observer on a sphere. Refraction is included only through the −0.833° of sunrise, and the Earth's flattening and the observer's height are neglected. The error is a few kilometres.
