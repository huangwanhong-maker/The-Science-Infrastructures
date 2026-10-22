# observation: Instruments and Field Work

*Status: planned.* Edge.

## 1. Purpose

This directory covers practical astronomy and oceanography with real instruments: controlling the instruments, acquiring measurements, and keeping observation logs. It handles only the contact with the physical world. The algorithms that reduce raw measurements, such as calibration, astrometry and sight reduction, belong in `core/domains/`.

## 2. Planned Modules

| Module | Scope |
|---|---|
| `instruments` | Control of telescopes, mounts and cameras through INDI or ASCOM |
| `field` | Sextant sights, tide-pole readings, observation logs |
| `sensors` | Tide, weather and ocean sensors, including home-built ones |
