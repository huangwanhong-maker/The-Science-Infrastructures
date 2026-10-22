# observation: Instruments and Field Work

*Status: draft* (`sensors`); the other modules are planned. Edge.

## 1. Purpose

This directory covers practical astronomy and oceanography with real instruments: controlling the instruments, acquiring measurements, and keeping observation logs. It handles only the contact with the physical world. The algorithms that reduce raw measurements, such as calibration, astrometry and sight reduction, belong in `domains/`.

## 2. Modules

| Module | Scope | Status |
|---|---|---|
| `instruments` | Control of telescopes, mounts and cameras through INDI or ASCOM | planned |
| `field` | Sextant sights, tide-pole readings, observation logs | planned |
| [`sensors`](sensors/) | Tide, weather and ocean sensors, including home-built ones: a store of their readings in SI units ([store.py](sensors/store.py)), ingestion from files and loggers ([ingest.py](sensors/ingest.py)), a simulated sensor as the template for a device ([simulate.py](sensors/simulate.py)), and an observations product the Earth Atlas draws | draft |

## 3. Dependencies

`sensors` writes its products through the data library's writer of observations (`data/products/observation_set.py`) and index (`data/products/product_index.py`), and keeps its readings in `data/cache/observations/`. It uses the Python standard library only. Nothing in `foundations/` or `domains/` depends on it.
