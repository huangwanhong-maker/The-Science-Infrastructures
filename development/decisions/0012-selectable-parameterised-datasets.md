# 0012. Selectable, Parameterised Datasets

- **Status:** Accepted
- **Date:** 2026-10-22
- **Decided by:** the user ("as a general library, we should allow to choose which data and how much of it"). Claude designed the mechanism.

## 1. Context

The first datasets were fixed files, such as the DE440 kernel and the HYG catalogue. The catalogues the project needs next are far larger than any one use of them:

- Gaia DR3 has 1.8 billion sources;
- galaxy surveys run to tens of millions of objects.

Archives serve such catalogues through queries, not files. Each application needs a different subset: the brightest stars, the stars within 100 pc, the galaxies within 11 Mpc, and so on.

## 2. Decision

1. **Two kinds of manifest.** A manifest in `data/registry/` is either:
   - `files`, which are downloaded and verified against recorded checksums;
   - or `tap`, which is a query to an archive through the IVOA Table Access Protocol.
2. **Parameters and profiles.** A `tap` manifest declares typed parameters, with bounds and descriptions, and named profiles. A retrieval chooses a profile and may override any parameter. Together these choose *which* part of a catalogue and *how much* of it.
3. **Cache layout.** Each choice gets its own folder:
   - `data/cache/<id>/<profile>/` for a named profile;
   - `data/cache/<id>/custom-<hash>/` when parameters are overridden.
4. **Provenance.** `retrieval.json` records:
   - the exact query;
   - the parameters;
   - the service;
   - the row count, the checksum and the time.
5. **The TAP client.** `data/catalogs/tap.py` is a client written with the standard library only. It always sends `MAXREC`, because services otherwise truncate results without notice: ESA Gaia returned 50 000 rows for a query of 200 000.
6. **Products and the index.** Products are built from cached datasets by `data/products/`. They are listed in `data/cache/products/index.json`, through which applications discover them, for example to let users choose a star catalogue.

## 3. Alternatives Considered

- **Download whole catalogues.** Rejected: Gaia alone is terabytes.
- **Archive-specific clients (astroquery).** Rejected for now: they add a large dependency for a protocol that a hundred lines can speak. They may be adopted later behind the same manifests.
- **Hard-coded queries in applications.** Rejected: queries would escape the registry, provenance would be lost, and subsets could not be shared between applications.

## 4. Consequences

- One mechanism serves ESA Gaia, CDS VizieR and SIMBAD, and will serve any other TAP archive, such as the NASA Exoplanet Archive or MAST.
- Query results may change if an archive revises a catalogue. The provenance record keeps the query and checksum of each retrieval.
- Archives other than TAP services (ERDDAP, OPeNDAP, STAC) will need further kinds of manifest, following the same pattern.
