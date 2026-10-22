# 0010. Domain-Organised Fields and Applications

- **Status:** Accepted
- **Date:** 2026-10-22
- **Decided by:** the user, who asked for:
  - `core/domains/` to hold "the domains we investigate integratively";
  - `applications/` to have the same folders;
  - the field to be renamed "meteorology";
  - geography to cover all of its branches.

  Claude designed the placement.

## 1. Context

- **The fields were invisible.** `core/domains/` listed its fields and bridges only in a README table, because directories were created only together with their first content ([0002](0002-layered-architecture-with-pure-core.md)). The fields under investigation therefore did not appear in the tree.
- **Geography was missing.** The laboratory's name, confirmed on 2026-10-22, includes geography, but geography was not one of the fields.
- **One field was named differently.** The atmospheric field was called `atmosphere`, whereas the other fields carried the names of their disciplines.
- **Applications were split by use, not by subject.** `applications/` was divided into `labs/` and `tools/`. That division says how an application is used, not what it is about.

## 2. Decision

1. **Fields.** `core/domains/` holds the fields investigated integratively: astronomy, geography, oceanography and meteorology. Geodesy is a supporting field.
2. **Geography.** Geography is taken in all its branches:
   - physical geography;
   - cartography and GIS;
   - remote sensing;
   - human geography;
   - other branches, as the work requires.
3. **Meteorology.** The field `atmosphere` is renamed `meteorology`, in the broad sense of atmospheric science. This amends [0007](0007-atmosphere-and-ocean-as-coupled-system.md).
4. **Bridges.**
   - The bridges stay as they were: tides, navigation, air_sea and climate.
   - Geography now joins navigation and climate.
   - Candidate new bridges are `coasts` and `hydrology`.
5. **Directories from the start.** Each field and bridge has a directory with a README from the start. This is an exception to the on-demand rule of 0002, because these directories define the shape of the project. Modules inside them are still created on demand.
6. **Applications mirror the domains.** `applications/` mirrors `core/domains/`, with one directory per field and bridge.
   - An application lives under the field or bridge it mainly serves.
   - An application that serves several fields belongs under the bridge that joins them.
   - Whether an application is a lab or a tool is recorded in its README, not by its directory.

## 3. Alternatives Considered

- **Keep fields as README rows only.** Rejected: the fields under investigation should be visible as the structure of the project.
- **Keep the name `atmosphere`.** Rejected by the user in favour of a discipline name.
- **Keep `applications/labs/` and `applications/tools/`.** Rejected. Organising by field keeps each application next to the science it serves. Being a lab or a tool is a property of the application, not a place.
- **Fold geodesy into geography.** Rejected. Geodesy is a distinct science, covering the gravity field and the Earth's rotation, and oceanography and meteorology also rest on it.

## 4. Consequences

- The scope statements, layouts, placement guide and theory subjects are updated.
- New datasets are planned for geography: elevation models, cartographic vector data, satellite imagery, and grids of population and land use.
- Human geography brings social data into the laboratory. Each dataset keeps its own licence; OpenStreetMap, for example, uses ODbL.
