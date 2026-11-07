# 0079. Data in Space and Time: One Contract for Custom Code, and Works Found Again

- **Status:** Accepted, provisionally (2026-11-07). Amends ADR 0078.
- **Date:** 2026-11-07
- **Decided by:**
  - **The request.** Having seen ADR 0078's workbench, the user wrote: "original urban workbench's some features seems disappear? Like ABM, Simulation, Oceanography, complex analysis and visualization"; "become delicate, elegant, creative"; and "make plug-in interfaces that allow us visualize spatial temporal data, and utilize spatial temporal data, with custom codes."
  - **The design.** Claude tested every layer and tool of the older plug-ins before changing anything, then designed what follows.

## 1. Context

**The features had not gone.** Every layer and tool of Street Oceanography, Human flow, Landscape Dynamics and Event Spacetime was tried in the new shell, each opening with no error. They were hidden: ADR 0078 ordered everything by question, so that no plug-in's name appeared anywhere, and it dropped from the list whatever was not of the place open. Someone who knows a work by its name could not find it.

**Custom code existed only per subject.** Earth Atlas plug-ins (ADR 0029) are given fields, and Economy Atlas plug-ins (ADR 0046) are given networks. Nothing was given what happened where and when, which is the most common form of data a city or a society leaves.

## 2. Decision

### 2.1 Works Are Named Again

- The lens list has two orders, **by question** and **by work**, over one list. *By work* names each plug-in with a line saying what it is.
- Beneath a place's name stand the works written for it, each opening its section.
- Nothing is hidden for being elsewhere. A layer or tool of another place says where it is, and asking for it goes there and turns it on.
- The search under the place's title (Ctrl K) finds ways of seeing as well as places.

### 2.2 A Record in Space and Time

A record is a field on a part of the Earth times a line of time, in the sense of ADR 0053, in one of five forms, by what it is sampled on:

| Form | Base | Default mark |
|---|---|---|
| `events` | points of space × time | blooms |
| `tracks` | curves in space × time | threads |
| `cells` | a partition of space × steps of time | columns |
| `grid` | a raster × steps of time | on the ground |
| `flows` | pairs of places × time | threads, risen between them |

Space is degrees on WGS 84, and time is milliseconds since 1970 in UTC (CLAUDE.md §4.1). A date with no zone in a file is read as UTC and said to be. A date that reads two ways (03/04/2020) is refused. What cannot be placed is left out, counted and reported by reason, never dropped in silence.

A drawing keeps ADR 0053's composite: restriction (the moment and its window), reduction (the value chosen), mark.

### 2.3 One Contract in Three Languages

A script is **given** the records shown, its parameters and the context of the view, and **answers** records (drawn), tables, curves and notes. As in ADR 0029 and 0046, the same names hold in each language:

- **JavaScript** (`visualization/js/gis/spacetime.js`), in a worker of the page;
- **Python** (`api/python/science/spacetime.py`), as a process, on numpy;
- **C++** (`api/cpp/spacetime_plugin.hpp`), compiled on first use against the foundations and the domains, as a process.

`science.networks.Host` is generalised, not copied: `spacetime.Host` differs only in its module, its kind and what a run is given. A plug-in directory reaches the same line of time from any layer by `ctx.st.show(record)`.

### 2.4 The Studio

A plug-in, *The studio*, under a sixth lens, *Yours*:

- **Data:** a file dropped anywhere, which opens the studio to take it (CSV, TSV, JSON lines, GeoJSON or a record of the contract), its columns told by name and content and asked only if they cannot be, or a record of the data library;
- **Code:** a script opened, written or kept, with a form built from the parameters it declares, and its results drawn.

### 2.5 Marks that Live in Time

Time moves by a uniform on the graphics card, so that scrubbing a hundred thousand events costs the page nothing (`spacetime-marks.js`).

- **Blooms:** an event opens with a ring widening about it and fades, leaving dust if asked.
- **Threads:** fine strips of a stated width on the screen, drawn up to the moment, the front running between samples. Among a city's buildings a thread stays faintly seen where one stands in front of it.

A line of time (`src/shell/timeline.js`) stands above the clock: how much happens along it, a cursor, how long a mark stays, and optionally the Sun following it. A record stands among a city's buildings when it is about the open place at a city's scale, and on the Earth's canvas otherwise.

### 2.6 Who May Run What

A Python or C++ script runs with the authority of whoever started the server, and any page in the same browser can post to `127.0.0.1`. So the server writes and runs scripts only:

- when the request says it is the workbench's (a custom header, which another site's page cannot send without a permission this server never grants), with a JSON content type, and with an `Origin` matching its `Host`;
- at this computer's own address, unless started with `--token`;
- from one directory, under names that cannot leave it.

ADR 0029's rule stands: only files in the application's own directory are run. Each of these refusals has a check.

## 3. Results

- **Old features:** none was broken, across the 29 layers and tools of the four older plug-ins, each opened in the new shell.
- **A found fault in a table.** An empty latitude read as 0, since `Number('')` is 0 in JavaScript: the row stood on the equator and stretched its record's box to reach it, so that going to the record flew half the world away. An empty coordinate now leaves the row out, counted on the record's card, and a check holds it.
- **A found fault in compiling.** A C++ plug-in compiled against the workbench's 72 include directories did not finish in 291 s, since each header is looked for in each directory across WSL's slow mount. The host now passes only those a program's headers are found in, and it compiles in 23 s.
- **The examples are held to what they compute.** `hex_bins.cpp` counts every event once. `density.py`'s surface integrates to 395.8 of 400 events, the rest spread beyond its grid. `track_speed.js` reads a degree of the equator an hour as 111.2 km/h.

Pictures are in the [gallery's §17](../../docs/gallery.md).

## 4. Consequences

- A hover on a record of the library names no region: the catalogues' own region names are not under the data library's naming policy (ADR 0069), so they are not carried.
- `ctx.st` joins what a layer is given.
- A layer that dims the city for itself gives it back when turned off, and its dimming is not remembered between sessions. Event Spacetime's foam had been leaving every later picture of the city half seen.

## 5. Open

- A record is sent to a server script as JSON, which is slow beyond a few hundred thousand rows. A binary form is the next step.
- A record shown is not kept between sessions.
- Only one raster can lie on the ground at a time, as ADR 0078's drape allows.
- The editor is a plain text area: no highlighting or completion.
