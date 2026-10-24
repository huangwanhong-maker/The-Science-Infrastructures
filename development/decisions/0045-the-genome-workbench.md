# 0045. The Genome Workbench: the Genome Atlas Renamed, its Toolbox, the User's Files and a Shell Built from Registries

- **Status:** Accepted, provisionally; amends 0043
- **Date:** 2026-10-24
- **Decided by:** the user, on the Genome Atlas of ADR 0043: "The UI need largely improve... Many textbox, buttons, etc are near in raw style. And I think the lab will be highly large, and it is a lab, not atalas, it better rename. And in the tab, add an option, I don't know how to name it, maybe 'toolbox'. And it should can read from file(s). Anyways, the UI and interactive logics need largely revise, please take efforts to do this... be advanced feeling, be creative, and consider extenablity, as we may add many algorithms, many features in the future." Asked to choose a name, they chose "Genome Workbench". Claude designed the shell, the kit, the Toolbox and the file library, wrote them, had agents write the engine's tool operations and move the older views onto the kit, and wrote this record.

## 1. Context

The Genome Atlas of ADR 0043 had nine views, each built of raw elements: checkboxes, radios, bare number inputs and buttons styled one by one. Nothing in it could be given the user's own data; every sequence was a genome of the catalogue, a simulation, or letters pasted into a box. Its views were a fixed list in the page, so a new method meant a new page, written by hand. The user expects many algorithms and features to be added to it.

An atlas is a book of maps. What the user describes is a place of instruments, used on one's own material; "lab" is retired with "laboratory" (CLAUDE.md §1; ADR 0026), and of the names offered the user chose *workbench*.

## 2. Decision

1. **The name.** The application is the Genome Workbench, in `applications/genome-workbench/`, still on port 8890; its engine's namespace is `workbench` and its operations are registered by `BENCH_OP`. The engine's build stays in `data/cache/engines/genome/`. Records written before this one keep the old name and paths (CLAUDE.md §6.3).
2. **A shell built from registries** (`src/core/registry.js`, `src/shell.js`).
   - **What it registers.** Views, tools and commands are declared, and the shell builds itself from them:
     - the rail of views, in four groups (Explore, Analyse, Instruments, System);
     - a command palette (Ctrl K or /) that finds views, tools, genomes, the user's files and commands, and reads a range as a position, a word of DNA as a motif and a name as a gene;
     - a dataset switcher;
     - an activity tray of the requests and jobs in flight;
     - a status bar;
     - day and night themes;
     - files dropped anywhere on the page.
   - **The address.** The page's state is written there as `#<view>&...` (the older `view=` is still read), so that any view, tool and position can be linked to.
   - **What a new piece costs.** A new view, tool or command is one declaration. Nothing else changes.
3. **A kit and its tokens** (`src/ui/kit.js`, `src/ui/icons.js`, `styles.css`).
   - **The components:**
     - buttons and icon buttons, segmented choices, selects, number fields with steppers and units, sliders and switches;
     - a sequence editor that counts what it holds, reads files and the library, and accepts drops;
     - drop zones, cards, tabs, metrics, callouts, empty states, sortable tables, key–value lists, popovers, menus, dialogs and toasts.
   - **The tokens.** The colours are tokens: mulberry on blush paper by day, on plum by night. The type is Source Serif 4 throughout, with hairlines.
   - **The views** compose the kit on a common bench (a column of controls beside a stage of results), not raw elements.
4. **The Toolbox** (`src/views/toolbox.js`, `src/tools/`).
   - **A tool is a declaration:** `{id, title, group, icon, keywords, summary, refs, inputs, example, run, render}`.
     - Its `inputs` are a typed schema (`sequence`, `sequences`, `file`, `genome`, `number`, `choice`, `select`, `switch`, `text`, with groups and conditions), and `src/ui/form.js` builds the form from it.
     - Its `render` draws the result with the shared renderers of `src/ui/results.js`: alignments on a canvas, tables with their CSV, trees with their Newick, figures, and downloads of every result.
   - **The groups.** Tools are grouped by subject: sequences, alignment, annotation, indexes, reads, assembly, expression, phylogeny, populations, statistics.
   - **Runs** are kept per tool, so that an earlier result can be brought back.
   - **Examples are real sequences.** They are read on request from the RefSeq genomes the recipe fetched (the cytochrome b and COX1 genes of five primates, stretches of phiX174 and lambda), never typed in.
5. **The user's files** (`files.py`, an edge, with the server's `/api/files`).
   - **Storage.** Files are kept under `data/cache/genome-workbench/files/<hash>/`, by the hash of their content.
   - **Recognition.** Each is recognised by its content (FASTA, FASTQ, GenBank, VCF, Newick, SAM or text), unpacked if gzipped, and summarised.
   - **Genomes.** FASTA and GenBank files join the catalogue as genomes of the group "files".
   - **Into the engine.** A request names a file as `file:<id>` in a parameter called `files`, `*_file` or `*_files`; only the server turns the id into a path. The engine thus reads only files whose paths are passed in (CLAUDE.md §3.3), and the cluster sends such requests only to local workers.
   - **Limits.** An upload is limited to 4 GB locally and 64 MB on a public host.
6. **The engine's tool operations** (`tool.*`).
   - They take sequences as text or as a file and compute with the domains alone.
   - Where a tool's computation is an existing view's on other inputs, the core is shared, not copied.
   - A general method a tool needed and the domains lacked was written in its domain, with its check (§4).
7. **Two activities, two places.**
   - The analysis views keep their simulations with their truth, for learning and research.
   - The Toolbox serves use, on the user's own material.
   - Each view links to the tools that take the user's files for the same question.

## 3. Alternatives Considered

- **Keeping the name "Genome Atlas".** Rejected: the user found it wrong for what the application is becoming.
- **A second application for the tools.** Rejected: the tools and the views share the engine, the genomes, the files and the user's attention; one shell with groups scales better than two pages.
- **A front-end framework with a build step** (React, Svelte). Rejected: the infrastructure's web applications use ES modules without a build step (CLAUDE.md §7). The registries and the kit give the extensibility a framework would, in a few hundred lines that are ours.
- **Forms written by hand for each tool.** Rejected: dozens of tools are expected, and a schema makes each one a declaration, uniform in look and behaviour, and readable by the palette and the file library (`toolsFor(kind)`).
- **Reading files in the page only.** Rejected for anything large: reads and genomes of gigabytes belong on the server, beside the engine, which reads them by path. The page still reads small files directly into a sequence editor.

## 4. Consequences

- **Adding to it.**
  - A new algorithm is an engine operation and a tool declaration.
  - A new view is a module registered in `src/app.js`.
  - A new command is one `command({...})`.
- **Files.** The user's files stay on their computer, in the data cache, which is never committed.
- **The older views** were moved onto the kit; their behaviour, parameters and drawings are kept.
- **The documents.** CLAUDE.md §§1, 5, 7 and 9, the README, the applications' index, the domains' READMEs and the NCBI registry name the Genome Workbench.
- **Open.**
  - The palette could also search the contents of the user's files (genes in a GenBank file).
  - Tools could be chained into pipelines, a tool's result becoming a file of the library.
  - Plug-ins on the Earth Atlas's pattern (ADR 0027, 0029) could add tools without editing the application.
