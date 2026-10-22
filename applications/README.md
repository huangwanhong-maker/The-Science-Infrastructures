# applications: Laboratories and Tools

*Status: draft.* This directory is an edge at the top of the architecture. It may depend on everything, and nothing depends on it.

## 1. Purpose

Applications are complete programs that people use. They come in two kinds:

- **Laboratories (labs)**, for learning and research.
- **Tools**, for practical work.

## 2. Organisation

`applications/` mirrors [`core/domains/`](../core/domains/), with one directory per field and per bridge ([ADR 0010](../development/decisions/0010-domain-organised-fields-and-applications.md)).

- **Placement.** An application lives under the field or bridge it mainly serves.
- **Applications that serve several fields** belong under the bridge that joins those fields. If no bridge fits, a new bridge is probably needed in `core/domains/`.
- **Lab or tool.** Each application records whether it is a lab or a tool in its README. The directory structure does not encode this.
- **Dependencies.** Each application declares what it depends on.

## 3. Index

| Directory | Planned applications |
|---|---|
| [`astronomy/`](astronomy/) | **Celestial Atlas (tool, working)**; orbital sandbox (lab); ephemeris calculator (tool) |
| [`geography/`](geography/) | **Earth Atlas (tool, working)**; projection lab, terrain and watershed lab, remote-sensing lab (labs); land-use atlas (tool) |
| [`oceanography/`](oceanography/) | Rotating tank (lab); ocean-profile explorer (tool) |
| [`meteorology/`](meteorology/) | Weather-map lab and atmospheric-sounding lab (labs) |
| [`geodesy/`](geodesy/) | Foucault pendulum (lab) |
| [`tides/`](tides/) | Tide-predicting machine (lab); tide predictor (tool) |
| [`navigation/`](navigation/) | Sight-reduction trainer (lab); route planner (tool) |
| [`air_sea/`](air_sea/) | Air–sea exchange lab (lab) |
| [`climate/`](climate/) | Energy-balance lab and Milankovitch lab (labs) |
