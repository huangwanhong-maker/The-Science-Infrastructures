# assets: Static Resources

## 1. Purpose

This directory holds small static files used across the project: images, fonts, icons, style tokens and small reference tables.

## 2. Rules

- Large or external datasets do not belong here. Register them in `data/` and fetch them through it.
- For every asset the project did not create, record its source and its licence.

## 3. Contents

| Path | Content | Status |
|---|---|---|
| [`gallery/`](gallery/) | Screenshots of the Earth Atlas and the Celestial Atlas, JPEG at 1600 × 1000, shown in [docs/gallery.md](../docs/gallery.md). Named `<application>-<number>-<subject>.jpg`. | draft |

The screenshots are the project's own, under its MIT licence, except `earth-15`, `earth-16` and `earth-17`. These show forecasts of Pangu-Weather, so they are under CC BY-NC-SA 4.0, the licence of its parameters. The data drawn in them keep their own licences, listed in the gallery's credits. To make them again, run the applications locally and open the addresses given in [docs/gallery.md §4](../docs/gallery.md).
