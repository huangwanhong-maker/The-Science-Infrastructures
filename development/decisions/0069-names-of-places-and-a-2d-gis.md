# 0069. The Names of Places, Areas in Dispute, and a 2D GIS

- **Status:** Accepted, provisionally (the policy's names and sources are the user's to confirm)
- **Date:** 2026-11-04
- **Amends:** [ADR 0014](0014-earth-atlas-separated-data-and-typed-layers.md) (the Earth Atlas's names), [ADR 0044](0044-value-in-space-and-time.md), [ADR 0046](0046-the-economy-atlas-data-and-plugins.md), [ADR 0068](0068-the-economy-atlas-as-a-scholarly-instrument.md)
- **Decided by:** the user: "some sensitive things should not appear on the map, like "TAIWAN". And please add 2D GIS module if lacks", and then: "Also other current sensitive naming in the world should avoid." Claude chose the policy, its sources and its treatment of areas in dispute, and wrote this record.

## 1. Context

The infrastructure named places as Natural Earth names them. Natural Earth names Taiwan as a country, and in Chinese and Japanese as the Republic of China. It lists Northern Cyprus and Somaliland as countries. Its polygons follow de facto control: Crimea is drawn in Russia, the Golan Heights in Israel, the southern Kurils in Russia, Aksai Chin in China, and Arunachal Pradesh and Jammu and Kashmir in India. It lists Jerusalem as Israel's capital and Laayoune as Morocco's. The Economy Atlas's map labelled its largest economies by these names, its cards and tooltips used them, and a choropleth coloured each area in dispute with its administering country's data. The gazetteer, which Cosmopoiesis draws on its globe, listed these places in the same way. The Earth Atlas labelled the seas whose names are in dispute.

The Economy Atlas's map was also its own code: the projection's paths and arcs, zoom and pan, the pointer's readout, labels placed without overlap, and the figure composed for publication. Any other two-dimensional map would have to write all of it again. The visualization layer had a GIS of cities in three dimensions (`gis`) and virtual globes (`globe`), but no map of the sphere in the plane.

## 2. Decision

### 2.1 A Policy of Names, Kept as Data

The policy is a file of the data library, `data/products/naming.json`, read by `data/products/naming.py`. It states its rules:

- **The ordinary rule.** Places are named by Natural Earth's long names, except where a name is in dispute or out of date.
- **Names in dispute or out of date.** These take the short name of the UN's M49, the usage of international statistics.
- **Places M49 omits or names in dispute.** These are named as the organisations whose statistics the infrastructure draws name them: the WTO, the OECD, the World Bank and the IMF.
- **Former names.** Names current before are kept as aliases, so that data of one's own written with them still resolve.

Each entry gives its source:

| Code | Name | Labelled on a map | Source |
|---|---|---|---|
| TWN | Chinese Taipei (中华台北, 中華台北) | no | WTO and OECD |
| HKG, MAC | Hong Kong, China; Macao, China | yes | WTO and OECD; M49 |
| PSE, PSX | West Bank and Gaza | no | World Bank and IMF |
| XKX, KOS | Kosovo (UNSCR 1244) | no | Security Council resolution 1244 |
| ESH, SAH | Western Sahara | no | M49 |
| FLK | Falkland Islands (Malvinas) | no | M49 |
| CYN, SOL, KAS, IOT | no name | no | not recognised by the UN, of indeterminate status, or in transfer |
| TUR, CZE, SWZ, CPV, CIV, MKD, MMR, TLS | Türkiye, Czechia, Eswatini, Cabo Verde, Côte d'Ivoire, North Macedonia, Myanmar, Timor-Leste | yes | M49 |

Every product that names countries applies the policy:

- the trade flows;
- the geography: its countries' names in every language, and whether a map may label them;
- the gazetteer;
- the regions;
- the Economy Atlas's tables and datasets, where an alias resolves to its code.

The Economy Atlas's server gives the policy to its page at `/api/names`.

### 2.2 What a Map Does Not Draw

- **Territories.** A map draws no name for a territory whose status is in dispute, and no name at all for one the policy leaves unnamed. Its data, where there are data, carry the policy's name in cards and tables.
- **Seas.** A map draws no name for the five seas whose names are in dispute: the Sea of Japan, the Persian Gulf, the Gulf of Mexico, the South China Sea and the Korea Strait.
- **Areas in dispute.** Natural Earth's 28 areas in dispute or breakaway (`ne_50m_admin_0_breakaway_disputed_areas`, added to the registry) are drawn over the countries. They are hatched on the paper, in no claimant's colour, and unnamed. Whatever colours the countries does not colour them, and a tooltip says only that the area is in dispute. Neither the countries' polygons nor their data are cut; the overlay is drawn over them.
- **Places.** A place in an area in dispute is given no country. No place there, or in a territory the map does not label, is listed as a national capital. Two places are named in the policy itself: Jerusalem has no country and is not listed as a capital (General Assembly resolution 181, Security Council resolution 478). Nicosia is the capital of Cyprus, although Natural Earth's point for it lies north of the buffer zone.

### 2.3 The Names in Prose

- **Taiwan as an economy.** In the theory note of the economics domain and in the gallery's captions, the economy is Chinese Taipei.
- **Seas in prose.** Where captions and notes named a sea in dispute, they now describe it: "the sea between Japan and the mainland", "the Gulf between Arabia and Iran".
- **Geographic names.** The island of Taiwan and the Taiwan Strait keep their names, which are geographic.

### 2.4 A 2D GIS in the Visualization Layer

The visualization layer has a module of maps of the sphere in the plane, `visualization/js/map`:

- `projection.js`: the geography domain's projections in the plane of a drawing, with outlines, graticules, polygons, polylines, and great-circle arcs cut at the map's edge;
- `view.js`: zoom about a point, fitting to points and centring;
- `labels.js`: names placed by priority without overlap, after the naming policy;
- `widget.js`: an interactive map, with its base, the areas in dispute, an overlay for an application's layers, zoom, pan, picking, and the readout of projection, place and zoom;
- `figure.js`: a map composed as a figure for publication, as SVG with its fonts or as PNG.

The pure parts are checked in Node by the visualization suite. The Economy Atlas's map now stands on the module, and what remains in its `map.js` is what names and draws economies.

## 3. Alternatives Considered

- **Natural Earth's editions by point of view** (its Chinese, Ukrainian and other editions). Each takes one country's view of its own borders. The policy instead takes the usage of international statistics, and leaves the areas in dispute in no one's colour.
- **Cutting the areas in dispute out of the countries' polygons.** This changes the countries' geometry and their areas. An overlay leaves the data as they are.
- **"Taiwan Province of China" or "Taiwan, China".** Chinese Taipei is the usage of the WTO and the OECD, from whom the infrastructure's trade data come. The other names resolve as aliases.
- **Dual names for seas**, such as "Sea of Japan (East Sea)". These still name the dispute; drawing no label is quieter.

## 4. Consequences

- **Rebuilding.** The policy shows only in products built after it: the trade flows, the geography, the gazetteer and the regions. A recipe does not judge a product stale when its sources are unchanged, so these are rebuilt by hand once. The Earth Atlas's and the Economy Atlas's recipes now fetch the areas in dispute. Cosmopoiesis reads the rebuilt gazetteer; its own recipe, in its own repository, is unchanged.
- **The Earth Atlas.** It follows the policy's names and labels. It does not yet draw the areas in dispute as a layer, since it colours no country by data; that layer is planned.
- **Changing a name.** A name the user wants otherwise is changed in one file, and every product and map follows when rebuilt.
