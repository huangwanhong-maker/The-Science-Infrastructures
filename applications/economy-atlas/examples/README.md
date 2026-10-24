# examples: Data of One's Own, to Try

*Status: draft (2026-10-24).* Small files in the formats the Economy Atlas reads ([datasets.py](../datasets.py), ADR 0046), to try **+ your data** on the map and in spacetime. **Every number here is invented**, and so are the places and the currencies, so that no one mistakes them for a record.

| File | What it is | Where |
|---|---|---|
| `towns-places.csv` | eight invented towns on the North Sea and the Baltic: `id,name,latitude,longitude,group` | the map: give it with the flows |
| `towns-flows.csv` | invented flows of goods between them in 2024 and 2025: `from,to,value,year` | the map |
| `codes-flows.csv` | invented flows between six of the library's economies, named by their ISO alpha-3 codes, so that no places file is needed | the map |
| `rates.csv` | an invented record of three currencies, held at their policy rates and exchanged with the aurum (AUR), its base, month by month through 2025: `month,currency,rate,policy` | spacetime |

From the page: **+ your data**, drop the files (for the towns, both at once), give a title, **Add these data**. From the command line:

```
python applications/economy-atlas/datasets.py add flows "Towns" applications/economy-atlas/examples/towns-places.csv applications/economy-atlas/examples/towns-flows.csv
```
