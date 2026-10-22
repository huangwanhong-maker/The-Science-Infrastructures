# climate

*Status: draft. Bridge, height 4.* The climate system: the coupled ocean and atmosphere, its energy balance, its zones, its variability, and its astronomical forcing.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `zones` | Climatic zones and their classification: Köppen–Geiger after Peel et al. (2007) and Beck et al. (2018), with dry months after Gaussen | draft (JavaScript: `climate.js`) |
| `energy_balance` | Energy-balance models | planned |
| `coupled` | The coupled ocean–atmosphere system | planned |
| `variability` | Variability such as ENSO and the monsoons | planned |
| `orbital_forcing` | Insolation and Milankovitch cycles | planned |
| later | The cryosphere and the land surface | planned |

## 2. Dependencies

| Domain | Used for |
|---|---|
| meteorology | The atmosphere |
| oceanography | The ocean |
| geography | The land surface |
| astronomy | Orbital elements, insolation |
| air-sea | Coupling between ocean and atmosphere |

Foundations: physis; calx.

## 3. Code

`js/index.js` exports `climate.js`. The checks (fifteen constructed climates, each meeting one rule; invariance of the class under the group (S₆ × S₆) ≀ Z/2 acting on the months; every class reachable) are run by:

```
node domains/climate/tests/run-tests.mjs
```

The comparison of our classification with the published maps of Beck et al. reads the data library's products, and is the Earth Atlas's.

## 4. Notes

- [The Köppen–Geiger classification of climates, and the climate diagram](theory/koppen-geiger-classification.md)

## 5. Used By

The Earth Atlas.
