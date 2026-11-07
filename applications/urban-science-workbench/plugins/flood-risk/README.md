# Flood risk

*A city's flood risk read from its own assessment, with the hazards named for any place.* Status: **draft**. Activity: research. [ADR 0078](../../../../development/decisions/0078-the-whole-earth-in-the-urban-workbench.md).

## 1. Purpose

Freetown is the one city with an open assessment of this depth: 1 m flood maps, every building modelled, and its loss tables. The plug-in draws them on the city, recomputes what it can, and says how far they agree.

## 2. Layers

- **How deep each flood reaches:** a map of the assessment by peril (river, surface water) and return period (20 to 1500 years), taken at about 4 m, on the ground.
- **The buildings it reaches:** those of the 153,435 modelled buildings standing in over 15 cm, on H3 cells by the value reached. What stands in the water, not the loss.
- **What each ward stands to lose:** the assessment's annual loss by ward, per square kilometre.
- **Where the Council records floods:** the City Council's own places and areas, a check of the model by what the city has seen.

## 3. Tools

- **What a flood costs the city:** the annual loss and exceedance curves, its figures beside what is recomputed here.
- **The hazards named for this place:** ThinkHazard!'s levels for an administrative unit, anywhere.

## 4. Engine

`engine/main.cpp`, on urban-science's `exposure`, tyche's `exceedance` and physis's H3: `summary`, `exposed`. It reads the product `data/products/freetown_flood.py` writes.

## 5. What Was Found

- The catalogue's length is not published in the tables. The sum of the event losses over the annual loss is 10,000.0, so it is 10,000 years.
- Recomputed from its 1,943 events, the annual loss is $808,904.91, the assessment's own figure to the cent.
- The published occurrence curve is the events in order, to 1 part in 10⁷, at every rank a return period falls on.
- From the dozen published points alone, the area within them is $762,246, 5.8 % below; holding the last loss beyond 1 in 1,000 years adds $82,630, ending 4.4 % above. The tail is an assumption that decides the sign of the error.
- A 1-in-100-year river flood puts 7,293 buildings, $384 million, in over 15 cm.

## 6. Limits

The values are the assessment's replacement values, in its USD. No vulnerability is applied here: the losses shown are the assessment's own. Surface water and river floods are drawn separately, as it publishes them.
