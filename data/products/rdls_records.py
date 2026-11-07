"""The Freetown Flood Risk Assessment as a record of the Risk Data Library Standard (RDLS 1.0, GFDRR), the catalogue entry by which others can find
what the workbench draws. The schema is fetched by the manifest `rdls-schema`, not kept here (CC BY-SA 4.0).

The record is of the dataset as a whole: its required fields, its coverage, its resources and who made it. `validate` holds a value against the
parts of JSON Schema the record needs: $ref, required, properties, items, type, enum, minItems, minLength. It does not evaluate the schema's
conditional rules (allOf, anyOf), so that passing here is necessary and not sufficient for conformance; a full validator may be run over the
file where one is installed.

    python data/products/rdls_records.py

writes data/cache/products/rdls-freetown.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

CACHE = Path(__file__).resolve().parent.parent / "cache"
TYPES = {"string": str, "array": list, "object": dict, "number": (int, float), "integer": int, "boolean": bool}


def validate(value, schema: dict, root: dict | None = None, at: str = "$") -> list[str]:
    """What is wrong with a value against a schema, as far as this reader goes."""
    root = root or schema
    if "$ref" in schema:
        target = root
        for part in schema["$ref"].lstrip("#/").split("/"):
            target = target[part]
        return validate(value, {**target, **{k: v for k, v in schema.items() if k != "$ref"}}, root, at)
    wrong: list[str] = []
    kind = schema.get("type")
    if isinstance(kind, str) and kind in TYPES and not isinstance(value, TYPES[kind]):
        return [f"{at}: not of type {kind}"]
    if "enum" in schema and value not in schema["enum"]:
        wrong.append(f"{at}: {value!r} is not in its list")
    if isinstance(value, str) and len(value) < schema.get("minLength", 0):
        wrong.append(f"{at}: too short")
    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            wrong.append(f"{at}: too few items")
        if "items" in schema:
            for k, item in enumerate(value):
                wrong += validate(item, schema["items"], root, f"{at}[{k}]")
    if isinstance(value, dict):
        for key in schema.get("required", []):
            if key not in value:
                wrong.append(f"{at}: {key} is missing")
        for key, sub in schema.get("properties", {}).items():
            if key in value:
                wrong += validate(value[key], sub, root, f"{at}.{key}")
    return wrong


def freetown() -> dict:
    base = "https://datacatalog.worldbank.org/search/dataset/0066934"
    bank = {"name": "The World Bank", "url": "https://www.worldbank.org"}
    jba = {"name": "JBA Risk Management", "url": "https://www.jbarisk.com"}
    return {
        "id": "freetown-flood-risk-assessment",
        "title": "Freetown, Sierra Leone Flood Risk Assessment",
        "description": "River and surface-water flood depths on 1 m grids at return periods of 20 to 1500 years, 153,436 modelled buildings, and the event loss tables and exceedance curves of a catalogue of 10,000 simulated years.",
        "risk_data_type": ["hazard", "exposure", "loss"],
        "publisher": bank,
        "contact_point": bank,
        "creator": jba,
        "spatial": {"scale": "urban", "countries": ["SLE"], "bbox": [-13.3032, 8.3417, -13.1151, 8.5038]},
        "license": "CC-BY-4.0",
        "resources": [
            {"id": "hazard-maps", "title": "Flood hazard maps", "description": "Depth at 1 m, by peril and return period, in WGS 84.", "access_url": base, "format": "tif", "coordinate_system": "EPSG:4326"},
            {"id": "exposure", "title": "Modelled exposure", "description": "Buildings at their places, with their values and occupancies.", "access_url": base, "format": "csv"},
            {"id": "results", "title": "Probabilistic analysis results", "description": "Event loss tables, annual averages and exceedance curves, by occupancy and for the whole city.", "access_url": base, "format": "csv"},
            {"id": "loss-maps", "title": "Flood hazard and loss maps", "description": "Annual losses by ward, and the places the City Council records as flooding.", "access_url": base, "format": "shp"},
        ],
        "attributions": [
            {"id": "city", "entity": {"name": "Freetown City Council", "url": "https://fcc.gov.sl"}, "role": "contributor"},
        ],
    }


def main() -> int:
    schema = json.loads((CACHE / "rdls-schema" / "rdls_schema.json").read_text(encoding="utf-8"))
    record = freetown()
    wrong = validate(record, schema)
    if wrong:
        print("\n".join(wrong), file=sys.stderr)
        return 1
    out = CACHE / "products" / "rdls-freetown.json"
    out.write_text(json.dumps({"datasets": [record]}, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"a record of RDLS, to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
