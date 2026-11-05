"""Seed regime-cases.jsonl from the cases a run printed before it was lost (WSL restarted).

The log keeps each case's value, person-hours, person-hours at level F, run time and the regimes chosen; the per-epoch
rewards were not printed, so the seeded records carry the totals and mark themselves "from_log".
"""
import json
import re
import sys

log, journal = sys.argv[1], sys.argv[2]
pat = re.compile(r"^\s+scale ([0-9.]+), day (\d+), (\S+)\s+value (-?[0-9.]+) \(([0-9.]+) person-hours, ([0-9.]+) at level F\) in ([0-9.]+) s; regimes((?: \d)+)(?: \((\d+) restarts\))?\s*$")
seen = set()
try:
    for line in open(journal, encoding="utf-8"):
        if line.strip():
            j = json.loads(line)
            seen.add((j["scale"], j["day"], j["policy"]))
except FileNotFoundError:
    pass
added = 0
with open(journal, "a", encoding="utf-8") as out:
    for line in open(log, encoding="utf-8"):
        m = pat.match(line)
        if not m:
            continue
        scale, day, policy = float(m[1]), float(m[2]), m[3]
        if (scale, day, policy) in seen:
            continue
        rec = {"scale": scale, "day": day, "policy": policy,
               "actions": [float(a) for a in m[8].split()],
               "value": float(m[4]), "restarts": int(m[9] or 0), "seconds": float(m[7]),
               "totals": {"person_hours": float(m[5]), "los_f_hours": float(m[6])},
               "from_log": True}
        out.write(json.dumps(rec) + "\n")
        seen.add((scale, day, policy))
        added += 1
print(f"seeded {added} cases; {len(seen)} in the journal")
