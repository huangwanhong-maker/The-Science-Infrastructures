// A new script in each language: the smallest thing that runs, for a record of events.

export const TEMPLATES = {
  py: `"""What this does, its source, and its limits."""

import numpy as np


def create(api):
    class Mine(api.SpaceTimeAnalysis):
        id = "mine"
        title = "The largest events"
        description = "The events shown, keeping those above a value."
        params = [api.param("share", "Share kept", 0.1, minimum=0.001, maximum=1)]

        def run(self, records, params, context):
            r = next((q for q in records if q.kind == "events"), None)
            if r is None:
                raise ValueError("This needs events: show a record of events first.")
            name = next(iter(r.values), None)
            v = r.values[name] if name else r.time
            keep = v >= np.nanquantile(v, 1 - params["share"])
            return [
                api.events(f"The largest of {r.title}", r.time[keep], r.lat[keep], r.lon[keep],
                           values={k: c[keep] for k, c in r.values.items()}, units=r.units),
                api.note(f"{int(keep.sum()):,} of {r.time.size:,} events."),
            ]

    return Mine()
`,
  js: `// What this does, its source, and its limits.

export function create(api) {
  return new class extends api.SpaceTimeAnalysis {
    id = 'mine';
    title = 'The largest events';
    description = 'The events shown, keeping those above a value.';
    params = [api.param('share', 'Share kept', 0.1, { min: 0.001, max: 1 })];

    run(records, params, context) {
      const r = records.find((q) => q.kind === 'events');
      if (!r) throw new Error('This needs events: show a record of events first.');
      const name = Object.keys(r.values)[0], v = name ? r.values[name] : r.time;
      const line = [...v].filter(Number.isFinite).sort((a, b) => a - b)[Math.floor(v.length * (1 - params.share))];
      const keep = v.flatMap((x, k) => (x >= line ? [k] : []));
      const of = (c) => keep.map((k) => c[k]);
      return [
        api.events('The largest of ' + r.title, { time: of(r.time), lat: of(r.lat), lon: of(r.lon), values: Object.fromEntries(Object.entries(r.values).map(([k, c]) => [k, of(c)])), units: r.units }),
        api.note(keep.length + ' of ' + r.time.length + ' events.'),
      ];
    }
  }();
}
`,
  cpp: `// What this does, its source, and its limits.

#include "spacetime_plugin.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace st = science::spacetime;

struct Mine : st::Analysis {
  science::Json describe() const override {
    return science::Json::object({{"id", "mine"}, {"title", "Every n-th event"}, {"description", "A thinned copy of the events shown."},
                                  {"params", science::Json::array({st::param("n", "Keep one in", 10, "integer", 1, 1000)})}});
  }
  std::vector<science::Json> run(const std::vector<st::Record>& records, const st::Params& p, const st::Context&, st::Progress&) override {
    const auto found = std::find_if(records.begin(), records.end(), [](const st::Record& r) { return r.kind == "events"; });
    if (found == records.end()) throw std::invalid_argument("This needs events: show a record of events first.");
    st::Events out;
    for (std::size_t k = 0; k < found->events.size(); k += static_cast<std::size_t>(p.integer("n", 10)))
      out.push(found->events.time[k], found->events.lat[k], found->events.lon[k]);
    return {out.to_json("Thinned: " + found->title), st::note(std::to_string(out.size()) + " events kept.")};
  }
};

SCIENCE_SPACETIME_PLUGIN(Mine)
`,
};
