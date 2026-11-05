// urban-science's movement: walkers as agents on the walking network (theory/movement.md; ADR 0071, ADR 0075).
//   network.hpp      the walking network, its directed segments, widths, grades, crossings and signals
//   route.hpp        the logit over walks: routes, targets, the walk sums z and π, the expected flow
//   walkers.hpp      the agent-based model: streams of trips, the step, admission, signals, tallies, positions
//   calibration.hpp  departures fitted to the counters each way, hour by hour (NNLS)
//   regime.hpp       a street's regime as a POMDP on the walkers, solved online by POMCP; the exactly solvable corridor
//   record.hpp       a run written as a record of simulated events
#pragma once

#include <urban/movement/calibration.hpp>
#include <urban/movement/network.hpp>
#include <urban/movement/record.hpp>
#include <urban/movement/regime.hpp>
#include <urban/movement/route.hpp>
#include <urban/movement/walkers.hpp>
