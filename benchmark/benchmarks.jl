using BenchmarkTools
using PositionVelocityTime

include("fixtures.jl")

const SUITE = BenchmarkGroup()

const GALILEO_INPUT = make_galileo_input()
const GPS_INPUT = make_gps_input()

# Sanity-warm the inputs so we can supply prev_pvt for warm-start benchmarks.
const GALILEO_PREV = calc_pvt(GALILEO_INPUT)
const GPS_PREV = calc_pvt(GPS_INPUT)

SUITE["calc_pvt"] = BenchmarkGroup()

for (system_label, input, prev_pvt) in (
    ("GalileoE1B", GALILEO_INPUT, GALILEO_PREV),
    ("GPSL1", GPS_INPUT, GPS_PREV),
)
    SUITE["calc_pvt"][system_label] = BenchmarkGroup()
    n_total = num_satellites(input)
    sat_counts = n_total == 4 ? (4,) : (4, n_total)
    for n in sat_counts
        subset = take_satellites(input, n)
        group = BenchmarkGroup()
        group["cold"] = @benchmarkable calc_pvt($subset)
        group["warm"] = @benchmarkable calc_pvt($subset, $prev_pvt)
        SUITE["calc_pvt"][system_label]["$(n)sats"] = group
    end
end

# The mixed-constellation case, which is what signal groups exist for: pooled into one
# vector its element type is abstract, and every per-satellite call inside the solve is
# a dynamic dispatch. This is the comparison that carries the whole change, so it is
# built through `combine_inputs` and runs on both revisions.
const MIXED_INPUT = combine_inputs((gps = GPS_INPUT, galileo = GALILEO_INPUT))
const MIXED_PREV = calc_pvt(MIXED_INPUT)
const MIXED_SATS = num_satellites(GPS_INPUT) + num_satellites(GALILEO_INPUT)

SUITE["calc_pvt"]["GPSL1+GalileoE1B"] = BenchmarkGroup()
mixed_group = BenchmarkGroup()
mixed_group["cold"] = @benchmarkable calc_pvt($MIXED_INPUT)
mixed_group["warm"] = @benchmarkable calc_pvt($MIXED_INPUT, $MIXED_PREV)
SUITE["calc_pvt"]["GPSL1+GalileoE1B"]["$(MIXED_SATS)sats"] = mixed_group

# Two 6.0-only measurements, absent from the base revision's run (AirspeedVelocity
# reports a benchmark only present on one side as new rather than as a regression):
# the documented flat-vector bridge, whose cost has to stay bounded to the collection
# pass, and that collection pass on its own — the half that specialises on the group
# shape, and the half a consumer running its own estimator over the measurement model
# pays.
if HAS_SIGNAL_GROUPS
    MIXED_STATES = make_mixed_states()
    mixed_group["cold, via signal_groups"] =
        @benchmarkable calc_pvt(PositionVelocityTime.signal_groups($MIXED_STATES))

    SUITE["collect_measurements"] = BenchmarkGroup()
    SUITE["collect_measurements"]["GPSL1"] =
        @benchmarkable PositionVelocityTime.collect_measurements($GPS_INPUT)
    SUITE["collect_measurements"]["GPSL1+GalileoE1B"] =
        @benchmarkable PositionVelocityTime.collect_measurements($MIXED_INPUT)
end
