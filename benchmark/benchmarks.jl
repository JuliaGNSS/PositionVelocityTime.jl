using BenchmarkTools
using PositionVelocityTime
using PositionVelocityTime: SignalGroup

include("fixtures.jl")

const SUITE = BenchmarkGroup()

const GALILEO_GROUP = make_galileo_group()
const GPS_GROUP = make_gps_group()

# Sanity-warm the inputs so we can supply prev_pvt for warm-start benchmarks.
const GALILEO_PREV = calc_pvt(GALILEO_GROUP)
const GPS_PREV = calc_pvt(GPS_GROUP)

SUITE["calc_pvt"] = BenchmarkGroup()

for (system_label, group, prev_pvt) in (
    ("GalileoE1B", GALILEO_GROUP, GALILEO_PREV),
    ("GPSL1", GPS_GROUP, GPS_PREV),
)
    SUITE["calc_pvt"][system_label] = BenchmarkGroup()
    n_total = length(group.satellites)
    sat_counts = n_total == 4 ? (4,) : (4, n_total)
    for n in sat_counts
        subset = take_satellites(group, n)
        bgroup = BenchmarkGroup()
        bgroup["cold"] = @benchmarkable calc_pvt($subset)
        bgroup["warm"] = @benchmarkable calc_pvt($subset, $prev_pvt)
        SUITE["calc_pvt"][system_label]["$(n)sats"] = bgroup
    end
end

# The mixed-constellation case, which is what signal groups exist for: pooled into one
# vector its element type was abstract, and every per-satellite call inside the solve
# was a dynamic dispatch. Benchmarked both ways a receiver can hand it over — named
# groups (the fast path, types known at compile time) and `signal_groups` over a flat
# vector (the documented bridge, inference-blind by construction) — so the cost of the
# bridge stays visible and bounded to the collection pass.
const MIXED_GROUPS = (gps = GPS_GROUP, galileo = GALILEO_GROUP)
const MIXED_STATES = [collect(GPS_GROUP.satellites); collect(GALILEO_GROUP.satellites)]
const MIXED_PREV = calc_pvt(MIXED_GROUPS)

SUITE["calc_pvt"]["GPSL1+GalileoE1B"] = BenchmarkGroup()
mixed_group = BenchmarkGroup()
mixed_group["cold"] = @benchmarkable calc_pvt($MIXED_GROUPS)
mixed_group["warm"] = @benchmarkable calc_pvt($MIXED_GROUPS, $MIXED_PREV)
mixed_group["cold, via signal_groups"] =
    @benchmarkable calc_pvt(PositionVelocityTime.signal_groups($MIXED_STATES))
SUITE["calc_pvt"]["GPSL1+GalileoE1B"]["$(length(MIXED_STATES))sats"] = mixed_group

# The function barrier itself, separately from the solver behind it: this is the half
# that specialises on the group shape, and the half a consumer running its own
# estimator over the measurement model pays.
SUITE["collect_measurements"] = BenchmarkGroup()
SUITE["collect_measurements"]["GPSL1"] =
    @benchmarkable PositionVelocityTime.collect_measurements($GPS_GROUP)
SUITE["collect_measurements"]["GPSL1+GalileoE1B"] =
    @benchmarkable PositionVelocityTime.collect_measurements($MIXED_GROUPS)
