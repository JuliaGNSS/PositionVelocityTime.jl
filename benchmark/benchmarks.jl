using BenchmarkTools
using PositionVelocityTime

include("fixtures.jl")

const SUITE = BenchmarkGroup()

const GALILEO_INPUT = make_galileo_input()
const GPS_INPUT = make_gps_input()

# Sanity-warm the inputs so we can supply prev_pvt for warm-start benchmarks.
const GALILEO_PREV = something(calc_pvt(GALILEO_INPUT).solution)
const GPS_PREV = something(calc_pvt(GPS_INPUT).solution)

SUITE["calc_pvt"] = BenchmarkGroup()

for (system_label, input, prev_pvt) in
    (("GalileoE1B", GALILEO_INPUT, GALILEO_PREV), ("GPSL1", GPS_INPUT, GPS_PREV))
    SUITE["calc_pvt"][system_label] = BenchmarkGroup()
    n_total = length(input)
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
# a dynamic dispatch.
const MIXED_INPUT = (gps = GPS_INPUT, galileo = GALILEO_INPUT)
const MIXED_PREV = something(calc_pvt(MIXED_INPUT).solution)
const MIXED_SATS = length(GPS_INPUT) + length(GALILEO_INPUT)

SUITE["calc_pvt"]["GPSL1+GalileoE1B"] = BenchmarkGroup()
mixed_group = BenchmarkGroup()
mixed_group["cold"] = @benchmarkable calc_pvt($MIXED_INPUT)
mixed_group["warm"] = @benchmarkable calc_pvt($MIXED_INPUT, $MIXED_PREV)
SUITE["calc_pvt"]["GPSL1+GalileoE1B"]["$(MIXED_SATS)sats"] = mixed_group

# The in-place solve: the same epochs, overwriting one solution with a reused workspace,
# which allocates nothing.
SUITE["calc_pvt!"] = BenchmarkGroup()
for (label, input, prev_pvt) in
    (("GPSL1", GPS_INPUT, GPS_PREV), ("GPSL1+GalileoE1B", MIXED_INPUT, MIXED_PREV))
    group = BenchmarkGroup()
    workspace = PVTWorkspace()
    solution = PVTSolution()
    cold = PVTSolution()
    group["cold"] = @benchmarkable calc_pvt!($solution, $workspace, $input, $cold)
    group["warm"] = @benchmarkable calc_pvt!($solution, $workspace, $input, $prev_pvt)
    SUITE["calc_pvt!"][label] = group
end

# The collection pass on its own: the half that specialises on the group shape, and the
# half a consumer running its own estimator over the measurement model pays.
SUITE["collect_measurements"] = BenchmarkGroup()
SUITE["collect_measurements"]["GPSL1"] =
    @benchmarkable PositionVelocityTime.collect_measurements($GPS_INPUT)
SUITE["collect_measurements"]["GPSL1+GalileoE1B"] =
    @benchmarkable PositionVelocityTime.collect_measurements($MIXED_INPUT)
