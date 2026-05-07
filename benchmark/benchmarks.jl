using BenchmarkTools
using PositionVelocityTime

include("fixtures.jl")

const SUITE = BenchmarkGroup()

const GALILEO_STATES = make_galileo_states()
const GPS_STATES = make_gps_states()

# Sanity-warm the inputs so we can supply prev_pvt for warm-start benchmarks.
const GALILEO_PREV = calc_pvt(GALILEO_STATES)
const GPS_PREV = calc_pvt(GPS_STATES)

SUITE["calc_pvt"] = BenchmarkGroup()

for (system_label, all_states, prev_pvt) in (
    ("GalileoE1B", GALILEO_STATES, GALILEO_PREV),
    ("GPSL1", GPS_STATES, GPS_PREV),
)
    SUITE["calc_pvt"][system_label] = BenchmarkGroup()
    n_total = length(all_states)
    sat_counts = n_total == 4 ? (4,) : (4, n_total)
    for n in sat_counts
        states = all_states[1:n]
        group = BenchmarkGroup()
        group["cold"] = @benchmarkable calc_pvt($states)
        group["warm"] = @benchmarkable calc_pvt($states, $prev_pvt)
        SUITE["calc_pvt"][system_label]["$(n)sats"] = group
    end
end
