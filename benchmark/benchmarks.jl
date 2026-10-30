using BenchmarkTools
using PositionVelocityTime

include("fixtures.jl")

const SUITE = BenchmarkGroup()

const GALILEO_INPUT = make_galileo_input()
const GPS_INPUT = make_gps_input()

# The observations were recorded in 2021; validate seed fixes and pin GPS rollover.
function benchmark_seed(input)
    result = calc_pvt(input; approximate_year = 2021)
    result.status === PositionVelocityTime.solved || error("benchmark input did not solve")
    result.solution
end
const GALILEO_PREV = benchmark_seed(GALILEO_INPUT)
const GPS_PREV = benchmark_seed(GPS_INPUT)

SUITE["calc_pvt"] = BenchmarkGroup()

for (system_label, input, prev_pvt) in
    (("GalileoE1B", GALILEO_INPUT, GALILEO_PREV), ("GPSL1", GPS_INPUT, GPS_PREV))
    SUITE["calc_pvt"][system_label] = BenchmarkGroup()
    n_total = length(input)
    sat_counts = n_total == 4 ? (4,) : (4, n_total)
    for n in sat_counts
        subset = take_satellites(input, n)
        group = BenchmarkGroup()
        group["cold"] = @benchmarkable calc_pvt($subset; approximate_year = 2021)
        group["warm"] = @benchmarkable calc_pvt($subset, $prev_pvt; approximate_year = 2021)
        SUITE["calc_pvt"][system_label]["$(n)sats"] = group
    end
end

# The mixed-constellation case, which is why observation types are kept in separate collections: pooled into one
# vector its element type is abstract, and every per-satellite call inside the solve is
# a dynamic dispatch.
const MIXED_INPUT = (gps = GPS_INPUT, galileo = GALILEO_INPUT)
const MIXED_PREV = benchmark_seed(MIXED_INPUT)
const MIXED_SATS = length(GPS_INPUT) + length(GALILEO_INPUT)

SUITE["calc_pvt"]["GPSL1+GalileoE1B"] = BenchmarkGroup()
mixed_group = BenchmarkGroup()
mixed_group["cold"] = @benchmarkable calc_pvt($MIXED_INPUT; approximate_year = 2021)
mixed_group["warm"] =
    @benchmarkable calc_pvt($MIXED_INPUT, $MIXED_PREV; approximate_year = 2021)
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
    group["cold"] = @benchmarkable calc_pvt!(
        $solution,
        $workspace,
        $input,
        $cold;
        approximate_year = 2021,
    )
    group["warm"] = @benchmarkable calc_pvt!(
        $solution,
        $workspace,
        $input,
        $prev_pvt;
        approximate_year = 2021,
    )
    SUITE["calc_pvt!"][label] = group
end

# The collection pass on its own: the half that specialises on the group shape, and the
# half a consumer running its own estimator over the measurement model pays.
SUITE["prepare_observations"] = BenchmarkGroup()
SUITE["prepare_observations"]["GPSL1"] =
    @benchmarkable PositionVelocityTime.prepare_observations(
        $GPS_INPUT;
        approximate_year = 2021,
    )
SUITE["prepare_observations"]["GPSL1+GalileoE1B"] =
    @benchmarkable PositionVelocityTime.prepare_observations(
        $MIXED_INPUT;
        approximate_year = 2021,
    )

# The same mixed solve after preparation, isolating the shared numerical model.
const MIXED_PREPARATION =
    PositionVelocityTime.prepare_observations(MIXED_INPUT; approximate_year = 2021)
const MIXED_PREPARED = first(MIXED_PREPARATION)
const MIXED_PARAMETERS = last(MIXED_PREPARATION)
SUITE["prepared calc_pvt"] =
    @benchmarkable calc_pvt($MIXED_PREPARED; ionospheric_model = $MIXED_PARAMETERS)
