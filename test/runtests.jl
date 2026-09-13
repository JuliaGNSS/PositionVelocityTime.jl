
using Test, PositionVelocityTime, GNSSDecoder, AstroTime, GNSSSignals, Geodesy, Dates, LinearAlgebra
using Unitful: Hz, m, s, °, ustrip
# `Dictionaries.Dictionary` is how GNSSDecoder keys its paged records (BeiDou B1C's
# per-GNSS BGTO sets, the almanacs); the BeiDou tests build one directly, and it also
# backs the satellites of a `PositionVelocityTime.SignalGroup`.
using Dictionaries: Dictionaries, Dictionary
# The fixtures build satellites as flat vectors, so most of the suite routes them
# through `signal_groups` — the documented bridge — rather than naming its groups.
# `SignalGroup` itself stays fully qualified everywhere: `Tracking` exports a type of
# the same name, and `test/tracking_ext.jl` does `using Tracking`, so the deliberate
# name collision is live in this very session.
using PositionVelocityTime: signal_groups

# Two shorthands the suite uses wherever it exercises the measurement-model surface
# directly rather than through `calc_pvt`. They live here rather than in `fixtures.jl`,
# which stays strictly the data shared with the benchmark script.
#
# `measurement_rows` is the flat `SatelliteMeasurement` rows behind a vector of
# `SatelliteState`s — what the collection pass hands the solver. `approximate_year` is
# pinned for the same reason every `calc_pvt` call here pins it: the fixtures were
# recorded in 2021 and GPS L1 C/A's week number is 10 bits.
measurement_rows(states; approximate_year = 2021) = first(
    PositionVelocityTime.collect_measurements(signal_groups(states); approximate_year),
)
# `carrier_hz` is a ranging signal's carrier as the plain `Float64` in Hz that
# `ionospheric_delay` takes — the `center_frequency` field of a measurement row.
carrier_hz(system) = ustrip(Hz, get_center_frequency(system))

include("aqua.jl")
include("fixtures.jl")
include("sat_time.jl")
include("sat_position.jl")
include("signal_groups.jl")
include("pvt.jl")
include("dop.jl")
include("cnav.jl")
include("gps_l2c.jl")
include("galileo_e5a.jl")
include("galileo_e5b_e6b.jl")
include("beidou.jl")
include("inter_frequency_bias.jl")
include("pvt_iono_tropo.jl")
include("get_week.jl")
include("tracking_ext.jl")
include("ionosphere.jl")
include("troposphere.jl")
include("pvt_integration.jl")
