
using Test, PositionVelocityTime, GNSSDecoder, AstroTime, GNSSSignals, Geodesy, Dates, LinearAlgebra
using Unitful: Hz, m, s, °, ustrip
# `Dictionaries.Dictionary` is how GNSSDecoder keys its paged records (BeiDou B1C's
# per-GNSS BGTO sets, the almanacs); the BeiDou tests build one directly, and it also
# backs the satellites of a `PositionVelocityTime.SignalGroup`.
using Dictionaries: Dictionaries, Dictionary
# The fixtures build one flat vector of `SatelliteState`s per ranging signal. This is
# how the suite turns one of them into the `SignalGroup` that `calc_pvt` takes; an epoch
# spanning several signals is written as several groups — a NamedTuple where the names
# carry meaning, a bare tuple (numbered `group1`, `group2`, …) where they do not.
#
# `SignalGroup` stays fully qualified here and everywhere below: `Tracking` exports a
# type of the same name, and `test/tracking_ext.jl` does `using Tracking`, so the
# deliberate collision is live in this very session.
function signal_group(signal, states)
    ids = unique(get_signal_id(state.system) for state in states)
    all(id -> id === get_signal_id(signal), ids) || error(
        "signal_group takes the satellites of one ranging signal, got $ids — " *
        "put each signal in its own group",
    )
    PositionVelocityTime.SignalGroup(
        signal,
        Dictionary([state.decoder.prn for state in states], states),
    )
end
# The signal read off the satellites, for the usual case of a non-empty group. A group
# that may be empty has to be given its signal explicitly — which is exactly why a
# `SignalGroup` carries one rather than deriving it from its satellites.
signal_group(states) = signal_group(first(states).system, states)

# Two more shorthands, for the tests that exercise the measurement-model surface
# directly rather than through `calc_pvt`. They live here rather than in `fixtures.jl`,
# which stays strictly the data shared with the benchmark script.
#
# `measurement_rows` is the flat `SatelliteMeasurement` rows behind one epoch's groups —
# what the collection pass hands the solver. `approximate_year` is pinned for the same
# reason every `calc_pvt` call here pins it: the fixtures were recorded in 2021 and
# GPS L1 C/A's week number is 10 bits.
measurement_rows(groups; approximate_year = 2021) =
    first(PositionVelocityTime.collect_measurements(groups; approximate_year))
# `carrier_hz` is a ranging signal's carrier as the plain `Float64` in Hz that
# `ionospheric_delay` takes — the `center_frequency` field of a measurement row.
carrier_hz(system) = ustrip(Hz, get_center_frequency(system))
# The inferred return type of `f(::types...)`. `Base.infer_return_type` says this in one
# call but only from Julia 1.11; `Base.return_types` is the spelling that also works on
# the 1.10 this package supports, and `only` is the assertion that the signature picks
# out exactly one method — which is what makes the two equivalent here.
inferred_return_type(f, types) = only(Base.return_types(f, types))

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
