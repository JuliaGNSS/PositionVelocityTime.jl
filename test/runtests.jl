
using Test, PositionVelocityTime, GNSSDecoder, AstroTime, GNSSSignals, Geodesy, Dates, LinearAlgebra
using Unitful: Hz, m, s, °, ustrip
# `Dictionaries.Dictionary` is how GNSSDecoder keys its paged records (BeiDou B1C's
# per-GNSS BGTO sets, the almanacs); the BeiDou tests build one directly, and it also
# backs the satellites of a `PositionVelocityTime.SignalGroup`.
using Dictionaries: Dictionaries, Dictionary
using PositionVelocityTime: SVector
# The fixtures build one flat vector of `SatelliteState`s per ranging signal. This is
# how the suite turns one of them into the `SignalGroup` that `calc_pvt` takes; an epoch
# spanning several signals is written as several groups — a NamedTuple where the names
# carry meaning, a bare tuple (numbered `group1`, `group2`, …) where they do not.
#
# `SignalGroup` stays fully qualified here and everywhere below, as it must wherever
# `Tracking`, which exports a type of the same name, is loaded too (as in
# `test/integration`).
function signal_group(signal, states)
    ids = unique(get_signal_id(state.signal) for state in states)
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
signal_group(states) = signal_group(first(states).signal, states)

# Two more shorthands, for the tests that exercise the measurement-model surface
# directly rather than through `calc_pvt`. They live here rather than in `fixtures.jl`,
# which stays strictly the data shared with the benchmark script.
#
# `measurements_of` gives the `SatelliteMeasurement`s behind one epoch's groups —
# what the collection pass hands the solver. `approximate_year` is pinned for the same
# reason every `calc_pvt` call here pins it: the fixtures were recorded in 2021 and
# GPS L1 C/A's week number is 10 bits.
measurements_of(groups; approximate_year = 2021) =
    first(PositionVelocityTime.collect_measurements(groups; approximate_year))
# `has_offset` is whether GNSSDecoder reports a broadcast time offset toward `target`
# (at any time of week: whether there is one does not depend on it); `steering` evaluates
# it at the time of week `t` of the decoder's week.
has_offset(decoder, target) =
    !isnothing(GNSSDecoder.get_broadcast_time_offset(decoder, target, 0.0))
steering(decoder, target, t) =
    something(GNSSDecoder.get_broadcast_time_offset(decoder, target, t))
# `ephemeris` is a decoder's broadcast ephemeris and clock, the record the propagator and
# the clock model read; `with_fields` is a copy of one with some of its fields replaced.
ephemeris(decoder) = something(GNSSDecoder.get_broadcast_ephemeris(decoder))
with_fields(ephemeris::BroadcastEphemeris; fields...) = BroadcastEphemeris(
    (get(fields, name, getfield(ephemeris, name)) for name in fieldnames(BroadcastEphemeris))...)
# `gps_epoch` and `bdt_epoch` date a time of week in week 0 of GPS and BeiDou Time, the
# epoch the ionospheric models take; the Klobuchar models read only its time of week.
gps_epoch(time_of_week) = get_week_start_epoch(GPST(), 0) + time_of_week
bdt_epoch(time_of_week) = get_week_start_epoch(BDT(), 0) + time_of_week
# `carrier_hz` is a ranging signal's carrier as the plain `Float64` in Hz that
# `ionospheric_delay` takes — the `center_frequency` field of a measurement.
carrier_hz(system) = ustrip(Hz, get_center_frequency(system))
# The inferred return type of `f(::types...)`. `Base.infer_return_type` says this in one
# call but only from Julia 1.11; `Base.return_types` is the spelling that also works on
# the 1.10 this package supports, and `only` is the assertion that the signature picks
# out exactly one method — which is what makes the two equivalent here.
inferred_return_type(f, types) = only(Base.return_types(f, types))

# Relabels a SatelliteState's decoder PRN. Nothing in the least-squares solve depends on
# the PRN — it identifies a satellite, as the `satellites`-dict key and for the layout's
# distinct-satellite count — so this isolates that identity from the measurements.
function with_prn(state, prn)
    new_decoder = GNSSDecoder.GNSSDecoderState(state.decoder; prn)
    SatelliteState(;
        decoder = new_decoder,
        signal = state.signal,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end


const ONE_CLOCK_COLUMN = PositionVelocityTime.BiasColumns([1], 1, [0], 0)

# The pseudorange at receiver time `t`, as `predict_pseudoranges!` models it: the satellite
# taken at its transmit time `t_t = t − ρ/c` (found by fixed-point iteration), the receiver
# at `user_at(t)` with clock bias `bias_at(t)` (m), and the satellite clock offset
# `sat_clock_at(t_t)` (m) subtracted.
function modelled_pseudorange(sat_at, sat_clock_at, user_at, bias_at, t)
    c = PositionVelocityTime.SPEED_OF_LIGHT
    ρ = 0.0
    t_t = t
    for _ in 1:10
        t_t = t - ρ / c
        ρ = PositionVelocityTime.predict_pseudoranges!(
            zeros(1), [sat_at(t_t)], [user_at(t); bias_at(t)], ONE_CLOCK_COLUMN)[1] -
            sat_clock_at(t_t)
    end
    return ρ, t_t
end

central_difference(f, t, h) = (f(t + h) - f(t - h)) / 2h

include("aqua.jl")
include("api.jl")
include("fixtures.jl")
# One file per source file, in the source's order.
include("collection.jl")
include("solution.jl")
include("bias_layout.jl")
include("measurement_model.jl")
include("estimation.jl")
include("satellite_time.jl")
include("satellite_position.jl")
include("atmosphere/klobuchar.jl")
include("atmosphere/ntcm_g.jl")
include("atmosphere/bdgim.jl")
include("atmosphere/troposphere.jl")
include("atmosphere/delays.jl")
include("solver.jl")
# Across the source files: each supported signal end to end, and the allocation-free
# contract. A real capture through Tracking.jl is `test/integration`'s, in its own
# environment.
include("signals/cnav.jl")
include("signals/gps_l2c.jl")
include("signals/galileo_e5a.jl")
include("signals/galileo_e5b_e6b.jl")
include("signals/beidou.jl")
include("allocation_free.jl")
