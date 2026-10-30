
using Test,
    PositionVelocityTime, GNSSDecoder, AstroTime, GNSSSignals, Geodesy, Dates, LinearAlgebra
using Unitful: Hz, m, s, °, ustrip
# Dictionaries preserve the decoder PRNs and the concrete observation type.
using Dictionaries: Dictionaries, Dictionary
using PositionVelocityTime: SVector

function signal_group(signal, states)
    all(state -> get_signal_id(state.signal) === get_signal_id(signal), states) ||
        error("put each ranging signal in its own observation collection")
    Dictionary([state.decoder.prn for state in states], states)
end
signal_group(states) = signal_group(first(states).signal, states)

# Successful fixtures expose their solution after checking the explicit outcome.
# Failure tests call the public functions directly and inspect both result fields.
function solved_pvt(args...; kwargs...)
    result = calc_pvt(args...; kwargs...)
    @test result.status === PositionVelocityTime.solved
    result.solution
end
function solved_pvt!(args...; kwargs...)
    result = calc_pvt!(args...; kwargs...)
    @test result.status === PositionVelocityTime.solved
    result.solution
end

# Two more shorthands, for the tests that exercise the measurement-model surface
# directly rather than through `calc_pvt`. They live here rather than in `fixtures.jl`,
# which stays strictly the data shared with the benchmark script.
#
# `measurements_of` gives the `PreparedObservation`s behind one epoch's groups —
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
    (
        get(fields, name, getfield(ephemeris, name)) for
        name in fieldnames(BroadcastEphemeris)
    )...,
)
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

# Relabels a SignalObservation's decoder PRN. Nothing in the least-squares solve depends on
# the PRN — it identifies a satellite, as the `observations`-dict key and for the layout's
# distinct-satellite count — so this isolates that identity from the measurements.
function with_prn(state, prn)
    new_decoder = GNSSDecoder.GNSSDecoderState(state.decoder; prn)
    SignalObservation(;
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
    for _ = 1:10
        t_t = t - ρ / c
        ρ =
            PositionVelocityTime.predict_pseudoranges!(
                zeros(1),
                [sat_at(t_t)],
                [user_at(t); bias_at(t)],
                ONE_CLOCK_COLUMN,
            )[1] - sat_clock_at(t_t)
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
