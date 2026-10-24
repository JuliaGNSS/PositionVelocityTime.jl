# The collection pass: satellite *states* in, flat satellite *measurements* out.
#
# Everything that depends on a satellite's navigation-data type — decoder, ranging signal,
# propagator, clock polynomial, broadcast time offsets — is evaluated here, once per
# satellite, into the one concrete `SatelliteMeasurement` row the solver works on. This
# is `calc_pvt`'s function barrier: upstream specialises on the group shape (and is
# small), downstream compiles once for every constellation mix.

"""
    SatelliteState{CP<:Real,D<:GNSSDecoderState,S<:AbstractGNSSSignal}

One satellite's decoder state together with its code and carrier measurements.

# Fields
- `decoder::GNSSDecoderState`: the decoded navigation data.
- `signal::AbstractGNSSSignal`: the ranging signal the code phase and Doppler were
  measured on (e.g. `GPSL1CA()`, `GalileoE1C()`); may be a dataless pilot whose data
  `decoder` decoded from the band's data component.
- `code_phase::CP`: code phase in chips.
- `carrier_doppler`: carrier Doppler in Hz.
- `carrier_phase::CP`: carrier phase in radians, as `Tracking.get_carrier_phase`
  (default `0.0`).

# Constructors
    SatelliteState(; decoder, signal, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, signal, tracked_satellite)

The second reads the measurements off a `Tracking.TrackedSat`; Tracking.jl (9.1 on)
defines it, so it is available once `Tracking` is loaded.
"""
@kwdef struct SatelliteState{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    signal::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
end

"""
    PositionVelocityTime.SignalGroup(signal, satellites)

One epoch's measurements of the satellites tracked on one ranging `signal`: the
counterpart of `Tracking.SignalGroup`, which holds their loop state, with the same group
and PRN keys. Unexported, so `using Tracking, PositionVelocityTime` does not clash.

Within a group every satellite has one concrete `SatelliteState` type, so the collection
pass dispatches statically on it (see [`collect_measurements`](@ref)).

# Fields
- `signal::AbstractGNSSSignal`: the ranging signal of every satellite in the group, known
  even when the group is empty.
- `satellites`: a `Dictionary{Int,<:SatelliteState}` keyed by PRN, or an
  `AbstractVector{<:SatelliteState}`. Iteration order is significant (see
  [`calc_pvt`](@ref)).

Every satellite's `state.signal` must be `signal`; this is not checked.
"""
struct SignalGroup{S<:AbstractGNSSSignal,C}
    signal::S
    satellites::C
end

"""
    PositionVelocityTime.SignalGroups{N}

The shape [`calc_pvt`](@ref) takes: a `NamedTuple` of `N` [`SignalGroup`](@ref)s, as
`Tracking.SignalGroups`, so one set of group keys runs through tracking and the fix.
"""
const SignalGroups{N} = NamedTuple{<:Any,<:NTuple{N,SignalGroup}}

# As `Tracking._normalize_signal_groups`: a NamedTuple passes through, a single group
# becomes `:default`, and a bare tuple is numbered `group1`, `group2`, … (not named after
# the signals, which two groups may share).
_normalize_signal_groups(groups::NamedTuple{<:Any,<:Tuple{Vararg{SignalGroup}}}) = groups
_normalize_signal_groups(group::SignalGroup) = (default = group,)
@generated function _normalize_signal_groups(groups::Tuple{Vararg{SignalGroup}})
    names = ntuple(i -> Symbol(:group, i), length(groups.parameters))
    :(NamedTuple{$names}(groups))
end

"""
    PositionVelocityTime.signal_groups(track_state, decoders) -> NamedTuple of SignalGroup

Build one epoch's [`SignalGroups`](@ref) from a receiver's own per-signal state. Tracking.jl
(9.1 on) defines the method, taking a `Tracking.TrackState` and the receiver's decoders.

There is no method taking a flat `Vector{SatelliteState}`: grouping a pooled vector each
epoch would bring back the dynamic dispatch groups exist to remove. Build the groups where
the satellites are tracked.
"""
function signal_groups end

# Anything else is refused at the entry point, with directions, rather than failing deep
# inside the collection pass.
_normalize_signal_groups(groups) = throw_not_signal_groups(
    "every group must be a `PositionVelocityTime.SignalGroup`, got a `$(typeof(groups))`",
)

function throw_not_signal_groups(what)
    throw(
        ArgumentError(
            string(
                "`calc_pvt` takes signal groups, ",
                what,
                ". Group the satellites by their ranging signal:\n\n",
                "    using PositionVelocityTime: SignalGroup\n",
                "    calc_pvt((gps     = SignalGroup(GPSL1CA(),   gps_states),\n",
                "              galileo = SignalGroup(GalileoE1B(), galileo_states)))\n\n",
                "A single group needs no NamedTuple around it, and each group's ",
                "satellites may be a `Dictionary` keyed by PRN or a plain vector. ",
                "With Tracking.jl 9.1 or later loaded, `PositionVelocityTime.signal_groups(",
                "track_state, decoders)` builds them from a whole `TrackState`.",
            ),
        ),
    )
end

"""
    PositionVelocityTime.SUPPORTED_TIME_SYSTEMS

The GNSS time systems this package solves with, `(GPST(), GST(), BDT())`, in the order of
each [`SatelliteMeasurement`](@ref)'s `broadcast_time_offsets` tuple.
"""
const SUPPORTED_TIME_SYSTEMS = (GPST(), GST(), BDT())

"""
    PositionVelocityTime.SupportedTimeSystem

`Union{GPST,GST,BDT}`: the GNSS time systems a satellite can be on, which are exactly
[`SUPPORTED_TIME_SYSTEMS`](@ref). The type of [`SatelliteMeasurement`](@ref)'s
`time_system` field, and of the keys and the reference system of a [`PVTSolution`](@ref).

A closed `Union` rather than the abstract `GNSSSignals.TimeSystem`, so every call on a
time system splits into static branches instead of a dynamic dispatch, which a
`juliac --trim` build rejects, and the row stays one concrete type. A time system added
later joins `SUPPORTED_TIME_SYSTEMS`, and this follows.
"""
const SupportedTimeSystem = Union{map(typeof, SUPPORTED_TIME_SYSTEMS)...}

# Position of a time system in `SUPPORTED_TIME_SYSTEMS`, or `nothing`; a constant
# whenever the argument type is known.
time_offset_index(target::GNSSSignals.TimeSystem) =
    findfirst(==(target), SUPPORTED_TIME_SYSTEMS)

"""
    BroadcastTimeOffset

One constellation's broadcast offset toward one target GNSS time system, flattened into
an `isbits` value: every [`SatelliteMeasurement`](@ref) carries one per entry of
[`SUPPORTED_TIME_SYSTEMS`](@ref). It keeps the polynomial, so
[`calc_steering_offset`](@ref) can evaluate it at any time of week.

# Fields
- `available::Bool`: the flattened [`time_offset_available`](@ref). When `false` the
  coefficients are `NaN`, so the steering offset is too.
- `A_0`, `A_1`, `A_2::Float64`: the polynomial's constant, rate and acceleration terms.
- `t_0::Float64`: the reference time of week, or `NaN` where the message carries none
  (BeiDou D1/D2, whose argument is the time of week itself).
- `week_term::Float64`: `SECONDS_PER_WEEK * (own week − WN_0)`; `0.0` where `t_0` is `NaN`.
- `leap_second_offset_to_target::Float64`:
  `leap_second_offset_to_gpst(own) − leap_second_offset_to_gpst(target)`, which
  [`calc_steering_offset`](@ref) takes back out.
"""
struct BroadcastTimeOffset
    available::Bool
    A_0::Float64
    A_1::Float64
    A_2::Float64
    t_0::Float64
    week_term::Float64
    leap_second_offset_to_target::Float64
end

# No offset toward that system. NaN coefficients make a misuse propagate a NaN; a NaN
# `t_0` alone would not, as it marks BeiDou D1/D2's missing reference epoch.
const NO_TIME_OFFSET = BroadcastTimeOffset(false, NaN, NaN, NaN, NaN, 0.0, NaN)

"""
    SatelliteMeasurement

One satellite's contribution to one epoch: the concrete, parameter-free row the solver
works on. [`collect_measurements`](@ref) has already evaluated everything that depends on
the message family, the decoder or the signal, so a mixed-constellation epoch is one
`Vector{SatelliteMeasurement}`.

# Fields
- `prn::Int`: the satellite's PRN, unique only within its constellation.
- `signal_id::Symbol`: `get_signal_id` of the ranging signal (`:GPSL1CA`, …), half of the
  `(signal, PRN)` key of `PVTSolution.satellites`.
- `band_id::Symbol`: `get_band_id` of the ranging signal (`:L1`, `:L5`, …), which groups
  the inter-frequency biases.
- `time_system::SupportedTimeSystem`: `get_time_system` of the ranging signal, which
  groups the receiver clock biases.
- `leap_second_offset_to_gpst::Float64`: [`leap_second_offset_to_gpst`](@ref) of
  `time_system`, precomputed.
- `transmit_time::Float64`: the transmit time ([`calc_transmit_time`](@ref)) as a
  seconds-of-week count **on this satellite's own system's scale**.
- `position::SVector{3,Float64}`: satellite ECEF position at `transmit_time` (m).
- `velocity::SVector{3,Float64}`: satellite ECEF velocity at `transmit_time` (m/s).
- `carrier_doppler::Float64`: the measured carrier Doppler (Hz, unitless).
- `center_frequency::Float64`: the ranging signal's carrier (Hz, unitless) — the
  Doppler wavelength, and the `1/f²` of every ionospheric model.
- `satellite_clock_drift::Float64`: [`calc_satellite_clock_drift`](@ref) at
  `transmit_time` (s/s).
- `week::Int`: the absolute week number of `time_system`, with the GPS L1 C/A
  1024-week rollover already resolved against `approximate_year`.
- `system_start_epoch::TAITime`: the TAI epoch of that system's week 0
  ([`week_start_epoch`](@ref)).
- `broadcast_time_offsets::NTuple{3,BroadcastTimeOffset}`: this satellite's broadcast
  offsets toward [`SUPPORTED_TIME_SYSTEMS`](@ref), in that order.
"""
struct SatelliteMeasurement
    prn::Int
    signal_id::Symbol
    band_id::Symbol
    time_system::SupportedTimeSystem
    leap_second_offset_to_gpst::Float64
    transmit_time::Float64
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
    carrier_doppler::Float64
    center_frequency::Float64
    satellite_clock_drift::Float64
    week::Int
    system_start_epoch::TAITime
    broadcast_time_offsets::NTuple{3,BroadcastTimeOffset}
end

# One satellite state → one flat row: every decoder- and signal-typed call of the solve.
function satellite_measurement(state::SatelliteState, approximate_year::Integer)
    decoder = state.decoder
    signal = state.signal
    transmit_time = calc_transmit_time(state)
    position_and_velocity = calc_satellite_position_and_velocity(decoder, transmit_time)
    time_system = get_time_system(signal)
    leap_second_offset = leap_second_offset_to_gpst(time_system)
    week = get_week(decoder; approximate_year)
    SatelliteMeasurement(
        decoder.prn,
        get_signal_id(signal),
        get_band_id(signal),
        time_system,
        leap_second_offset,
        transmit_time,
        position_and_velocity.position,
        position_and_velocity.velocity,
        ustrip(Hz, state.carrier_doppler),
        ustrip(Hz, get_center_frequency(signal)),
        calc_satellite_clock_drift(decoder, transmit_time),
        week,
        system_start_epoch(signal),
        broadcast_time_offsets(decoder, week, leap_second_offset),
    )
end

"""
    collect_measurements(groups; approximate_year = year(now(UTC)))
        -> (Vector{SatelliteMeasurement}, ionospheric_parameters)

Flatten the [`SignalGroups`](@ref) of one epoch into the solver's input: the
[`SatelliteMeasurement`](@ref) rows, and the ionospheric parameters
[`select_ionospheric_parameters`](@ref) would pick from the same satellites.

A satellite is kept only if its navigation data is complete (checked first, which
guarantees its health bit is decoded) and it reports healthy. The rows are in group
order × within-group order (see [`calc_pvt`](@ref), also for `approximate_year`).
"""
function collect_measurements(groups; approximate_year::Integer = year(now(UTC)))
    measurements = SatelliteMeasurement[]
    # Reserve the unfiltered count (a row is ~310 bytes) rather than regrow.
    sizehint!(measurements, sum(group -> length(group.satellites),
        _normalize_signal_groups(groups); init = 0))
    measurements, collect_measurements!(measurements, groups; approximate_year)
end

"""
    collect_measurements!(measurements::Vector{SatelliteMeasurement}, groups;
                          approximate_year = year(now(UTC)))
        -> ionospheric_parameters

[`collect_measurements`](@ref) into `measurements` (emptied first), returning only the
ionospheric parameters, an `isbits` `Union` returned unboxed. Allocates nothing once its
buffers have grown. It is the pass [`calc_pvt!`](@ref) runs.
"""
collect_measurements!(
    measurements::Vector{SatelliteMeasurement},
    groups;
    approximate_year::Integer = year(now(UTC)),
) = _collect_measurements!(measurements, groups, approximate_year)

# The collection pass proper. Positional, for the solver's own call.
function _collect_measurements!(measurements, groups, approximate_year::Integer)
    normalised = _normalize_signal_groups(groups)
    empty!(measurements)
    # In group order. Over the `Tuple`, not the `NamedTuple`, whose `map` would not
    # specialise on the closure and dispatch each group dynamically.
    candidates = map(values(normalised)) do group
        collect_group!(measurements, group, approximate_year)
    end
    select_from_ionospheric_candidates(merge_all_ionospheric_candidates(candidates))
end

# One group's satellites, appended to `measurements`; returns its ionospheric candidates.
# Specialises on the group's `SatelliteState` type, so every call in the loop is static.
function collect_group!(measurements, group::SignalGroup, approximate_year::Integer)
    candidates = NO_IONOSPHERIC_CANDIDATES
    for state in group.satellites
        is_decoding_completed_for_positioning(state.decoder) || continue
        is_sat_healthy(state.decoder) || continue
        push!(measurements, satellite_measurement(state, approximate_year))
        candidates = update_ionospheric_candidates(candidates, state.decoder)
    end
    candidates
end
