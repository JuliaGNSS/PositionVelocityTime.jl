# The collection pass: satellite *states* in, flat satellite *measurements* out.
#
# Everything that depends on a satellite's navigation-data type — decoder, ranging signal,
# propagator, clock polynomial, broadcast time offsets — is evaluated here, once per
# satellite, into the one concrete `SatelliteMeasurement` the solver works on. This
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
    SatelliteState(decoder, signal, tracked_sat)

The second reads the measurements off a `Tracking.TrackedSat`; Tracking.jl defines
it, so it is available once `Tracking` is loaded.
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

# A NamedTuple passes through, a single group becomes `:default`, and a bare tuple is
# numbered `group1`, `group2`, … (not named after the signals, which two groups may share).
_normalize_signal_groups(groups::NamedTuple{<:Any,<:Tuple{Vararg{SignalGroup}}}) = groups
_normalize_signal_groups(group::SignalGroup) = (default = group,)
@generated function _normalize_signal_groups(groups::Tuple{Vararg{SignalGroup}})
    names = ntuple(i -> Symbol(:group, i), length(groups.parameters))
    :(NamedTuple{$names}(groups))
end

"""
    PositionVelocityTime.signal_groups(track_state, decoders) -> NamedTuple of SignalGroup

Build one epoch's signal groups, a `NamedTuple` of [`SignalGroup`](@ref)s, from a
receiver's own per-signal state. Tracking.jl defines the method, taking a
`Tracking.TrackState` and the receiver's decoders.

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
                "With Tracking.jl loaded, `PositionVelocityTime.signal_groups(",
                "track_state, decoders)` builds them from a whole `TrackState`.",
            ),
        ),
    )
end

"""
    PositionVelocityTime.SUPPORTED_TIME_SYSTEMS

The GNSS time systems this package solves with, `(GPST(), GST(), BDT())`: those a
satellite can be on, and those a [`SatelliteMeasurement`](@ref) carries the broadcast time
offsets toward.
"""
const SUPPORTED_TIME_SYSTEMS = (GPST(), GST(), BDT())

"""
    PositionVelocityTime.SupportedTimeSystem

`Union{GPST,GST,BDT}`: the GNSS time systems a satellite can be on, which are exactly
[`SUPPORTED_TIME_SYSTEMS`](@ref). The type of [`SatelliteMeasurement`](@ref)'s
`time_system` field, and of the keys and the reference system of a [`PVTSolution`](@ref).

A closed `Union` rather than the abstract `GNSSSignals.TimeSystem`, so every call on a
time system splits into static branches instead of a dynamic dispatch, which a
`juliac --trim` build rejects, and a `SatelliteMeasurement` stays one concrete type. A time
system added later joins `SUPPORTED_TIME_SYSTEMS`, and this follows; a measurement then
carries a field for the offset toward it, which [`broadcast_time_offset`](@ref) returns.
"""
const SupportedTimeSystem = Union{map(typeof, SUPPORTED_TIME_SYSTEMS)...}

"""
    SatelliteMeasurement

One satellite's contribution to one epoch: the concrete, parameter-free record the solver
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
  1024-week rollover already resolved against `approximate_year`; its start is
  GNSSSignals' `get_week_start_epoch(measurement.time_system, measurement.week)`.
- `time_offset_to_gpst`, `time_offset_to_gst`, `time_offset_to_bdt`: this satellite's
  broadcast time offsets toward GPS, Galileo and BeiDou Time (s), evaluated at
  `transmit_time` by GNSSDecoder's `get_broadcast_time_offset`, or `nothing` where none
  is broadcast; [`broadcast_time_offset`](@ref) picks one by target. Fields of their own
  rather than a tuple, whose `Union` elements would be stored boxed.
"""
struct SatelliteMeasurement
    prn::Int
    signal_id::Symbol
    band_id::Symbol
    time_system::SupportedTimeSystem
    transmit_time::Float64
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
    carrier_doppler::Float64
    center_frequency::Float64
    satellite_clock_drift::Float64
    week::Int
    time_offset_to_gpst::Union{Nothing,Float64}
    time_offset_to_gst::Union{Nothing,Float64}
    time_offset_to_bdt::Union{Nothing,Float64}
end

# One satellite state → one measurement: every decoder- and signal-typed call of the solve.
function satellite_measurement(state::SatelliteState, approximate_year::Integer)
    decoder = state.decoder
    signal = state.signal
    ephemeris = something(get_broadcast_ephemeris(decoder))
    transmit_time = correct_clock(ephemeris, group_delay(decoder, signal),
        calc_uncorrected_transmit_time(state))
    pos_and_vel = calc_satellite_position_and_velocity(ephemeris, transmit_time)
    SatelliteMeasurement(
        decoder.prn,
        get_signal_id(signal),
        get_band_id(signal),
        get_time_system(signal),
        transmit_time,
        pos_and_vel.position,
        pos_and_vel.velocity,
        ustrip(Hz, state.carrier_doppler),
        ustrip(Hz, get_center_frequency(signal)),
        calc_satellite_clock_drift(ephemeris, transmit_time),
        something(get_broadcast_week(decoder; approximate_year)),
        get_broadcast_time_offset(decoder, GPST(), transmit_time),
        get_broadcast_time_offset(decoder, GST(), transmit_time),
        get_broadcast_time_offset(decoder, BDT(), transmit_time),
    )
end

"""
    collect_measurements(groups; approximate_year = year(now(UTC)))
        -> (Vector{SatelliteMeasurement}, ionospheric_parameters)

Flatten the signal groups of one epoch (see [`calc_pvt`](@ref)) into the solver's input: the
[`SatelliteMeasurement`](@ref)s, and the ionospheric parameters
[`select_ionospheric_parameters`](@ref) would pick from the same satellites.

A satellite is kept only if its navigation data is complete (checked first, which
guarantees its health bit is decoded) and it reports healthy. The measurements are in
group order × within-group order (see [`calc_pvt`](@ref), also for `approximate_year`).
"""
function collect_measurements(groups; approximate_year::Integer = year(now(UTC)))
    measurements = SatelliteMeasurement[]
    # Reserve the unfiltered count (a measurement is ~170 bytes) rather than regrow.
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
    normalized = _normalize_signal_groups(groups)
    empty!(measurements)
    # In group order. Over the `Tuple`, not the `NamedTuple`, whose `map` would not
    # specialise on the closure and dispatch each group dynamically.
    candidates = map(values(normalized)) do group
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
