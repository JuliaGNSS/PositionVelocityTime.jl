# The collection pass: satellite *states* in, flat satellite *measurements* out.
#
# Everything that depends on a satellite's navigation-data type — the decoder, the
# ranging signal, the propagator, the clock polynomial, the broadcast time offsets —
# is evaluated here, once per satellite, and reduced to the single concrete
# [`SatelliteMeasurement`](@ref) the solver works on. This is the function
# barrier of `calc_pvt`: everything upstream of it specialises on the group shape
# (and is small), everything downstream compiles exactly once no matter which
# constellations the epoch mixes.

"""
    SatelliteState{CP<:Real,D<:GNSSDecoderState,S<:AbstractGNSSSignal}

Combines the GNSS decoder state with code and carrier phase measurements for a single satellite.

# Fields
- `decoder::GNSSDecoderState`: GNSS decoder state containing decoded navigation data
- `signal::AbstractGNSSSignal`: the ranging signal the code phase and carrier Doppler
  were measured on (e.g., `GPSL1CA()`, `GalileoE1C()`). It may be a dataless pilot,
  whose navigation data `decoder` decoded from the data component of the same band.
- `code_phase::CP`: Code phase measurement in chips
- `carrier_doppler`: Carrier Doppler frequency in Hz
- `carrier_phase::CP`: Carrier phase measurement in radians, matching
  `Tracking.get_carrier_phase` (default: `0.0`)

# Constructors
    SatelliteState(; decoder, signal, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, signal, sat_state)

The second constructor extracts code phase, carrier Doppler, and carrier phase from a
`Tracking` satellite state (`Tracking.TrackedSat`). Tracking.jl defines it (from
Tracking 9.1 on), so it is available once `Tracking` is loaded; this package does not
depend on Tracking.
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

A per-epoch measurement snapshot of the satellites tracked on one ranging `signal`.

Deliberately the same name as `Tracking.SignalGroup`, and deliberately **not
exported**: `using Tracking, PositionVelocityTime` does not clash, call sites read
`PositionVelocityTime.SignalGroup(...)`, and Julia prints the unexported type fully
qualified so an error message distinguishes the two by itself. The correspondence is
exact — Tracking's group holds the *loop state* of a set of satellites, this one holds
one epoch's *measurements* of the same satellites, with the same group keys and, when
built through Tracking.jl, the same dictionary keys.

Groups are the unit of type stability: within a group every satellite shares one
concrete `SatelliteState` type, so the per-satellite calls of
[`collect_measurements`](@ref) are statically dispatched. Pooling several groups into
one vector — what `calc_pvt` used to take — is what made them dynamic.

# Fields
- `signal::AbstractGNSSSignal`: the ranging signal of every satellite in the group.
  Carried on the group itself rather than read off a satellite, so the group's signal
  is known even when it holds none.
- `satellites`: the satellites, as a `Dictionary{Int,<:SatelliteState}` keyed by PRN
  (mirroring Tracking, and structurally preventing a duplicate PRN within the group)
  or as an `AbstractVector{<:SatelliteState}` (so a consumer can reuse one buffer per
  group). Iteration order is significant — see [`calc_pvt`](@ref).

Every satellite's own `state.signal` must be `signal`. This is not checked per epoch;
a group built by Tracking.jl ([`signal_groups`](@ref)) satisfies it by
construction.
"""
struct SignalGroup{S<:AbstractGNSSSignal,C}
    signal::S
    satellites::C
end

"""
    PositionVelocityTime.SignalGroups{N}

Type alias for the storage shape [`calc_pvt`](@ref) takes: a `NamedTuple` of `N`
[`SignalGroup`](@ref)s. Mirrors `Tracking.SignalGroups`, so a receiver can carry one
set of group keys through tracking, decoding and the fix.
"""
const SignalGroups{N} = NamedTuple{<:Any,<:NTuple{N,SignalGroup}}

# Normalization, mirroring `Tracking._normalize_signal_groups`: a NamedTuple passes
# through, a single group becomes the `:default` group (so a one-constellation solve
# stays a one-liner), and a bare tuple is numbered by position — `group1`, `group2`, …
# Numbered rather than named after the signals because group order is significant
# (see `calc_pvt`) and because two groups may legitimately share one signal id.
_normalize_signal_groups(groups::NamedTuple{<:Any,<:Tuple{Vararg{SignalGroup}}}) = groups
_normalize_signal_groups(group::SignalGroup) = (default = group,)
@generated function _normalize_signal_groups(groups::Tuple{Vararg{SignalGroup}})
    names = ntuple(i -> Symbol(:group, i), length(groups.parameters))
    :(NamedTuple{$names}(groups))
end

"""
    PositionVelocityTime.signal_groups(track_state, decoders) -> NamedTuple of SignalGroup

Build one epoch's [`SignalGroups`](@ref) from a receiver's own per-signal state. The
method that does this is defined by Tracking.jl (from Tracking 9.1 on) and takes a
`Tracking.TrackState` together with that receiver's decoders; see it for the details.

There is deliberately **no** method taking a flat `Vector{SatelliteState}`. Groups are
not a wrapper around the old input — they are the shape a receiver already has, and
converting a pooled vector back into them at every epoch would reintroduce, in the
conversion, exactly the dynamic dispatch the grouping exists to remove. Build the groups
where the satellites are tracked, and hand the same groups to [`calc_pvt`](@ref).
"""
function signal_groups end

# A flat vector of satellite states is what `calc_pvt` took before 6.0. It has no
# meaningful normalization — the signals it pools are a runtime property, so any
# grouping derived here would be inference-blind — and silently accepting it would hide
# that. Say so instead, and point at what to build. The same directions answer anything
# else that is not signal groups, above all the NamedTuple of bare state vectors
# `(gps = gps_states, …)` a 5.x call site is most likely to be migrated to first, which
# would otherwise surface as a `FieldError` deep inside the collection pass.
_normalize_signal_groups(states::AbstractVector{<:SatelliteState}) =
    throw_not_signal_groups("not a vector of `SatelliteState`s")
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
    PositionVelocityTime.CANDIDATE_HUB_SYSTEMS

The GNSS time systems a clock collapse can be anchored on, in the fixed order
[`decide_bias_layout`](@ref) tries them: `(GPST(), GST(), BDT())`. Every
[`SatelliteMeasurement`](@ref) carries its broadcast time offset toward each.
"""
const CANDIDATE_HUB_SYSTEMS = (GPST(), GST(), BDT())

"""
    PositionVelocityTime.SupportedTimeSystem

`Union{GPST,GST,BDT}`: the GNSS time systems a satellite can be on, which are exactly
[`CANDIDATE_HUB_SYSTEMS`](@ref) — every system GNSSSignals defines is a candidate hub.
The type of [`SatelliteMeasurement`](@ref)'s `time_system` field, and of the keys and
the reference system of a [`PVTSolution`](@ref).

A closed `Union` rather than the abstract `GNSSSignals.TimeSystem` so that every call on
a time system — hashing it into `inter_system_biases`, `get_leap_second_offset`, a
hub's broadcast offset — is split into statically dispatched branches instead of being a
dynamic dispatch, which a `juliac --trim` build rejects. The measurement it types stays one
concrete type, so the solver still compiles once for every constellation mix. A time
system GNSSSignals adds later has to be added to `CANDIDATE_HUB_SYSTEMS` anyway (a
measurement carries a field for the offset toward it, which
[`broadcast_time_offset`](@ref) returns), and this follows.
"""
const SupportedTimeSystem = Union{map(typeof, CANDIDATE_HUB_SYSTEMS)...}

"""
    SatelliteMeasurement

One satellite's contribution to one epoch, reduced to a single concrete type with no
type parameters — what the whole PVT solver works on.

Everything that depends on the navigation-message family, the decoder or the ranging
signal has already been evaluated by [`collect_measurements`](@ref): the transmit
time, the propagated orbit, the clock drift, the group-delay-corrected observables,
the classification keys and the broadcast time offsets. What remains is plain data —
numbers, `Symbol`s and the `time_system` singleton — stored inline in a vector, so a
mixed-constellation epoch is one `Vector{SatelliteMeasurement}` and the solver
compiles once for every mix.

# Fields
- `prn::Int`: the satellite's PRN, unique only within its constellation.
- `signal_id::Symbol`: `get_signal_id` of the ranging signal (`:GPSL1CA`, …), the
  per-signal half of the `(signal, PRN)` key of `PVTSolution.satellites`.
- `band_id::Symbol`: `get_band_id` of the ranging signal (`:L1`, `:L5`, …), the
  grouping key of the receiver inter-frequency biases.
- `time_system::SupportedTimeSystem`: `get_time_system` of the ranging signal
  (`GPST()`, `GST()`, `BDT()`), the grouping key of the receiver clock biases. A small
  `Union` of singletons (see [`SupportedTimeSystem`](@ref)), so the measurement stays concrete
  and every call on the field splits statically.
- `transmit_time::Float64`: the corrected transmit time ([`calc_transmit_time`](@ref)) as a
  seconds-of-week count **on this satellite's own system's scale**.
- `position::SVector{3,Float64}`: satellite ECEF position at `transmit_time` (m).
- `velocity::SVector{3,Float64}`: satellite ECEF velocity at `transmit_time` (m/s).
- `carrier_doppler::Float64`: the measured carrier Doppler (Hz, unitless).
- `center_frequency::Float64`: the ranging signal's carrier (Hz, unitless) — the
  Doppler wavelength, and the `1/f²` of every ionospheric model.
- `clock_drift::Float64`: [`calc_satellite_clock_drift`](@ref) at `transmit_time` (s/s).
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
    clock_drift::Float64
    week::Int
    time_offset_to_gpst::Union{Nothing,Float64}
    time_offset_to_gst::Union{Nothing,Float64}
    time_offset_to_bdt::Union{Nothing,Float64}
end

# One satellite state → one measurement. Every decoder- and signal-typed call of the
# whole solve happens here.
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
        something(get_week(decoder; approximate_year)),
        get_broadcast_time_offset(decoder, GPST(), transmit_time),
        get_broadcast_time_offset(decoder, GST(), transmit_time),
        get_broadcast_time_offset(decoder, BDT(), transmit_time),
    )
end

"""
    collect_measurements(groups; approximate_year = year(now(UTC)))
        -> (Vector{SatelliteMeasurement}, ionospheric_parameters)

Flatten the [`SignalGroups`](@ref) of one epoch into the solver's input: the
per-satellite [`SatelliteMeasurement`](@ref)s, and the single constellation-wide
ionospheric model [`select_ionospheric_parameters`](@ref) would pick for them (folded
into the same pass — both need exactly the satellites that survive the health gate).

A satellite is kept only if its full navigation-data set is decoded and it reports
healthy; completeness is checked first, which is what guarantees the health bit has
been decoded at all. Everything else about it — orbit, clock, observables,
classification keys, broadcast time offsets — is evaluated here, so nothing downstream
touches a decoder again.

This is the function barrier: it specialises on the group shape (and stays small),
while the solver it feeds sees one concrete element type for every constellation mix.
The order of the measurements is group order × within-group order, and it is significant — see
[`calc_pvt`](@ref).

`approximate_year` resolves the GPS L1 C/A 1024-week rollover; see [`calc_pvt`](@ref).
"""
function collect_measurements(groups; approximate_year::Integer = year(now(UTC)))
    measurements = SatelliteMeasurement[]
    # An epoch's satellite count is known before a single measurement is built, and one is
    # ~170 bytes, so reserving the whole vector up front saves the geometric regrowth
    # (five reallocations and copies for a dozen satellites). The reservation is the
    # unfiltered count — an unhealthy satellite leaves a little slack, never a regrowth.
    sizehint!(measurements, sum(group -> length(group.satellites),
        _normalize_signal_groups(groups); init = 0))
    collect_measurements!(measurements, groups; approximate_year)
end

"""
    collect_measurements!(measurements::Vector{SatelliteMeasurement}, groups;
                          approximate_year = year(now(UTC)))
        -> (measurements, ionospheric_parameters)

[`collect_measurements`](@ref) into an existing vector: `measurements` is emptied and
refilled with this epoch's measurements. What it held is overwritten.

The returned tuple pairs the (non-`isbits`) vector with a `Union` of coefficient sets,
which Julia returns boxed: one small allocation per call. A consumer that must not
allocate calls [`collect_measurement_rows!`](@ref), which returns the ionospheric
correction alone.
"""
collect_measurements!(
    measurements::Vector{SatelliteMeasurement},
    groups;
    approximate_year::Integer = year(now(UTC)),
) = measurements, collect_measurement_rows!(measurements, groups; approximate_year)

"""
    collect_measurement_rows!(measurements::Vector{SatelliteMeasurement}, groups;
                              approximate_year = year(now(UTC)))
        -> ionospheric_parameters

[`collect_measurements!`](@ref) without the vector in the result: `measurements` is
emptied and refilled with this epoch's measurements, and only the ionospheric parameters are
returned. That is a `Union` of `isbits` coefficient sets, which Julia returns unboxed, so
once `measurements` has held as many measurements before, this allocates nothing. It is the
collection pass [`calc_pvt!`](@ref) runs, for a consumer that runs its own estimator over
the same measurements in a loop that must not allocate.
"""
collect_measurement_rows!(
    measurements::Vector{SatelliteMeasurement},
    groups;
    approximate_year::Integer = year(now(UTC)),
) = _collect_measurements!(measurements, groups, approximate_year)

# The collection pass proper. Positional, for the solver's own call.
function _collect_measurements!(measurements, groups, approximate_year::Integer)
    normalized = _normalize_signal_groups(groups)
    empty!(measurements)
    # `map` over the groups visits them in order, which is what makes the measurement order
    # group order × within-group order. Over their `Tuple` rather than the `NamedTuple`
    # itself: `map(f, ::NamedTuple)` only passes `f` on, so Julia does not specialise it
    # on the closure, and every group's call would be a dynamic dispatch.
    candidates = map(values(normalized)) do group
        collect_group!(measurements, group, approximate_year)
    end
    select_from_ionospheric_candidates(merge_all_ionospheric_candidates(candidates))
end

# One group's satellites, appended to `measurements`; returns the group's ionospheric
# candidates. Specialises on the group's concrete `SatelliteState` type, so every call
# in the loop — the health gate, the propagator, the clock polynomial, the coefficient
# accessors — is statically dispatched.
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
