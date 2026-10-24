# The collection pass: satellite *states* in, flat satellite *measurements* out.
#
# Everything that depends on a satellite's navigation-data type — the decoder, the
# ranging signal, the propagator, the clock polynomial, the broadcast time offsets —
# is evaluated here, once per satellite, and reduced to the single concrete
# [`SatelliteMeasurement`](@ref) row the solver works on. This is the function
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
    SatelliteState(decoder, signal, tracked_satellite)

The second constructor extracts code phase, carrier Doppler, and carrier phase from
`tracked_satellite`, a `Tracking.TrackedSat`. Tracking.jl defines it (from
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
one vector would make them dynamic.

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

# Normalisation, mirroring `Tracking._normalize_signal_groups`: a NamedTuple passes
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

There is deliberately **no** method taking a flat `Vector{SatelliteState}`: groups are
the shape a receiver already has, and converting a pooled vector into them at every
epoch would reintroduce, in the conversion, exactly the dynamic dispatch the grouping
exists to remove. Build the groups
where the satellites are tracked, and hand the same groups to [`calc_pvt`](@ref).
"""
function signal_groups end

# Anything that is not signal groups is refused at the entry point, with directions to
# what to build, rather than surfacing as a `FieldError` deep inside the collection
# pass: a flat vector of satellite states has no meaningful normalisation (the signals
# it pools are a runtime property, so any grouping derived here would be
# inference-blind), nor does a NamedTuple of bare state vectors.
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

The GNSS time systems this package solves with: `(GPST(), GST(), BDT())`. Each
[`SatelliteMeasurement`](@ref)'s `broadcast_time_offsets` tuple holds its broadcast offsets
toward these, in this order, so an offset toward a system is a tuple index rather than a
dictionary lookup.
"""
const SUPPORTED_TIME_SYSTEMS = (GPST(), GST(), BDT())

"""
    PositionVelocityTime.SupportedTimeSystem

`Union{GPST,GST,BDT}`: the GNSS time systems a satellite can be on, which are exactly
[`SUPPORTED_TIME_SYSTEMS`](@ref). The type of [`SatelliteMeasurement`](@ref)'s
`time_system` field, and of the keys and the reference system of a [`PVTSolution`](@ref).

A closed `Union` rather than the abstract `GNSSSignals.TimeSystem` so that every call on
a time system — hashing it into `inter_system_biases`, `leap_second_offset_to_gpst`, an
offset's index — is split into statically dispatched branches instead of being a
dynamic dispatch, which a `juliac --trim` build rejects. The row it types stays one
concrete type, so the solver still compiles once for every constellation mix. A time
system GNSSSignals adds later has to be added to `SUPPORTED_TIME_SYSTEMS` anyway (the
`broadcast_time_offsets` tuple is indexed by it), and this follows.
"""
const SupportedTimeSystem = Union{map(typeof, SUPPORTED_TIME_SYSTEMS)...}

# Position of a time system in `SUPPORTED_TIME_SYSTEMS`, or `nothing` for a system no
# broadcast offset can target. Every element is a distinct singleton, so this folds to
# a constant whenever the argument type is known.
time_offset_index(target::GNSSSignals.TimeSystem) =
    findfirst(==(target), SUPPORTED_TIME_SYSTEMS)

"""
    BroadcastTimeOffset

One constellation's broadcast offset toward one target GNSS time system, flattened out
of the navigation message into an `isbits` row so that evaluating it costs no dispatch
and no decoder field access. One per entry of [`SUPPORTED_TIME_SYSTEMS`](@ref) is
carried by every [`SatelliteMeasurement`](@ref).

The polynomial *coefficients* are stored rather than an evaluated offset, so that
[`calc_steering_offset`](@ref) can evaluate it at any time of week.

# Fields
- `available::Bool`: whether this satellite's message carries a usable offset toward
  the target at all — the flattened [`time_offset_available`](@ref). Every other field
  is meaningless when `false`, and the coefficients are then `NaN`, so that
  [`calc_steering_offset`](@ref) of an unavailable offset is `NaN`.
- `A_0`, `A_1`, `A_2::Float64`: the broadcast polynomial's constant, rate and
  acceleration terms.
- `t_0::Float64`: the offset's reference time of week, or `NaN` where the message
  carries none (BeiDou D1/D2, whose two-term offset has no reference epoch — there the
  polynomial argument is the time of week itself).
- `week_term::Float64`: `SECONDS_PER_WEEK * (own week − WN_0)`, precomputed, so the
  polynomial argument needs no week lookup. `0.0` where `t_0` is `NaN`.
- `leap_second_offset_to_target::Float64`:
  `leap_second_offset_to_gpst(own) − leap_second_offset_to_gpst(target)`, the *defined*
  whole-second part `GNSSDecoder.GNSSTimeOffset.A_0` folds in and
  [`calc_steering_offset`](@ref) takes back out (see there).
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

# The "this satellite broadcasts nothing toward that system" row. The coefficients and
# the leap-second offset are NaN so that a misuse propagates a NaN rather than a
# plausible zero. A NaN `t_0` alone would not: `calc_steering_offset` reads it as the BeiDou
# D1/D2 "no reference epoch" marker and evaluates the polynomial on `t` instead.
const NO_TIME_OFFSET = BroadcastTimeOffset(false, NaN, NaN, NaN, NaN, 0.0, NaN)

"""
    SatelliteMeasurement

One satellite's contribution to one epoch, reduced to a single concrete type with no
type parameters — the row the whole PVT solver works on.

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
  `Union` of singletons (see [`SupportedTimeSystem`](@ref)), so the row stays concrete
  and every call on the field splits statically.
- `leap_second_offset_to_gpst::Float64`: [`leap_second_offset_to_gpst`](@ref) of
  `time_system`, precomputed, so the solver reads a number rather than branching on
  the time system per satellite.
- `transmit_time::Float64`: the transmit time ([`calc_transmit_time`](@ref)) as a
  seconds-of-week count **on this satellite's own system's scale**.
- `position::SVector{3,Float64}`: satellite ECEF position at `transmit_time` (m).
- `velocity::SVector{3,Float64}`: satellite ECEF velocity at `transmit_time` (m/s).
- `carrier_doppler::Float64`: the measured carrier Doppler (Hz, unitless).
- `center_frequency::Float64`: the ranging signal's carrier (Hz, unitless) — the
  Doppler wavelength, and the `1/f²` of every ionospheric model.
- `satellite_clock_drift::Float64`: [`calc_satellite_clock_drift`](@ref) at `transmit_time`
  (s/s).
- `week::Int`: the absolute week number of `time_system`, with the GPS L1 C/A
  1024-week rollover already resolved against `approximate_year`.
- `system_start_epoch::TAITime`: the TAI epoch of that system's week 0, which dates the
  reported fix and the day of year of the tropospheric model ([`day_of_year`](@ref)).
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

# One satellite state → one flat row. Every decoder- and signal-typed call of the
# whole solve happens here.
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
per-satellite [`SatelliteMeasurement`](@ref) rows, and the single constellation-wide
ionospheric parameters [`select_ionospheric_parameters`](@ref) would pick for them (folded
into the same pass — both need exactly the satellites that survive the health gate).

A satellite is kept only if its full navigation-data set is decoded and it reports
healthy; completeness is checked first, which is what guarantees the health bit has
been decoded at all. Everything else about it — orbit, clock, observables,
classification keys, broadcast time offsets — is evaluated here, so nothing downstream
touches a decoder again.

This is the function barrier: it specialises on the group shape (and stays small),
while the solver it feeds sees one concrete element type for every constellation mix.
The flat row order is group order × within-group order, and it is significant — see
[`calc_pvt`](@ref).

`approximate_year` resolves the GPS L1 C/A 1024-week rollover; see [`calc_pvt`](@ref).
"""
function collect_measurements(groups; approximate_year::Integer = year(now(UTC)))
    measurements = SatelliteMeasurement[]
    # An epoch's satellite count is known before a single row is built, and a row is
    # ~310 bytes, so reserving the whole vector up front saves the geometric regrowth
    # (five reallocations and copies for a dozen satellites). The reservation is the
    # unfiltered count — an unhealthy satellite leaves a little slack, never a regrowth.
    sizehint!(measurements, sum(group -> length(group.satellites),
        _normalize_signal_groups(groups); init = 0))
    measurements, collect_measurements!(measurements, groups; approximate_year)
end

"""
    collect_measurements!(measurements::Vector{SatelliteMeasurement}, groups;
                          approximate_year = year(now(UTC)))
        -> ionospheric_parameters

[`collect_measurements`](@ref) into an existing vector: `measurements` is emptied and
refilled with this epoch's rows, and only the ionospheric parameters are returned. They are
a `Union` of `isbits` coefficient sets, which Julia returns unboxed, so once `measurements`
has held as many rows before, this allocates nothing. It is the collection pass
[`calc_pvt!`](@ref) runs, for a consumer that runs its own estimator over the same rows in a
loop that must not allocate.
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
    # `map` over the groups visits them in order, which is what makes the flat row order
    # group order × within-group order. Over their `Tuple` rather than the `NamedTuple`
    # itself: `map(f, ::NamedTuple)` only passes `f` on, so Julia does not specialise it
    # on the closure, and every group's call would be a dynamic dispatch.
    candidates = map(values(normalised)) do group
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
