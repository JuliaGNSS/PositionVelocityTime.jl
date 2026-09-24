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
    PositionVelocityTime.SignalGroup(signal, satellites)

A per-epoch measurement snapshot of the satellites tracked on one ranging `signal`.

Deliberately the same name as `Tracking.SignalGroup`, and deliberately **not
exported**: `using Tracking, PositionVelocityTime` does not clash, call sites read
`PositionVelocityTime.SignalGroup(...)`, and Julia prints the unexported type fully
qualified so an error message distinguishes the two by itself. The correspondence is
exact — Tracking's group holds the *loop state* of a set of satellites, this one holds
one epoch's *measurements* of the same satellites, with the same group keys and, when
built through the Tracking extension, the same dictionary keys.

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

Every satellite's own `state.system` must be `signal`. This is not checked per epoch;
a group built by the Tracking extension ([`signal_groups`](@ref)) satisfies it by
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
method that does this lives in the `Tracking` extension and takes a `Tracking.TrackState`
together with that receiver's decoders; see it for the details.

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
                "With `Tracking` loaded, `PositionVelocityTime.signal_groups(",
                "track_state, decoders)` builds them from a whole `TrackState`.",
            ),
        ),
    )
end

"""
    PositionVelocityTime.CANDIDATE_HUB_SYSTEMS

The GNSS time systems a clock collapse can be anchored on, in the fixed order
[`decide_bias_layout`](@ref) tries them: `(GPST(), GST(), BDT())`. Each
[`SatelliteMeasurement`](@ref)'s `time_offsets` tuple is indexed by this order, so a
hub's broadcast offset is a tuple index rather than a dictionary lookup.
"""
const CANDIDATE_HUB_SYSTEMS = (GPST(), GST(), BDT())

# Position of a time system in `CANDIDATE_HUB_SYSTEMS`, or `nothing` for a system no
# broadcast offset can target. Every element is a distinct singleton, so this folds to
# a constant whenever the argument type is known.
hub_system_index(target::GNSSSignals.TimeSystem) =
    findfirst(==(target), CANDIDATE_HUB_SYSTEMS)

"""
    BroadcastTimeOffset

One constellation's broadcast offset toward one target GNSS time system, flattened out
of the navigation message into an `isbits` row so that evaluating it costs no dispatch
and no decoder field access. One per entry of [`CANDIDATE_HUB_SYSTEMS`](@ref) is
carried by every [`SatelliteMeasurement`](@ref).

The polynomial *coefficients* are stored rather than an evaluated offset because
[`calc_hub_range_offsets`](@ref) evaluates one system-wide offset at **each
satellite's own transmit time**; an evaluated scalar per satellite would change that
(negligibly in metres, but visibly against the tolerances the steering tests pin).

# Fields
- `available::Bool`: whether this satellite's message carries a usable offset toward
  the target at all — the flattened [`time_offset_available`](@ref). Every other field
  is meaningless when `false`.
- `A_0`, `A_1`, `A_2::Float64`: the broadcast polynomial's constant, rate and
  acceleration terms.
- `t_0::Float64`: the offset's reference time of week, or `NaN` where the message
  carries none (BeiDou D1/D2, whose two-term offset has no reference epoch — there the
  polynomial argument is the time of week itself).
- `week_term::Float64`: `SECONDS_PER_WEEK * (own week − WN_0)`, precomputed, so the
  polynomial argument needs no week lookup. `0.0` where `t_0` is `NaN`.
- `count_anchor::Float64`: `time_scale_offset_to_gpst(own) −
  time_scale_offset_to_gpst(target)`, the *defined* whole-second part
  `GNSSDecoder.GNSSTimeOffset.A_0` folds in and [`calc_steering_offset`](@ref) takes
  back out (see there).
"""
struct BroadcastTimeOffset
    available::Bool
    A_0::Float64
    A_1::Float64
    A_2::Float64
    t_0::Float64
    week_term::Float64
    count_anchor::Float64
end

# The "this satellite broadcasts nothing toward that system" row. `t_0` is NaN so that
# a misuse of the coefficients propagates a NaN rather than a plausible zero.
const NO_TIME_OFFSET = BroadcastTimeOffset(false, 0.0, 0.0, 0.0, NaN, 0.0, 0.0)

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
  per-signal half of the `(signal, PRN)` key of `PVTSolution.sats`.
- `band_id::Symbol`: `get_band_id` of the ranging signal (`:L1`, `:L5`, …), the
  grouping key of the receiver inter-frequency biases.
- `time_system::GNSSSignals.TimeSystem`: `get_time_system` of the ranging signal
  (`GPST()`, `GST()`, `BDT()`), the grouping key of the receiver clock biases. An
  abstract field, but every instance is a singleton, so this is one pointer.
- `count_offset_to_gpst::Float64`: [`time_scale_offset_to_gpst`](@ref) of
  `time_system`, precomputed — reading it off the abstract `time_system` field would
  be a dynamic dispatch per satellite.
- `time::Float64`: the corrected transmit time ([`calc_corrected_time`](@ref)) as a
  seconds-of-week count **on this satellite's own system's scale**.
- `position::SVector{3,Float64}`: satellite ECEF position at `time` (m).
- `velocity::SVector{3,Float64}`: satellite ECEF velocity at `time` (m/s).
- `carrier_doppler::Float64`: the measured carrier Doppler (Hz, unitless).
- `center_frequency::Float64`: the ranging signal's carrier (Hz, unitless) — the
  Doppler wavelength, and the `1/f²` of every ionospheric model.
- `clock_drift::Float64`: [`calc_satellite_clock_drift`](@ref) at `time` (s/s).
- `week::Int`: the absolute week number of `time_system`, with the GPS L1 C/A
  1024-week rollover already resolved against `approximate_year`.
- `system_start_epoch::TAIEpoch{Float64}`: the TAI epoch of that system's week 0,
  which dates the reported fix.
- `system_start_time::DateTime`: the same origin as a calendar date, for
  [`day_of_year`](@ref).
- `time_offsets::NTuple{3,BroadcastTimeOffset}`: this satellite's broadcast offsets
  toward [`CANDIDATE_HUB_SYSTEMS`](@ref), in that order.
"""
struct SatelliteMeasurement
    prn::Int
    signal_id::Symbol
    band_id::Symbol
    time_system::GNSSSignals.TimeSystem
    count_offset_to_gpst::Float64
    time::Float64
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
    carrier_doppler::Float64
    center_frequency::Float64
    clock_drift::Float64
    week::Int
    system_start_epoch::TAIEpoch{Float64}
    system_start_time::DateTime
    time_offsets::NTuple{3,BroadcastTimeOffset}
end

"""
    time_offset_available(measurement::SatelliteMeasurement, target) -> Bool
    time_offset_available(decoder::GNSSDecoderState, target) -> Bool

Whether a usable broadcast offset from this satellite's own GNSS time system to
`target` is available. Such an offset lets that constellation's measurements be
expressed on the target system's clock, which makes a fix possible when the geometry
is too weak to estimate an independent clock bias for it — see
[`decide_bias_layout`](@ref), which collapses onto a hub system this way. Toward GPS
Time that is Galileo's GGTO and BeiDou's BGTO; toward Galileo System Time it is
BeiDou's BGTO variant and the GGTO GPS itself broadcasts on CNAV/CNAV-2.

The measurement form reads the flag [`collect_measurements`](@ref) already flattened
into `measurement.time_offsets` (equivalently `measurement.time_offsets[i].available`),
and is `false` for any `target` outside [`CANDIDATE_HUB_SYSTEMS`](@ref), which no
message broadcasts an offset toward.

The decoder form is the underlying predicate, and delegates to GNSSDecoder's
`get_time_offset`, which screens every way a signal can fail to have one: it
broadcasts none toward `target` (GPS L1 C/A and Galileo E6-B broadcast none at all),
it does but has not decoded one yet, the decoded one names a different target system,
the ICD's "not available" sentinel is set (a `GNSS_ID` of 0; Galileo's all-ones GGTO),
or — on Galileo — the reference week has not arrived yet, since word type 10 can be
decoded before the week number is. Always `false` when `target` is the decoder's own
system: the offset from a scale to itself is not a broadcast quantity.
"""
function time_offset_available(
    measurement::SatelliteMeasurement,
    target::GNSSSignals.TimeSystem,
)
    index = hub_system_index(target)
    isnothing(index) ? false : measurement.time_offsets[index].available
end

time_offset_available(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem,
) = !isnothing(get_time_offset(decoder, target))

"""
    calc_steering_offset(offset::BroadcastTimeOffset, t) -> Float64

The broadcast *steering* offset between a constellation's time scale and the target
scale `offset` was built against, in seconds, evaluated at the own-scale time of week
`t`. Subtract it to convert a transmit time to the target system's time. Only defined
where `offset.available` (see [`time_offset_available`](@ref)).

Tens of nanoseconds: it is the residual between two atomic scales, not the whole
difference between their counts. GNSSDecoder's `GNSSTimeOffset.A_0` folds in the
*defined* whole-second offset as well — so that `t_target = t_own − Δt` holds for the
seconds — and `offset.count_anchor` takes that part back out, because
[`calc_time_scale_offsets`](@ref) has already applied it to the transmit times before
they were differenced. Both sides are the same
[`time_scale_offset_to_gpst`](@ref) expression, so the two halves compose exactly
rather than approximately.

!!! note "The subtraction costs about seven digits of the residual"

    Recovering a ~1e-8 s residual by subtracting 14 s from a `Float64` leaves roughly
    1.8e-15 s of rounding — one ULP at 14. That is 0.5 µm of range, so it is
    irrelevant to a fix, but it is why the tests here compare the steering term with
    an explicit tolerance rather than the default `≈`.
"""
function calc_steering_offset(offset::BroadcastTimeOffset, t)
    # A `NaN` `t_0` marks a message that carries no reference epoch at all (BeiDou
    # D1/D2); there the polynomial argument is the time of week itself. Where it is
    # present, `week_term` already carries `SECONDS_PER_WEEK * (own week − WN_0)`.
    Δτ = isnan(offset.t_0) ? t : t - offset.t_0 + offset.week_term
    offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2 - offset.count_anchor
end

# The broadcast offsets toward every candidate hub, flattened. `map` over the
# heterogeneous `CANDIDATE_HUB_SYSTEMS` tuple specialises per element, so each
# `get_time_offset` call is statically dispatched on both the decoder and the target.
function broadcast_time_offsets(decoder, week, count_offset_to_gpst)
    map(CANDIDATE_HUB_SYSTEMS) do target
        offset = get_time_offset(decoder, target)
        isnothing(offset) && return NO_TIME_OFFSET
        BroadcastTimeOffset(
            true,
            offset.A_0,
            offset.A_1,
            offset.A_2,
            isnothing(offset.t_0) ? NaN : Float64(offset.t_0),
            isnothing(offset.WN_0) ? 0.0 : SECONDS_PER_WEEK * (week - offset.WN_0),
            count_offset_to_gpst - time_scale_offset_to_gpst(target),
        )
    end
end

"""
    broadcast_time_offset(decoder, target; approximate_year = year(now(UTC)))
        -> BroadcastTimeOffset

The single [`BroadcastTimeOffset`](@ref) from `decoder`'s own GNSS time system toward
`target`, as [`collect_measurements`](@ref) flattens it onto every
[`SatelliteMeasurement`](@ref) — the standalone form, for a consumer evaluating a
broadcast time offset outside a fix. Unavailable (`available == false`) for a `target`
outside [`CANDIDATE_HUB_SYSTEMS`](@ref), and for any message that carries no such
offset (see [`time_offset_available`](@ref)).
"""
function broadcast_time_offset(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem;
    approximate_year::Integer = year(now(UTC)),
)
    index = hub_system_index(target)
    isnothing(index) && return NO_TIME_OFFSET
    broadcast_time_offsets(
        decoder,
        get_week(decoder; approximate_year),
        time_scale_offset_to_gpst(get_time_system(decoder)),
    )[index]
end

# One satellite state → one flat row. Every decoder- and signal-typed call of the
# whole solve happens here.
function satellite_measurement(state::SatelliteState, approximate_year::Integer)
    decoder = state.decoder
    system = state.system
    time = calc_corrected_time(state)
    position_and_velocity = calc_satellite_position_and_velocity(decoder, time)
    time_system = get_time_system(system)
    count_offset_to_gpst = time_scale_offset_to_gpst(time_system)
    week = get_week(decoder; approximate_year)
    SatelliteMeasurement(
        decoder.prn,
        get_signal_id(system),
        get_band_id(system),
        time_system,
        count_offset_to_gpst,
        time,
        get_sat_position(position_and_velocity),
        get_sat_velocity(position_and_velocity),
        ustrip(Hz, state.carrier_doppler),
        ustrip(Hz, get_center_frequency(system)),
        calc_satellite_clock_drift(decoder, time),
        week,
        system_start_epoch(system),
        get_system_start_time(system),
        broadcast_time_offsets(decoder, week, count_offset_to_gpst),
    )
end

"""
    collect_measurements(groups; approximate_year = year(now(UTC)))
        -> (Vector{SatelliteMeasurement}, ionospheric_correction)

Flatten the [`SignalGroups`](@ref) of one epoch into the solver's input: the
per-satellite [`SatelliteMeasurement`](@ref) rows, and the single constellation-wide
ionospheric model [`select_ionospheric_correction`](@ref) would pick for them (folded
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
    normalized = _normalize_signal_groups(groups)
    measurements = SatelliteMeasurement[]
    # An epoch's satellite count is known before a single row is built, and a row is
    # ~320 bytes, so reserving the whole vector up front saves the geometric regrowth
    # (five reallocations and copies for a dozen satellites). The reservation is the
    # unfiltered count — an unhealthy satellite leaves a little slack, never a regrowth.
    sizehint!(measurements, sum(group -> length(group.satellites), normalized; init = 0))
    # `map` over the NamedTuple visits the groups in order, which is what makes the
    # flat row order group order × within-group order.
    candidates = map(normalized) do group
        collect_group!(measurements, group, approximate_year)
    end
    measurements,
    select_from_ionospheric_candidates(merge_all_ionospheric_candidates(values(candidates)))
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
