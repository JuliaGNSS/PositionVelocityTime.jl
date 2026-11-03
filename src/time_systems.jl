# The GNSS time systems: broadcast offsets between them, and the whole seconds between
# their counts.

"""
    BroadcastInterSystemTimeOffsets

The broadcast inter-system time offsets of one epoch, gathered by
[`prepare_observations`](@ref) from all its observations and read with
[`broadcast_inter_system_time_offset`](@ref). Toward GPS Time that is Galileo's GGTO and
BeiDou's BGTO; toward Galileo System Time, BeiDou's BGTO variant and the GGTO GPS
broadcasts on CNAV/CNAV-2. An offset can tie two receiver clocks (see
[`decide_bias_layout`](@ref)).

An offset belongs to a pair of time systems, not to the satellite that broadcasts it, so
there is one per pair and direction: the steering residual `Δt = t_own − t_target` (s)
between the two atomic scales, as GNSSDecoder's `get_broadcast_inter_system_time_offset`
defines and signs it, without the whole seconds between their counts, which GNSSSignals'
`get_time_system_offset` applies to the transmit times. Each is evaluated at the
broadcasting satellite's transmit time, on that satellite's own scale. Both directions of
a pair are kept as broadcast, Galileo's GST → GPST beside GPS's GPST → GST: which to trust
is the constraint's choice, not the table's.

Where several satellites broadcast one direction, the one with the latest reference epoch
(GNSSDecoder's `get_broadcast_inter_system_time_offset_epoch`) wins, so the newer set wins
while one is being uploaded. An offset with a reference epoch beats one without (BeiDou
D1/D2 broadcasts none), and among equals the first observation stays.

`BroadcastInterSystemTimeOffsets()` holds none. Immutable and allocation-free.
"""
struct BroadcastInterSystemTimeOffsets
    # `[own, target]` by `time_system_index`, `NaN` where none is broadcast.
    offsets::SMatrix{3,3,Float64,9}
    # Each offset's reference epoch (s since TAITime's origin), `-Inf` where it has none.
    epochs::SMatrix{3,3,Float64,9}
end
BroadcastInterSystemTimeOffsets() = BroadcastInterSystemTimeOffsets(
    fill(NaN, SMatrix{3,3,Float64,9}), fill(-Inf, SMatrix{3,3,Float64,9}))

# A time system's row and column in `BroadcastInterSystemTimeOffsets`, in the order of
# `SUPPORTED_TIME_SYSTEMS`. Identity tests rather than methods, so a call on the
# `SupportedTimeSystem` field of a prepared observation needs no dispatch.
time_system_index(system::SupportedTimeSystem) =
    system === GPST() ? 1 : system === GST() ? 2 : 3

"""
    broadcast_inter_system_time_offset(offsets::BroadcastInterSystemTimeOffsets, own,
        target) -> Union{Nothing,Float64}

The broadcast offset `Δt = t_own − t_target` (s) from time system `own` to `target` (see
[`BroadcastInterSystemTimeOffsets`](@ref)), or `nothing` where none was broadcast, `own` is
`target`, or either is outside [`SUPPORTED_TIME_SYSTEMS`](@ref), as GNSSDecoder's
`get_broadcast_inter_system_time_offset` answers. The reverse direction is not
substituted.
"""
function broadcast_inter_system_time_offset(
    offsets::BroadcastInterSystemTimeOffsets,
    own::SupportedTimeSystem,
    target::SupportedTimeSystem,
)
    broadcast_inter_system_time_offset_at(offsets, time_system_index(own),
        time_system_index(target))
end

# `broadcast_inter_system_time_offset` by `time_system_index`. A caller holding both
# systems as `SupportedTimeSystem` elements indexes first: a call on two of them spans
# 3 × 3 signatures, more than Julia 1.10 union-splits, so it would dispatch dynamically
# and box its result.
function broadcast_inter_system_time_offset_at(offsets::BroadcastInterSystemTimeOffsets,
    own::Int, target::Int)
    Δt = offsets.offsets[own, target]
    isnan(Δt) ? nothing : Δt
end
broadcast_inter_system_time_offset(::BroadcastInterSystemTimeOffsets,
    ::GNSSSignals.TimeSystem, ::GNSSSignals.TimeSystem) = nothing

# `offsets` with `decoder`'s offsets toward every supported time system, evaluated at
# `transmit_time` on `own`, the decoder's system; each only where it is fresher than the
# one held (see `BroadcastInterSystemTimeOffsets`). Unrolled, so each target is a
# constant.
function update_inter_system_time_offsets(offsets::BroadcastInterSystemTimeOffsets,
    decoder, own, transmit_time)
    offsets = update_inter_system_time_offset(offsets, decoder, own, GPST(), transmit_time)
    offsets = update_inter_system_time_offset(offsets, decoder, own, GST(), transmit_time)
    update_inter_system_time_offset(offsets, decoder, own, BDT(), transmit_time)
end

function update_inter_system_time_offset(offsets::BroadcastInterSystemTimeOffsets,
    decoder, own, target, transmit_time)
    own === target && return offsets
    Δt = get_broadcast_inter_system_time_offset(decoder, target, transmit_time)
    isnothing(Δt) && return offsets
    reference = get_broadcast_inter_system_time_offset_epoch(decoder, target)
    epoch = isnothing(reference) ? -Inf : reference - TAITime(0, 0.0)
    i, j = time_system_index(own), time_system_index(target)
    isnan(offsets.offsets[i, j]) || epoch > offsets.epochs[i, j] || return offsets
    BroadcastInterSystemTimeOffsets(setindex(offsets.offsets, Δt, i, j),
        setindex(offsets.epochs, epoch, i, j))
end

# The two tables of consecutive collections as one, cell by cell under the same rule:
# `later`'s offset replaces `earlier`'s only where it is fresher.
function merge_inter_system_time_offsets(earlier::BroadcastInterSystemTimeOffsets,
    later::BroadcastInterSystemTimeOffsets)
    take = map((Δt, epoch, later_Δt, later_epoch) ->
            !isnan(later_Δt) && (isnan(Δt) || later_epoch > epoch),
        earlier.offsets, earlier.epochs, later.offsets, later.epochs)
    BroadcastInterSystemTimeOffsets(map(ifelse, take, later.offsets, earlier.offsets),
        map(ifelse, take, later.epochs, earlier.epochs))
end

"""
    calc_pseudorange(reference_time, transmit_time, correction = 0.0) -> Float64

The pseudorange (m) of a satellite that sent `transmit_time` (s, a time of week on its own
count) against `reference_time` (s, a time of week on the reference count): how far the
transmit time, moved onto the reference count by `correction`, trails the reference time,
folded across a week crossover (GNSSSignals' `fold_week_crossover`) and scaled by
[`SPEED_OF_LIGHT`](@ref).

`correction` is the seconds between the two counts, GNSSSignals'
`get_time_system_offset(own_system, reference_system)`: `14` from BeiDou Time onto GPS
Time, `0` between systems that count alike. It moves the pseudorange,
not the transmit time: the satellite is propagated from its transmit time against the
message's own `t_0e`, on the broadcasting system's scale, so shifting the transmit time
would move the 14 s into the ephemeris, about 55 km of along-track error at BeiDou MEO
velocities. A solve corrects toward its own reference system, whose count also dates its
epoch (`calc_pvt` toward the one of [`select_reference_system`](@ref)).

A crossover can fall between the two times: within the spread of transmit times, or, in a
mixed solve, between BeiDou times on the GPS count (`SOW + 14`, past 604800 for the first
14 s of a GPS week) and GPS times already wrapped. `reference_time` need not be a receive
time: a common one cancels into the receiver clock biases.
"""
calc_pseudorange(reference_time, transmit_time, correction = 0.0) =
    fold_week_crossover(reference_time - (correction + transmit_time)) * SPEED_OF_LIGHT
