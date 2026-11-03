# The GNSS time systems: the broadcast offsets between them.

"""
    BroadcastTimeOffsets

The broadcast offsets between the GNSS time systems of one epoch, gathered by
[`collect_measurements`](@ref) from all its observations and read with
[`broadcast_time_offset`](@ref). Toward GPS Time that is Galileo's GGTO and BeiDou's BGTO;
toward Galileo System Time, BeiDou's BGTO variant and the GGTO GPS broadcasts on
CNAV/CNAV-2. An offset can put a constellation's observations on another system's clock
(see [`decide_bias_layout`](@ref)).

An offset belongs to a pair of time systems, not to the satellite that broadcasts it, so
there is one per pair and direction: the steering residual `Δt = t_own − t_target` (s)
between the two atomic scales, as GNSSDecoder's `get_broadcast_time_offset` defines and
signs it, without the whole seconds between their counts, which GNSSSignals'
`get_time_system_offset` applies to the transmit times. Each is evaluated at the
broadcasting satellite's transmit time, on that satellite's own scale. Both directions of
a pair are kept as broadcast, Galileo's GST → GPST beside GPS's GPST → GST.

Where several satellites broadcast one direction, the one with the latest reference epoch
(GNSSDecoder's `get_broadcast_time_offset_epoch`) wins, so the newer set wins while one is
being uploaded. An offset with a reference epoch beats one without (BeiDou D1/D2 broadcasts
none), and among equals the first observation stays.

`BroadcastTimeOffsets()` holds none. Immutable and allocation-free.
"""
struct BroadcastTimeOffsets
    # `[own, target]` by `time_system_index`, `NaN` where none is broadcast.
    offsets::SMatrix{3,3,Float64,9}
    # Each offset's reference epoch (s since TAITime's origin), `-Inf` where it has none.
    epochs::SMatrix{3,3,Float64,9}
end
BroadcastTimeOffsets() = BroadcastTimeOffsets(
    fill(NaN, SMatrix{3,3,Float64,9}), fill(-Inf, SMatrix{3,3,Float64,9}))

# A time system's row and column in `BroadcastTimeOffsets`, in the order of
# `CANDIDATE_HUB_SYSTEMS`. Identity tests rather than methods, so a call on the
# `SupportedTimeSystem` field of a prepared observation needs no dispatch.
time_system_index(system::SupportedTimeSystem) =
    system === GPST() ? 1 : system === GST() ? 2 : 3

"""
    broadcast_time_offset(offsets::BroadcastTimeOffsets, own, target) -> Union{Nothing,Float64}

The broadcast offset `Δt = t_own − t_target` (s) from time system `own` to `target` (see
[`BroadcastTimeOffsets`](@ref)), or `nothing` where none was broadcast, `own` is
`target`, or either is outside [`CANDIDATE_HUB_SYSTEMS`](@ref), as GNSSDecoder's
`get_broadcast_time_offset` answers. The reverse direction is not substituted.
"""
function broadcast_time_offset(
    offsets::BroadcastTimeOffsets,
    own::SupportedTimeSystem,
    target::SupportedTimeSystem,
)
    broadcast_time_offset_at(offsets, time_system_index(own), time_system_index(target))
end

# `broadcast_time_offset` by `time_system_index`. A caller holding both systems as
# `SupportedTimeSystem` elements indexes first: a call on two of them spans 3 × 3
# signatures, more than Julia 1.10 union-splits, so it would dispatch dynamically and box
# its result.
function broadcast_time_offset_at(offsets::BroadcastTimeOffsets, own::Int, target::Int)
    Δt = offsets.offsets[own, target]
    isnan(Δt) ? nothing : Δt
end
broadcast_time_offset(::BroadcastTimeOffsets, ::GNSSSignals.TimeSystem,
    ::GNSSSignals.TimeSystem) = nothing

# `offsets` with `decoder`'s offsets toward every supported time system, evaluated at
# `transmit_time` on `own`, the decoder's system; each only where it is fresher than the
# one held (see `BroadcastTimeOffsets`). Unrolled, so each target is a constant.
function update_time_offsets(offsets::BroadcastTimeOffsets, decoder, own, transmit_time)
    offsets = update_time_offset(offsets, decoder, own, GPST(), transmit_time)
    offsets = update_time_offset(offsets, decoder, own, GST(), transmit_time)
    update_time_offset(offsets, decoder, own, BDT(), transmit_time)
end

function update_time_offset(offsets::BroadcastTimeOffsets, decoder, own, target,
    transmit_time)
    own === target && return offsets
    Δt = get_broadcast_time_offset(decoder, target, transmit_time)
    isnothing(Δt) && return offsets
    reference = get_broadcast_time_offset_epoch(decoder, target)
    epoch = isnothing(reference) ? -Inf : reference - TAITime(0, 0.0)
    i, j = time_system_index(own), time_system_index(target)
    isnan(offsets.offsets[i, j]) || epoch > offsets.epochs[i, j] || return offsets
    BroadcastTimeOffsets(setindex(offsets.offsets, Δt, i, j),
        setindex(offsets.epochs, epoch, i, j))
end

# The two tables of consecutive collections as one, cell by cell under the same rule:
# `later`'s offset replaces `earlier`'s only where it is fresher.
function merge_time_offsets(earlier::BroadcastTimeOffsets, later::BroadcastTimeOffsets)
    take = map((Δt, epoch, later_Δt, later_epoch) ->
            !isnan(later_Δt) && (isnan(Δt) || later_epoch > epoch),
        earlier.offsets, earlier.epochs, later.offsets, later.epochs)
    BroadcastTimeOffsets(map(ifelse, take, later.offsets, earlier.offsets),
        map(ifelse, take, later.epochs, earlier.epochs))
end
