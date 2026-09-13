module PositionVelocityTimeTrackingExt

using PositionVelocityTime: PositionVelocityTime, SatelliteState
using GNSSDecoder: GNSSDecoderState
using GNSSSignals: AbstractGNSSSignal
using Dictionaries: Dictionary
using Tracking: Tracking, get_code_phase, get_carrier_doppler, get_carrier_phase

# Both packages define a type named `SignalGroup`, deliberately (see the docstring of
# `PositionVelocityTime.SignalGroup`), and this file is the one place where both are in
# scope. So every `SignalGroup` below is written out in full, on both sides, with no
# exceptions — an unqualified one here would resolve to Tracking's export.

function PositionVelocityTime.SatelliteState(
    decoder::GNSSDecoderState,
    system::AbstractGNSSSignal,
    sat_state::Tracking.TrackedSat,
)
    SatelliteState(
        decoder,
        system,
        get_code_phase(sat_state),
        get_carrier_doppler(sat_state),
        get_carrier_phase(sat_state),
    )
end

"""
    PositionVelocityTime.SignalGroup(group::Tracking.SignalGroup, decoders;
                                     signal = first(group.signals))

The per-epoch measurement snapshot of one `Tracking.SignalGroup`: each tracked
satellite's loop state paired with the navigation data decoded from it, as a
`PositionVelocityTime.SignalGroup` [`calc_pvt`](@ref) can solve.

`decoders` is a `Dictionary` (or any `getindex`/`haskey` mapping) from the group's own
satellite keys — PRNs, the keys of `group.satellites` — to `GNSSDecoderState`s. Those
keys are preserved: the result is keyed exactly as the tracking group is, so a receiver
can index either by the same identifier. A satellite with no entry in `decoders` is
skipped rather than erroring, since tracking a satellite starts long before its
navigation message is complete.

`signal` is the ranging signal the pseudorange is generated on, and defaults to the
group's first signal — both Tracking's own estimator-driver signal and GNSSReceiver's
`RANGING_SIGNAL_INDEX = 1`. Pass it explicitly to range on another signal of the group.
"""
function PositionVelocityTime.SignalGroup(
    group::Tracking.SignalGroup,
    decoders;
    signal::AbstractGNSSSignal = first(group.signals),
)
    # The keys are selected first and the states mapped over them, rather than both
    # being pushed in one loop: `map` over the kept keys gives a concretely-typed
    # vector of states (one decoder type, one signal, one `TrackedSat` type per
    # group), and a concretely-typed group is the whole point of grouping.
    kept = [key for key in keys(group.satellites) if haskey(decoders, key)]
    states = map(key -> SatelliteState(decoders[key], signal, group.satellites[key]), kept)
    PositionVelocityTime.SignalGroup(signal, Dictionary(kept, states))
end

"""
    PositionVelocityTime.signal_groups(track_state::Tracking.TrackState, decoders)
        -> NamedTuple of PositionVelocityTime.SignalGroup

One epoch's measurements from a whole `Tracking.TrackState`, ready for
[`calc_pvt`](@ref). `decoders` is a `NamedTuple` with the same group keys as
`track_state.groups`, each holding that group's PRN → `GNSSDecoderState` mapping.

Both levels of key are preserved — the group names and, within each group, the
satellite keys — so the fix is indexed exactly as the tracking state is. Each group's
ranging signal is its first signal; see
[`PositionVelocityTime.SignalGroup`](@ref)`(::Tracking.SignalGroup, decoders)`.

Unlike `signal_groups` on a flat vector of states, this is fully inferable: the group
names and every group's concrete type come from `track_state`'s own type.
"""
function PositionVelocityTime.signal_groups(
    track_state::Tracking.TrackState,
    decoders::NamedTuple,
)
    map(track_state.groups, decoders) do group, group_decoders
        PositionVelocityTime.SignalGroup(group, group_decoders)
    end
end

end
