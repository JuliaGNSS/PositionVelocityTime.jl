module PositionVelocityTimeTrackingExt

using PositionVelocityTime: PositionVelocityTime, SatelliteState
using GNSSDecoder: GNSSDecoderState
using GNSSSignals: AbstractGNSSSignal
using Tracking: Tracking, get_code_phase, get_carrier_doppler, get_carrier_phase

# The name of the per-satellite tracking state changed between Tracking major versions:
# `SatState` up to Tracking 1, `TrackedSat` from Tracking 2 on. Bind to whichever exists so
# the extension supports both.
const TrackingSatState = isdefined(Tracking, :TrackedSat) ? Tracking.TrackedSat : Tracking.SatState

function PositionVelocityTime.SatelliteState(
    decoder::GNSSDecoderState,
    system::AbstractGNSSSignal,
    sat_state::TrackingSatState,
)
    SatelliteState(
        decoder,
        system,
        get_code_phase(sat_state),
        get_carrier_doppler(sat_state),
        get_carrier_phase(sat_state),
    )
end

end
