module PositionVelocityTimeTrackingExt

using PositionVelocityTime: PositionVelocityTime, SatelliteState
using GNSSDecoder: GNSSDecoderState
using GNSSSignals: AbstractGNSSSignal
using Tracking: Tracking, get_code_phase, get_carrier_doppler, get_carrier_phase

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

end
