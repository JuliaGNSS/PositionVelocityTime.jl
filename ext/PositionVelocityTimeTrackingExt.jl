module PositionVelocityTimeTrackingExt

using PositionVelocityTime: PositionVelocityTime, SatelliteState
using GNSSDecoder: GNSSDecoderState
using GNSSSignals: AbstractGNSSSignal
using Tracking:
    Tracking, get_code_phase, get_carrier_doppler, get_carrier_phase, estimate_cn0

function PositionVelocityTime.SatelliteState(
    decoder::GNSSDecoderState,
    system::AbstractGNSSSignal,
    sat_state::Tracking.TrackedSat;
    cn0 = tracked_cn0(sat_state, system),
    pseudorange_variance = nothing,
)
    SatelliteState(;
        decoder,
        system,
        code_phase = get_code_phase(sat_state),
        carrier_doppler = get_carrier_doppler(sat_state),
        carrier_phase = get_carrier_phase(sat_state),
        cn0,
        pseudorange_variance,
    )
end

# C/N₀ of the signal the pseudorange was generated on, so a tracking receiver gets the
# weighted solve without wiring anything up (pass `cn0` explicitly to override, or
# `cn0 = nothing` to opt out). `estimate_cn0` selects a signal of a multi-signal
# satellite by type, which is ambiguous when the satellite is tracked on two signals of
# that type and an error when it is tracked on none; neither is worth failing a fix
# over, so both fall back to the satellite's estimator-driver signal (`signals[1]`).
# A satellite whose estimator has not integrated a prompt yet reports `0.0dBHz`, which
# the variance model treats as "not reported" rather than as an unusably weak signal.
function tracked_cn0(sat_state::Tracking.TrackedSat, system::AbstractGNSSSignal)
    signals = Tracking.get_signals(sat_state)
    matching = count(sig -> Tracking.get_signal(sig) isa typeof(system), signals)
    estimate_cn0(sat_state, matching == 1 ? typeof(system) : 1)
end

end
