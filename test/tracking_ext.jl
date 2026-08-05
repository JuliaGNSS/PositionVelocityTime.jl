using Tracking
using Tracking:
    TrackedSat, get_code_phase, get_carrier_doppler, get_carrier_phase, estimate_cn0

@testset "Tracking extension: SatelliteState from a Tracking sat state" begin
    # The extension is only available once Tracking is loaded.
    @test Base.get_extension(PositionVelocityTime, :PositionVelocityTimeTrackingExt) !==
          nothing

    gpsl1 = GPSL1CA()
    decoder = GNSSDecoderState(gpsl1, 1)
    sat_state = TrackedSat(gpsl1, 1, 123.0, 1234.0Hz; carrier_phase = 0.25)

    state = SatelliteState(decoder, gpsl1, sat_state)

    @test state isa SatelliteState
    @test state.decoder === decoder
    @test state.system === gpsl1
    # Each measurement must be wired from the matching Tracking getter.
    @test state.code_phase == get_code_phase(sat_state)
    @test state.carrier_doppler == get_carrier_doppler(sat_state)
    @test state.carrier_phase == get_carrier_phase(sat_state)

    # C/N₀ comes along automatically, so a tracking receiver gets the weighted solve
    # without wiring anything up. This satellite's estimator has seen no prompt yet, so
    # Tracking reports `0.0dBHz` — which the variance model reads as "not reported"
    # rather than as an unusably weak signal, leaving this epoch unweighted.
    @test state.cn0 == estimate_cn0(sat_state, typeof(gpsl1))
    @test state.cn0 == 0.0dBHz
    @test !PositionVelocityTime.has_measurement_uncertainty(state)
    @test isnothing(state.pseudorange_variance)

    # Both are overridable for a caller that owns its own error model.
    overridden = SatelliteState(
        decoder,
        gpsl1,
        sat_state;
        cn0 = 42.0dBHz,
        pseudorange_variance = (3.0m)^2,
    )
    @test overridden.cn0 == 42.0dBHz
    @test overridden.pseudorange_variance == (3.0m)^2
    @test PositionVelocityTime.has_measurement_uncertainty(overridden)
    @test isnothing(SatelliteState(decoder, gpsl1, sat_state; cn0 = nothing).cn0)
end
