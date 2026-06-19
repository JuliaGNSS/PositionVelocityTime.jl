using Tracking
using Tracking: TrackedSat, get_code_phase, get_carrier_doppler, get_carrier_phase

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
end
