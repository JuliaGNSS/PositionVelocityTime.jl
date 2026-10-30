# GPS CNAV (L5, GPSCNAVData) and CNAV-2 (L1C-D, GPSCNAV2Data) support.
#
# No real CNAV capture is available, so the orbit is validated by self-consistency:
# the analytic velocity must match a finite difference of the position.
@testset "CNAV (GPS L5 / L1C-D)" begin
    l5 = GNSSDecoder.GPSL5IDecoderState(1)
    l1cd = GNSSDecoder.GPSL1C_DDecoderState(1)
    swap(state, d) = GNSSDecoder.GNSSDecoderState(state; data = d, raw_data = d)

    @testset "classification" begin
        @test GNSSSignals.get_time_system(GNSSSignals.GPSL5I()) == GNSSSignals.GPST()
        @test GNSSSignals.get_time_system(GNSSSignals.GPSL1C_D()) == GNSSSignals.GPST()
        @test !has_offset(l5, GPST())
        @test !has_offset(l1cd, GPST())
    end

    @testset "get_broadcast_week uses full broadcast WN (no rollover)" begin
        # approximate_year must be ignored for CNAV (full 13-bit WN).
        @test GNSSDecoder.get_broadcast_week(swap(l5, GNSSDecoder.GPSCNAVData(; WN = 2200)); approximate_year = 2042) ==
              2200
        @test GNSSDecoder.get_broadcast_week(swap(l1cd, GNSSDecoder.GPSCNAV2Data(; WN = 2200)); approximate_year = 2010) ==
              2200
    end

    @testset "Klobuchar parameters from CNAV α_0…β_3" begin
        for (state, T) in ((l5, GNSSDecoder.GPSCNAVData), (l1cd, GNSSDecoder.GPSCNAV2Data))
            d = T(; α_0 = 1e-8, α_1 = 2e-8, α_2 = 3e-8, α_3 = 4e-8,
                β_0 = 1e4, β_1 = 2e4, β_2 = 3e4, β_3 = 4e4)
            kp = GNSSDecoder.get_klobuchar_parameters(swap(state, d))
            @test kp isa GNSSDecoder.KlobucharParameters
            @test kp.α_0 == 1e-8
            @test kp.β_3 == 4e4
        end
        # Partially populated → nothing, not an error.
        @test GNSSDecoder.get_klobuchar_parameters(swap(l5, GNSSDecoder.GPSCNAVData(; α_0 = 1e-8))) === nothing
    end

    # The CNAV and CNAV-2 elements reach the propagator as GNSSDecoder's
    # `BroadcastEphemeris`, off its reference values (tested there). What is left here
    # is the propagator's own handling of the quasi-Keplerian rates the Keplerian
    # messages lack: the analytic velocity must match a finite difference of the
    # position with `A_dot` and `Δn_0_dot` set.
    @testset "CNAV velocity matches position finite difference (A_dot, Δn_0_dot)" begin
        keplerian = ephemeris(first(gps_l1_states(0.0Hz)).decoder)
        quasi = with_fields(keplerian; A_dot = 0.05, Δn_0_dot = 2e-14)
        position(t) =
            PositionVelocityTime.calc_satellite_position_and_velocity(quasi, t).position
        for t in keplerian.t_0e .+ (300.0, 2100.0)
            pv = PositionVelocityTime.calc_satellite_position_and_velocity(quasi, t)
            h = 0.5
            @test pv.velocity ≈ (position(t + h) .- position(t - h)) ./ (2h) atol = 1e-3
        end
    end

    # The group delays a CNAV/CNAV-2 range takes, pilots and the not-yet-broadcast terms
    # included, are GNSSDecoder's `get_group_delay`, tested there.

    # A pilot ranging signal has 0 Hz data rate; the transmit-time bit-count term must
    # use the decoder's data-signal rate instead, so ranging the band on the L1C pilot
    # gives the same (finite) uncorrected time as ranging on the L1C data component.
    @testset "pilot ranging signal: transmit time well-defined" begin
        d = GNSSDecoder.GPSCNAV2Data(; ITOW = 18, toi = 10)
        dec = GNSSDecoder.GNSSDecoderState(l1cd; data = d, raw_data = d,
            num_bits_since_time_reference = 50)
        common = (; decoder = dec, code_phase = 123.0,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)
        t_data =
            PositionVelocityTime.calc_uncorrected_transmit_time(SatelliteState(; signal = GNSSSignals.GPSL1C_D(), common...))
        t_pilot =
            PositionVelocityTime.calc_uncorrected_transmit_time(SatelliteState(; signal = GNSSSignals.GPSL1C_P(), common...))
        @test isfinite(t_pilot)
        @test t_pilot ≈ t_data
    end
end
