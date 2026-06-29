# GPS CNAV (L5, GPSCNAVData) and CNAV-2 (L1C-D, GPSL1C_DData) support.
#
# No real CNAV capture is available, so the orbit is validated by self-consistency
# (Q2): a CNAV ephemeris with zero deltas (ΔA from A_REF, no rates, ΔΩ_dot relative
# to Ω̇_REF) must reproduce the directly-broadcast Keplerian (LNAV) result exactly,
# and the analytic velocity must match a finite difference of the position.
@testset "CNAV (GPS L5 / L1C-D)" begin
    A_REF = 26_559_710.0
    Ωdot_REF = -2.6e-9 * π

    l5 = GNSSDecoder.GPSL5IDecoderState(1)
    l1cd = GNSSDecoder.GPSL1C_DDecoderState(1)
    swap(state, d) = GNSSDecoder.GNSSDecoderState(state; data = d, raw_data = d)

    @testset "classification" begin
        @test PositionVelocityTime.get_time_system(l5) == :GPS
        @test PositionVelocityTime.get_time_system(l1cd) == :GPS
        @test !PositionVelocityTime.ggto_available(l5)
        @test !PositionVelocityTime.ggto_available(l1cd)
    end

    @testset "get_week uses full broadcast WN (no rollover)" begin
        # approximate_year must be ignored for CNAV (full 13-bit WN).
        @test PositionVelocityTime.get_week(swap(l5, GNSSDecoder.GPSCNAVData(; WN = 2200)); approximate_year = 2042) ==
              2200
        @test PositionVelocityTime.get_week(swap(l1cd, GNSSDecoder.GPSL1C_DData(; WN = 2200)); approximate_year = 2010) ==
              2200
    end

    @testset "Klobuchar params from CNAV α_0…β_3" begin
        for (state, T) in ((l5, GNSSDecoder.GPSCNAVData), (l1cd, GNSSDecoder.GPSL1C_DData))
            d = T(; α_0 = 1e-8, α_1 = 2e-8, α_2 = 3e-8, α_3 = 4e-8,
                β_0 = 1e4, β_1 = 2e4, β_2 = 3e4, β_3 = 4e4)
            kp = PositionVelocityTime.klobuchar_params(swap(state, d))
            @test kp isa PositionVelocityTime.KlobucharParams
            @test kp.α_0 == 1e-8
            @test kp.β_3 == 4e4
        end
        # Partially populated → nothing, not an error.
        @test PositionVelocityTime.klobuchar_params(swap(l5, GNSSDecoder.GPSCNAVData(; α_0 = 1e-8))) === nothing
    end

    # Shared Keplerian element set (realistic GPS MEO).
    kep = (; M_0 = 0.3, e = 0.012, sqrt_A = 5153.65, Ω_0 = 1.0, i_0 = 0.96, ω = -0.5,
        Δn = 4.7e-9, Ω_dot = -8.1e-9, i_dot = 1.3e-10, C_uc = 1e-6, C_us = 2e-6,
        C_rc = 200.0, C_rs = 30.0, C_ic = -1e-7, C_is = 9e-8, t_0e = 300.0)
    lnav = GNSSDecoder.GPSL1CAData(; M_0 = kep.M_0, e = kep.e, sqrt_A = kep.sqrt_A, Ω_0 = kep.Ω_0,
        i_0 = kep.i_0, ω = kep.ω, Δn = kep.Δn, Ω_dot = kep.Ω_dot, i_dot = kep.i_dot,
        C_uc = kep.C_uc, C_us = kep.C_us, C_rc = kep.C_rc, C_rs = kep.C_rs,
        C_ic = kep.C_ic, C_is = kep.C_is, t_0e = kep.t_0e)

    @testset "CNAV orbit reduces to the equivalent Keplerian ($T)" for T in
                                                                       (GNSSDecoder.GPSCNAVData, GNSSDecoder.GPSL1C_DData)
        cnav = T(; M_0 = kep.M_0, e = kep.e, ΔA = kep.sqrt_A^2 - A_REF, A_dot = 0.0,
            Δn_0 = kep.Δn, Δn_0_dot = 0.0, Ω_0 = kep.Ω_0, i_0 = kep.i_0, ω = kep.ω,
            ΔΩ_dot = kep.Ω_dot - Ωdot_REF, i_dot = kep.i_dot, C_uc = kep.C_uc,
            C_us = kep.C_us, C_rc = kep.C_rc, C_rs = kep.C_rs, C_ic = kep.C_ic,
            C_is = kep.C_is, t_0e = kep.t_0e)
        dl = swap(GNSSDecoder.GPSL1CADecoderState(1), lnav)
        ref_state = T === GNSSDecoder.GPSCNAVData ? l5 : l1cd
        dc = swap(ref_state, cnav)
        for t in (300.0, 1500.0, 4000.0)
            pl = PositionVelocityTime.calc_satellite_position_and_velocity(dl, t)
            pc = PositionVelocityTime.calc_satellite_position_and_velocity(dc, t)
            @test pc.position ≈ pl.position atol = 1e-6
            @test pc.velocity ≈ pl.velocity atol = 1e-9
        end
    end

    @testset "CNAV velocity matches position finite difference (A_dot, Δn_0_dot)" begin
        cnav = GNSSDecoder.GPSCNAVData(; M_0 = 0.3, e = 0.012, ΔA = 1234.5, A_dot = 0.05,
            Δn_0 = 4.7e-9, Δn_0_dot = 2e-14, Ω_0 = 1.0, i_0 = 0.96, ω = -0.5,
            ΔΩ_dot = -8.1e-9 - Ωdot_REF, i_dot = 1.3e-10, C_uc = 1e-6, C_us = 2e-6,
            C_rc = 200.0, C_rs = 30.0, C_ic = -1e-7, C_is = 9e-8, t_0e = 300.0)
        dc = swap(l5, cnav)
        for t in (600.0, 2400.0)
            pv = PositionVelocityTime.calc_satellite_position_and_velocity(dc, t)
            h = 0.5
            fd =
                (PositionVelocityTime.calc_satellite_position(dc, t + h) .- PositionVelocityTime.calc_satellite_position(dc, t - h)) ./
                (2h)
            @test pv.velocity ≈ fd atol = 1e-3
        end
    end

    # Group delay is selected by the *ranging* signal, with values read from the
    # CNAV-2 message: one L1C-D decoder serves a range generated on the L1C data,
    # the L1C pilot, or C/A. This decouples ranging signal from nav source — a
    # jointly-tracked L1 band may range on the pilot while decoding L1C-D.
    @testset "group delay selected by ranging signal (CNAV-2 ISC_L1CD/L1CP/L1CA)" begin
        d = GNSSDecoder.GPSL1C_DData(; T_GD = 5.0e-9, ISC_L1CD = 1.0e-9, ISC_L1CP = 2.0e-9,
            ISC_L1CA = 3.0e-9)
        dc = swap(l1cd, d)
        @test PositionVelocityTime.correct_by_group_delay(dc, GNSSSignals.GPSL1C_D(), 100.0) ≈
              100.0 - 5.0e-9 + 1.0e-9
        @test PositionVelocityTime.correct_by_group_delay(dc, GNSSSignals.GPSL1C_P(), 100.0) ≈
              100.0 - 5.0e-9 + 2.0e-9
        @test PositionVelocityTime.correct_by_group_delay(dc, GNSSSignals.GPSL1CA(), 100.0) ≈
              100.0 - 5.0e-9 + 3.0e-9
    end

    # A pilot ranging signal has 0 Hz data rate; the transmit-time bit-count term must
    # use the decoder's data-signal rate instead, so ranging the band on the L1C pilot
    # gives the same (finite) uncorrected time as ranging on the L1C data component.
    @testset "pilot ranging signal: transmit time well-defined" begin
        d = GNSSDecoder.GPSL1C_DData(; ITOW = 18, toi = 10)
        dec = GNSSDecoder.GNSSDecoderState(l1cd; data = d, raw_data = d,
            num_bits_after_valid_syncro_sequence = 50)
        common = (; decoder = dec, code_phase = 123.0,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)
        t_data =
            PositionVelocityTime.calc_uncorrected_time(SatelliteState(; system = GNSSSignals.GPSL1C_D(), common...))
        t_pilot =
            PositionVelocityTime.calc_uncorrected_time(SatelliteState(; system = GNSSSignals.GPSL1C_P(), common...))
        @test isfinite(t_pilot)
        @test t_pilot ≈ t_data
    end
end
