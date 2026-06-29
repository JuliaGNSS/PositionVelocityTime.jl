# GPS L2C (CNAV on the L2 band) support.
#
# L2C shares the CNAV `GNSSDecoder.GPSCNAVData` container with L5-I; the two are
# distinguished only by the decoder's constants type (`GPSL2CMConstants` vs
# `GPSL5IConstants`). So the orbit/clock/week/ionosphere paths are the L5 ones,
# exercised here on an L2C state to confirm the classification, the L2 frequency
# band, the 50 Hz L2 CM symbol rate and the ISC_L2C group delay are wired up.
@testset "GPS L2C (CNAV on L2)" begin
    A_REF = 26_559_710.0
    Ωdot_REF = -2.6e-9 * π

    l2c = GNSSDecoder.GPSL2CMDecoderState(1)
    l5 = GNSSDecoder.GPSL5IDecoderState(1)
    swap(state, d) = GNSSDecoder.GNSSDecoderState(state; data = d, raw_data = d)
    state_of(d) = SatelliteState(; decoder = swap(l2c, d), system = GNSSSignals.GPSL2CM(),
        code_phase = 0.0, carrier_doppler = 0.0Hz, carrier_phase = 0.0)

    @testset "classification: GPS time system, L2 band, GPSL2CM signal id" begin
        st = state_of(GNSSDecoder.GPSCNAVData())
        @test PositionVelocityTime.get_time_system(l2c) == :GPS
        @test !PositionVelocityTime.ggto_available(l2c)
        @test PositionVelocityTime.get_frequency_band(st) == :L2
        @test PositionVelocityTime.get_signal_id(st) == :GPSL2CM
    end

    @testset "nav-message symbol rate is L2 CM (50 Hz), not L5-I (100 Hz)" begin
        d = GNSSDecoder.GPSCNAVData()
        @test PositionVelocityTime.nav_data_frequency(swap(l2c, d)) ==
              GNSSSignals.get_data_frequency(GNSSSignals.GPSL2CM())
        @test PositionVelocityTime.nav_data_frequency(swap(l5, d)) ==
              GNSSSignals.get_data_frequency(GNSSSignals.GPSL5I())
    end

    @testset "get_week uses full broadcast WN (no rollover)" begin
        # approximate_year must be ignored for CNAV (full 13-bit WN).
        @test PositionVelocityTime.get_week(swap(l2c, GNSSDecoder.GPSCNAVData(; WN = 2200));
            approximate_year = 2042) == 2200
    end

    @testset "Klobuchar params from CNAV α_0…β_3" begin
        d = GNSSDecoder.GPSCNAVData(; α_0 = 1e-8, α_1 = 2e-8, α_2 = 3e-8, α_3 = 4e-8,
            β_0 = 1e4, β_1 = 2e4, β_2 = 3e4, β_3 = 4e4)
        kp = PositionVelocityTime.klobuchar_params(swap(l2c, d))
        @test kp isa PositionVelocityTime.KlobucharParams
        @test kp.α_0 == 1e-8
        @test kp.β_3 == 4e4
        # Partially populated → nothing, not an error.
        @test PositionVelocityTime.klobuchar_params(swap(l2c, GNSSDecoder.GPSCNAVData(; α_0 = 1e-8))) ===
              nothing
    end

    # One CNAV decoder carries the whole ISC set (message type 30), so the group delay
    # is selected by the ranging signal: ISC_L2C for an L2C-M range, ISC_L5I5 for an L5
    # range (IS-GPS-200N §30.3.3.3.1.1). -T_GD + ISC.
    @testset "group delay selected by ranging signal (ISC_L2C vs ISC_L5I5)" begin
        dc = swap(l2c, GNSSDecoder.GPSCNAVData(; T_GD = 5.0e-9, ISC_L2C = 7.0e-9, ISC_L5I5 = 1.0e-9))
        @test PositionVelocityTime.correct_by_group_delay(dc, GNSSSignals.GPSL2CM(), 100.0) ≈
              100.0 - 5.0e-9 + 7.0e-9
        @test PositionVelocityTime.correct_by_group_delay(dc, GNSSSignals.GPSL5I(), 100.0) ≈
              100.0 - 5.0e-9 + 1.0e-9
    end

    # The orbit uses the shared CNAV quasi-Keplerian algorithm, so an L2C ephemeris with
    # zero deltas must reproduce the directly-broadcast Keplerian (LNAV) orbit exactly.
    @testset "CNAV orbit on L2C reduces to the equivalent Keplerian" begin
        kep = (; M_0 = 0.3, e = 0.012, sqrt_A = 5153.65, Ω_0 = 1.0, i_0 = 0.96, ω = -0.5,
            Δn = 4.7e-9, Ω_dot = -8.1e-9, i_dot = 1.3e-10, C_uc = 1e-6, C_us = 2e-6,
            C_rc = 200.0, C_rs = 30.0, C_ic = -1e-7, C_is = 9e-8, t_0e = 300.0)
        lnav = GNSSDecoder.GPSL1CAData(; M_0 = kep.M_0, e = kep.e, sqrt_A = kep.sqrt_A,
            Ω_0 = kep.Ω_0, i_0 = kep.i_0, ω = kep.ω, Δn = kep.Δn, Ω_dot = kep.Ω_dot,
            i_dot = kep.i_dot, C_uc = kep.C_uc, C_us = kep.C_us, C_rc = kep.C_rc,
            C_rs = kep.C_rs, C_ic = kep.C_ic, C_is = kep.C_is, t_0e = kep.t_0e)
        cnav = GNSSDecoder.GPSCNAVData(; M_0 = kep.M_0, e = kep.e, ΔA = kep.sqrt_A^2 - A_REF,
            A_dot = 0.0, Δn_0 = kep.Δn, Δn_0_dot = 0.0, Ω_0 = kep.Ω_0, i_0 = kep.i_0,
            ω = kep.ω, ΔΩ_dot = kep.Ω_dot - Ωdot_REF, i_dot = kep.i_dot, C_uc = kep.C_uc,
            C_us = kep.C_us, C_rc = kep.C_rc, C_rs = kep.C_rs, C_ic = kep.C_ic,
            C_is = kep.C_is, t_0e = kep.t_0e)
        dl = swap(GNSSDecoder.GPSL1CADecoderState(1), lnav)
        dc = swap(l2c, cnav)
        for t in (300.0, 1500.0, 4000.0)
            pl = PositionVelocityTime.calc_satellite_position_and_velocity(dl, t)
            pc = PositionVelocityTime.calc_satellite_position_and_velocity(dc, t)
            @test pc.position ≈ pl.position atol = 1e-6
            @test pc.velocity ≈ pl.velocity atol = 1e-9
        end
    end
end
