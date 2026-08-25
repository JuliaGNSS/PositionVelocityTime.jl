# Galileo E5a F/NAV (GNSSDecoder.GalileoE5aData) support.
#
# E5a is a separate frequency band from E1, so a satellite tracked on both yields
# two independent measurements. E5a F/NAV shares the Galileo Keplerian ephemeris,
# clock, NTCM-G ionosphere and GGTO layout with the E1B I/NAV message — only its
# broadcast group delay (E1–E5a) and carrier band (1176.45 MHz) differ.
@testset "Galileo E5a (F/NAV) support" begin

    # Build an E5a decoder carrying the same Keplerian ephemeris as an E1B fixture
    # satellite; the Galileo constants (μ, F, Ω̇ₑ) are identical, so the shared orbit
    # path must reproduce the E1B result bit-for-bit.
    base_e1b = galileo_e1b_states(0.0Hz)[1].decoder
    d = base_e1b.data
    e5a_state = GNSSDecoder.GalileoE5aDecoderState(2)
    e5a_eph = GNSSDecoder.GalileoE5aData(;
        WN = d.WN, TOW = d.TOW, t_0e = d.t_0e, M_0 = d.M_0, e = d.e, sqrt_A = d.sqrt_A,
        Ω_0 = d.Ω_0, i_0 = d.i_0, ω = d.ω, i_dot = d.i_dot, Ω_dot = d.Ω_dot, Δn = d.Δn,
        C_uc = d.C_uc, C_us = d.C_us, C_rc = d.C_rc, C_rs = d.C_rs, C_ic = d.C_ic,
        C_is = d.C_is, t_0c = d.t_0c, a_f0 = d.a_f0, a_f1 = d.a_f1, a_f2 = d.a_f2,
        BGD_E1_E5a = -3.5e-9,
        E5a_SHS = GNSSDecoder.signal_ok,
        E5a_DVS = GNSSDecoder.navigation_data_valid,
    )
    swap(data) = GNSSDecoder.GNSSDecoderState(e5a_state; data = data, raw_data = data)
    e5a = swap(e5a_eph)

    @testset "classification / signal id / week / start time" begin
        @test GNSSSignals.get_time_system(GNSSSignals.GalileoE5aI()) == GNSSSignals.GST()
        @test GNSSSignals.get_signal_id(GalileoE5aI()) == :GalileoE5aI
        @test PositionVelocityTime.get_week(e5a) == d.WN                          # full WN, no rollover
        @test PositionVelocityTime.system_start_epoch(GNSSSignals.GalileoE5aI()) ==
              PositionVelocityTime.system_start_epoch(GNSSSignals.GalileoE1B())
    end

    @testset "orbit matches the equivalent E1B Keplerian" begin
        for t in (132000.0, 132769.0, 135000.0)
            pe = PositionVelocityTime.calc_satellite_position_and_velocity(base_e1b, t)
            pa = PositionVelocityTime.calc_satellite_position_and_velocity(e5a, t)
            @test pa.position ≈ pe.position atol = 1e-6
            @test pa.velocity ≈ pe.velocity atol = 1e-9
        end
    end

    @testset "group delay uses the E1–E5a BGD, scaled to the E5a carrier" begin
        # Which BGD applies is per band, so within a band it does not depend on the
        # data/pilot split — E5a-I and E5a-Q get the same correction.
        #
        # How MUCH of it does depend on the band. The broadcast clock is referred to
        # the E1/E5a ionosphere-free combination, so a single-frequency user on the
        # combination's *second* band applies `(f_E1/f_E5a)²·BGD`, not the bare BGD
        # (OS SIS ICD Issue 2.2 §5.1.5). That factor is 1.79, so the difference is
        # ~0.8 of a BGD — a few ns, i.e. metre-level and per-satellite.
        e5a_scaling =
            (get_center_frequency(GalileoE1B) / get_center_frequency(GalileoE5aI()))^2
        @test e5a_scaling ≈ (1575.42 / 1176.45)^2
        @test PositionVelocityTime.correct_by_group_delay(e5a, GalileoE5aI(), 100.0) ≈
              100.0 - e5a_scaling * e5a.data.BGD_E1_E5a
        @test PositionVelocityTime.correct_by_group_delay(e5a, GalileoE5aQ(), 100.0) ≈
              100.0 - e5a_scaling * e5a.data.BGD_E1_E5a
        # F/NAV is decoded on E5a-I only, so E5a is the one band it can be asked about:
        # an E1 range is generated on E1-B and carries an I/NAV decoder. The
        # cross-band pairing is an error rather than a wrong number.
        @test_throws MethodError PositionVelocityTime.correct_by_group_delay(
            e5a, GalileoE1B(), 100.0)
        # The E1B I/NAV path applies the E1–E5b BGD instead — a different field.
        @test PositionVelocityTime.correct_by_group_delay(base_e1b, GalileoE1B(), 100.0) ≈
              100.0 - base_e1b.data.BGD_E1_E5b
    end

    @testset "NTCM-G coefficients (shared model with E1B)" begin
        @test PositionVelocityTime.ntcm_g_params(e5a) === nothing            # none broadcast above
        with_iono = swap(GNSSDecoder.GalileoE5aData(e5a_eph; a_i0 = 45.0, a_i1 = 0.1, a_i2 = 0.01))
        p = PositionVelocityTime.ntcm_g_params(with_iono)
        @test p isa PositionVelocityTime.NTCMGParams
        @test (p.a_i0, p.a_i1, p.a_i2, p.week_number) == (45.0, 0.1, 0.01, d.WN)
    end

    @testset "GGTO (shared word-type-10 layout with E1B)" begin
        @test !PositionVelocityTime.gpst_offset_available(e5a)
        g = swap(GNSSDecoder.GalileoE5aData(e5a_eph; A_0G = 5e-9, A_1G = 1e-15, t_0G = 100, WN_0G = 1134))
        @test PositionVelocityTime.gpst_offset_available(g)
        @test PositionVelocityTime.calc_gpst_offset(g, 132000.0) ≈
              5e-9 + 1e-15 * (132000.0 - 100 + 604800 * mod(d.WN - 1134, 64))
    end
end
