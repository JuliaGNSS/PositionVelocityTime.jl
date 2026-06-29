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
        broadcast_group_delay_e1_e5a = -3.5e-9,
        signal_health_e5a = GNSSDecoder.signal_ok,
        data_validity_status_e5a = GNSSDecoder.navigation_data_valid,
    )
    swap(data) = GNSSDecoder.GNSSDecoderState(e5a_state; data = data, raw_data = data)
    e5a = swap(e5a_eph)

    @testset "classification / signal id / week / start time" begin
        @test PositionVelocityTime.get_time_system(e5a) == :Galileo
        @test PositionVelocityTime.get_signal_id(SatelliteState(;
            decoder = e5a, system = GalileoE5aI(), code_phase = 0.0,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)) == :GalileoE5aI
        @test PositionVelocityTime.get_week(e5a) == d.WN                          # full WN, no rollover
        @test PositionVelocityTime.get_system_start_time(e5a) == PositionVelocityTime.get_system_start_time(base_e1b)
    end

    @testset "orbit matches the equivalent E1B Keplerian" begin
        for t in (132000.0, 132769.0, 135000.0)
            pe = PositionVelocityTime.calc_satellite_position_and_velocity(base_e1b, t)
            pa = PositionVelocityTime.calc_satellite_position_and_velocity(e5a, t)
            @test pa.position ≈ pe.position atol = 1e-6
            @test pa.velocity ≈ pe.velocity atol = 1e-9
        end
    end

    @testset "group delay uses the E1–E5a BGD (per band, signal-agnostic)" begin
        # Galileo BGD is per band, so it depends on the decoder, not the ranging
        # signal — pass the data and pilot components and expect the same result.
        @test PositionVelocityTime.correct_by_group_delay(e5a, GalileoE5aI(), 100.0) ≈
              100.0 - e5a.data.broadcast_group_delay_e1_e5a
        @test PositionVelocityTime.correct_by_group_delay(e5a, GalileoE5aQ(), 100.0) ≈
              100.0 - e5a.data.broadcast_group_delay_e1_e5a
        # The E1B I/NAV path applies the E1–E5b BGD instead — a different field.
        @test PositionVelocityTime.correct_by_group_delay(base_e1b, GalileoE1B(), 100.0) ≈
              100.0 - base_e1b.data.broadcast_group_delay_e1_e5b
    end

    @testset "NTCM-G coefficients (shared model with E1B)" begin
        @test PositionVelocityTime.ntcm_g_params(e5a) === nothing            # none broadcast above
        with_iono = swap(GNSSDecoder.GalileoE5aData(e5a_eph; a_i0 = 45.0, a_i1 = 0.1, a_i2 = 0.01))
        p = PositionVelocityTime.ntcm_g_params(with_iono)
        @test p isa PositionVelocityTime.NTCMGParams
        @test (p.a_i0, p.a_i1, p.a_i2, p.week_number) == (45.0, 0.1, 0.01, d.WN)
    end

    @testset "GGTO (shared word-type-10 layout with E1B)" begin
        @test !PositionVelocityTime.ggto_available(e5a)
        g = swap(GNSSDecoder.GalileoE5aData(e5a_eph; A_0G = 5e-9, A_1G = 1e-15, t_0G = 100, WN_0G = 1134))
        @test PositionVelocityTime.ggto_available(g)
        @test PositionVelocityTime.calc_ggto_offset(g, 132000.0) ≈
              5e-9 + 1e-15 * (132000.0 - 100 + 604800 * mod(d.WN - 1134, 64))
    end
end
