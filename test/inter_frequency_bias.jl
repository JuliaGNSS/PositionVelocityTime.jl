# Receiver inter-frequency bias (IFB) estimation: when satellites are processed on
# more than one frequency band, calc_pvt estimates one extra unknown per band beyond
# the reference (shared across constellations on that band).
@testset "inter-frequency bias" begin
    C = 299792458.0
    kw = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)

    # Build an L5 (E5a) measurement for a Galileo E1B satellite, transmit-time-consistent
    # with it (same ephemeris/clock, BGD matched, observables solved to the same epoch).
    # `ifb_shift_s` injects a uniform receiver L5 delay (seconds); `ggto` (A_0G, seconds)
    # makes the copy carry a broadcast Galileo–GPS Time Offset.
    function as_e5a(state; ifb_shift_s = 0.0, ggto = nothing)
        d = state.decoder.data
        ggto_fields = isnothing(ggto) ? (;) : (; A_0G = ggto, A_1G = 0.0, t_0G = 0, WN_0G = d.WN)
        e5a_data = GNSSDecoder.GalileoE5aData(; WN = d.WN, TOW = d.TOW, t_0e = d.t_0e,
            M_0 = d.M_0, e = d.e, sqrt_A = d.sqrt_A, Ω_0 = d.Ω_0, i_0 = d.i_0, ω = d.ω,
            i_dot = d.i_dot, Ω_dot = d.Ω_dot, Δn = d.Δn, C_uc = d.C_uc, C_us = d.C_us,
            C_rc = d.C_rc, C_rs = d.C_rs, C_ic = d.C_ic, C_is = d.C_is, t_0c = d.t_0c,
            a_f0 = d.a_f0, a_f1 = d.a_f1, a_f2 = d.a_f2,
            broadcast_group_delay_e1_e5a = d.broadcast_group_delay_e1_e5b,
            signal_health_e5a = GNSSDecoder.signal_ok,
            data_validity_status_e5a = GNSSDecoder.navigation_data_valid, ggto_fields...)
        dec = GNSSDecoder.GNSSDecoderState(GNSSDecoder.GalileoE5aDecoderState(state.decoder.prn);
            data = e5a_data, raw_data = e5a_data, num_bits_after_valid_syncro_sequence = 0)
        target = PositionVelocityTime.calc_uncorrected_time(state) + ifb_shift_s
        codefreq = Float64(GNSSSignals.get_code_frequency(GalileoE5aI()) / Hz)
        code_phase = (target - PositionVelocityTime.get_tow(dec)) * codefreq
        SatelliteState(; decoder = dec, system = GalileoE5aI(), code_phase = code_phase,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)
    end

    # Build a GPS L2C (CNAV) copy of a GPS L1 C/A satellite: the same orbit and clock in
    # the quasi-Keplerian CNAV form (ΔA from A_REF, ΔΩ̇ from Ω̇_REF, matching T_GD and
    # ISC_L2C = 0), transmit-time-consistent with the L1 original. `ifb_shift_s` injects
    # a uniform receiver L2 delay (seconds).
    function as_l2c(state; ifb_shift_s = 0.0)
        d = state.decoder.data
        A_REF = 26_559_710.0
        Ωdot_REF = -2.6e-9 * π
        l2c_data = GNSSDecoder.GPSCNAVData(;
            WN = PositionVelocityTime.get_week(state.decoder; approximate_year = 2021),
            TOW = d.TOW, t_0e = d.t_0e, M_0 = d.M_0, e = d.e, ΔA = d.sqrt_A^2 - A_REF,
            A_dot = 0.0, Δn_0 = d.Δn, Δn_0_dot = 0.0, Ω_0 = d.Ω_0, i_0 = d.i_0, ω = d.ω,
            i_dot = d.i_dot, ΔΩ_dot = d.Ω_dot - Ωdot_REF, C_uc = d.C_uc, C_us = d.C_us,
            C_rc = d.C_rc, C_rs = d.C_rs, C_ic = d.C_ic, C_is = d.C_is, t_0c = d.t_0c,
            a_f0 = d.a_f0, a_f1 = d.a_f1, a_f2 = d.a_f2, T_GD = d.T_GD, ISC_L2C = 0.0,
            l2_health = false)
        dec = GNSSDecoder.GNSSDecoderState(GNSSDecoder.GPSL2CMDecoderState(state.decoder.prn);
            data = l2c_data, raw_data = l2c_data, num_bits_after_valid_syncro_sequence = 0)
        target = PositionVelocityTime.calc_uncorrected_time(state) + ifb_shift_s
        codefreq = Float64(GNSSSignals.get_code_frequency(GPSL2CM()) / Hz)
        code_phase = (target - PositionVelocityTime.get_tow(dec)) * codefreq
        SatelliteState(; decoder = dec, system = GPSL2CM(), code_phase = code_phase,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)
    end

    @testset "user_position recovers per-band IFB and per-system clocks" begin
        user = ECEFfromLLA(wgs84)(LLA(50.1, 8.7, 120.0))
        ecef_from_enu = ECEFfromENU(user, wgs84)
        # 8 satellites at ~GNSS altitude in varied geometry.
        azels = [(0, 80), (60, 40), (120, 30), (180, 55), (240, 25), (300, 60),
            (30, 70), (150, 20)]
        sat(az, el) = Vector(ecef_from_enu(ENU(
            cosd(el) * sind(az) * 2.02e7, cosd(el) * cosd(az) * 2.02e7, sind(el) * 2.02e7)))
        sat_positions = reduce(hcat, [sat(az, el) for (az, el) in azels])

        # Two time systems (1=GPS, 2=Galileo) and two bands (ref + one extra, e.g. L5).
        # 1=GPS, 2=Galileo clock columns; ifb column 1 = extra band, 0 = reference.
        bias_columns = PositionVelocityTime.BiasColumns([1, 1, 1, 1, 2, 2, 2, 1], 2,
            [0, 0, 0, 1, 0, 0, 1, 1], 1)

        # Truth: position + tc_GPS = 30 m, tc_Gal = 80 m, ifb_extra = 17 m.
        ξ_true = [user[1], user[2], user[3], 30.0, 80.0, 17.0]
        ρ = zeros(8)
        PositionVelocityTime.calc_ρ_hat!(ρ, sat_positions, ξ_true, bias_columns)

        ξ, resid = PositionVelocityTime.user_position(sat_positions, ρ, bias_columns)
        @test ξ[1:3] ≈ [user[1], user[2], user[3]] atol = 1e-3
        @test ξ[4] ≈ 30.0 atol = 1e-4      # GPS clock
        @test ξ[5] ≈ 80.0 atol = 1e-4      # Galileo clock
        @test ξ[6] ≈ 17.0 atol = 1e-4      # inter-frequency bias of the extra band
        @test maximum(abs, resid) < 1e-6   # consistent data → zero residual

        # Drop the IFB column and the same data no longer fits: the extra-band
        # satellites carry an unmodeled 17 m offset.
        bias_columns_no_ifb =
            PositionVelocityTime.BiasColumns([1, 1, 1, 1, 2, 2, 2, 1], 2, zeros(Int, 8), 0)
        _, resid_no_ifb = PositionVelocityTime.user_position(sat_positions, ρ, bias_columns_no_ifb)
        @test maximum(abs, resid_no_ifb) > 1.0
    end

    @testset "band layout follows coverage-graph connectivity" begin
        bcl = PositionVelocityTime.band_ifb_layout
        # Disjoint: GPS only on L5, Galileo only on L1 ⇒ two components, no IFB column
        # (each band is the sole band of its component, so its bias folds into a clock).
        ifb, extra, ncomp = bcl([:GPS, :GPS, :Galileo, :Galileo], [:L5, :L5, :L1, :L1])
        @test ncomp == 2
        @test isempty(extra)
        @test all(==(0), ifb)
        # Connected via a constellation spanning both bands ⇒ one component, one IFB.
        _, extra2, ncomp2 = bcl([:GPS, :GPS, :Galileo, :Galileo], [:L1, :L5, :L1, :L5])
        @test ncomp2 == 1
        @test length(extra2) == 1
        # Single constellation on two bands ⇒ connected by its shared clock ⇒ one IFB.
        _, extra3, ncomp3 = bcl([:GPS, :GPS], [:L1, :L5])
        @test ncomp3 == 1
        @test length(extra3) == 1
    end

    @testset "layout count gate accounts for the extra band" begin
        # GPS-only (no GGTO collapse possible), dummy states: 3 position + 1 clock +
        # num_ifb unknowns. The decoder is only touched on the GGTO path, never reached
        # for a GPS-only constellation, so plain integers stand in for the states.
        decide = PositionVelocityTime.decide_bias_layout
        # Dual-band GPS ⇒ 1 IFB ⇒ needs 5 satellites; 4 is too few.
        @test decide(collect(1:5), fill(:GPS, 5), [:L1, :L5, :L1, :L5, :L1], zeros(5)) !==
              nothing
        @test decide(collect(1:4), fill(:GPS, 4), [:L1, :L5, :L1, :L5], zeros(4)) === nothing
        # Single-band GPS ⇒ no IFB ⇒ 4 satellites suffice.
        @test decide(collect(1:4), fill(:GPS, 4), fill(:L1, 4), zeros(4)) !== nothing
    end

    @testset "single-band fix reports no inter-frequency bias" begin
        # The L1-only GPS+Galileo fixtures share one band ⇒ no IFB unknown.
        kw = (; approximate_year = 2021, enable_ionospheric_correction = false,
            enable_tropospheric_correction = false)
        pvt = calc_pvt([gps_l1_states(0.0Hz); galileo_e1b_states(0.0Hz)]; kw...)
        @test isempty(pvt.inter_frequency_biases)
        @test get_num_used_sats(pvt) >= 4
    end

    # End-to-end through calc_pvt: a Galileo E1B (L1) + E5a (L5) two-band fix. The E5a
    # copies reproduce each E1B satellite's transmit time (same ephemeris/clock, BGD
    # matched, observables solved to the same epoch), so without an injected bias the
    # IFB is zero and the fix is unchanged; a uniform L5 delay is then recovered as the
    # inter-frequency bias rather than corrupting the position.
    @testset "calc_pvt estimates the IFB across the L1 and L5 bands" begin
        e1b = galileo_e1b_states(0.0Hz)
        ref = calc_pvt(e1b; kw...)                          # L1-only Galileo fix
        @test get_num_used_sats(ref) >= 4

        # Consistent L5 copies ⇒ band grouping triggers, IFB ≈ 0, fix unchanged. The
        # copies reproduce each E1B transmit time exactly, so the residuals are the
        # L1-only fixtures' own ~m-level noise (the baseline), not zero.
        base = maximum(abs, [info.residual for info in values(ref.sats)])
        pvt0 = calc_pvt([e1b; map(as_e5a, e1b)]; kw...)
        @test pvt0.reference_system == :Galileo
        @test haskey(pvt0.inter_frequency_biases, :L5)
        @test get_num_used_sats(pvt0) == 2 * get_num_used_sats(ref)
        @test norm(pvt0.position - ref.position) < 1e-2
        @test abs(pvt0.inter_frequency_biases[:L5]) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvt0.sats)]) ≈ base atol = 1e-2

        # A uniform 12 m receiver L5 delay (signals appear farther ⇒ earlier transmit
        # time) is absorbed by the IFB, leaving the position and the residuals at the
        # baseline — without the IFB unknown the L5 satellites would carry ~12 m
        # residuals instead.
        δ = 12.0
        pvtδ = calc_pvt([e1b; map(s -> as_e5a(s; ifb_shift_s = -δ / C), e1b)]; kw...)
        @test pvtδ.inter_frequency_biases[:L5] ≈ δ atol = 0.05
        @test norm(pvtδ.position - ref.position) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvtδ.sats)]) ≈ base atol = 0.05
    end

    # GPS two-band fix: L1 C/A + L2C. The L2C copies reproduce each L1 satellite's
    # transmit time exactly, so without an injected bias the L2 inter-frequency bias is
    # zero and the fix is unchanged; a uniform L2 delay is then recovered as the IFB
    # rather than corrupting the position. Exercises the new L2 band through calc_pvt.
    @testset "calc_pvt estimates the IFB across the GPS L1 and L2 bands" begin
        gps = gps_l1_states(0.0Hz)
        ref = calc_pvt(gps; kw...)                          # L1-only GPS fix
        @test get_num_used_sats(ref) >= 4
        base = maximum(abs, [info.residual for info in values(ref.sats)])

        # Consistent L2C copies ⇒ band grouping triggers, IFB ≈ 0, fix unchanged.
        pvt0 = calc_pvt([gps; map(as_l2c, gps)]; kw...)
        @test pvt0.reference_system == :GPS
        @test haskey(pvt0.inter_frequency_biases, :L2)
        @test get_num_used_sats(pvt0) == 2 * get_num_used_sats(ref)
        @test norm(pvt0.position - ref.position) < 1e-2
        @test abs(pvt0.inter_frequency_biases[:L2]) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvt0.sats)]) ≈ base atol = 1e-2

        # A uniform 12 m receiver L2 delay is absorbed by the IFB, leaving the position
        # and residuals at the baseline.
        δ = 12.0
        pvtδ = calc_pvt([gps; map(s -> as_l2c(s; ifb_shift_s = -δ / C), gps)]; kw...)
        @test pvtδ.inter_frequency_biases[:L2] ≈ δ atol = 0.05
        @test norm(pvtδ.position - ref.position) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvtδ.sats)]) ≈ base atol = 0.05
    end

    # Regression for the disjoint-band bug: GPS on L1 only + Galileo on L5 only makes a
    # band's IFB column collinear with the stranded constellation's clock. The
    # component-aware layout must not return a degenerate (huge, canceling) solution.
    @testset "disjoint-band coverage is observable, not degenerate" begin
        gps_l1 = gps_l1_states(0.0Hz)                       # GPS on L1 only
        # Connected all-L1 reference fix, for the true position and inter-system offset.
        connected = calc_pvt([gps_l1; galileo_e1b_states(0.0Hz)]; kw...)

        # No GGTO to reconnect the split ⇒ fold: no IFB column. The L5 copies reproduce
        # the Galileo transmit times exactly, so the fix is the *same* as the connected
        # all-L1 reference — a finite, well-conditioned solution. (We compare to the
        # reference rather than bounding the magnitudes: these fixtures' true
        # inter-system offset is ~2.4e6 m, ~8 ms between the recordings. The degenerate
        # solution the bug produced instead had a ~1e7 GDOP and a huge IFB that does not
        # match the reference at all.)
        gal_l5 = map(as_e5a, galileo_e1b_states(0.0Hz))     # Galileo on L5 only
        pvt = calc_pvt([gps_l1; gal_l5]; kw...)
        @test isempty(pvt.inter_frequency_biases)            # L5 IFB folded into the Galileo clock
        @test get_num_used_sats(pvt) == length(gps_l1) + length(gal_l5)
        @test pvt.inter_system_biases[:Galileo] ≈ connected.inter_system_biases[:Galileo] rtol = 1e-6
        @test isfinite(pvt.dop.GDOP)
        @test pvt.dop.GDOP ≈ connected.dop.GDOP rtol = 1e-6
        @test norm(pvt.position - connected.position) < 1e-3
        @test all(isfinite, [info.residual for info in values(pvt.sats)])

        # With a (correct) broadcast GGTO, the collapse reconnects the bands: a clean
        # inter-frequency bias is recovered and the inter-system bias comes from the GGTO.
        true_isb = connected.inter_system_biases[:Galileo]
        gal_l5_ggto = map(s -> as_e5a(s; ggto = -true_isb / C), galileo_e1b_states(0.0Hz))
        pvt_ggto = calc_pvt([gps_l1; gal_l5_ggto]; kw...)
        @test pvt_ggto.reference_system == :GPS              # Galileo collapsed onto GPS
        @test haskey(pvt_ggto.inter_frequency_biases, :L5)   # reconnected ⇒ IFB observable
        @test pvt_ggto.inter_system_biases[:Galileo] ≈ true_isb atol = 5
        @test abs(pvt_ggto.inter_frequency_biases[:L5]) < 5
        @test isfinite(pvt_ggto.dop.GDOP) && 0 < pvt_ggto.dop.GDOP < 1e4
        @test norm(pvt_ggto.position - connected.position) < 10
    end
end
