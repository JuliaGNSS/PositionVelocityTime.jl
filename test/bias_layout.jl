# Receiver inter-frequency bias (IFB) estimation: when satellites are processed on
# more than one frequency band, calc_pvt estimates one extra unknown per band beyond
# the reference (shared across constellations on that band).
@testset "inter-frequency bias" begin
    C = 299792458.0
    kw = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)

    # Decompose an uncorrected transmit time `target` (s) into a realistic measurement on
    # `system` with time-of-week `tow`: whole data symbols go into the bit count and only
    # the residual (< one data symbol) into the code phase — matching
    # `calc_uncorrected_time`, which reduces the shared code phase modulo one data symbol.
    # Packing the whole sub-TOW interval into `code_phase` (as a naive copy would) is
    # reduced away by that mod, so the copy must spread it across the bit count.
    function bits_and_code_phase(system, tow, target)
        datafreq = Float64(GNSSSignals.get_data_frequency(system) / Hz)
        codefreq = Float64(GNSSSignals.get_code_frequency(system) / Hz)
        elapsed = target - tow
        num_bits = floor(Int, elapsed * datafreq)
        code_phase = (elapsed - num_bits / datafreq) * codefreq
        (num_bits, code_phase)
    end

    # Build an L5 (E5a) measurement for a Galileo E1B satellite, transmit-time-consistent
    # with it (same ephemeris/clock, BGD matched, observables solved to the same epoch).
    # `ifb_shift_s` injects a uniform receiver L5 delay (seconds); `ggto` (A_0G, seconds)
    # makes the copy carry a broadcast Galileo–GPS Time Offset.
    #
    # "BGD matched" needs care: the E1 range applies `BGD(E1,E5b)` as broadcast while
    # the E5a range applies `(f_E1/f_E5a)²·BGD(E1,E5a)` (OS SIS ICD §5.1.5), so copying
    # the E1B field across verbatim would leave the two clock corrections differing by
    # 0.79 of a BGD — a few nanoseconds, which this test would then read out as a
    # spurious sub-metre L5 inter-frequency bias. Pre-dividing by the scaling is what
    # makes the copy carry the *same* satellite clock correction as its E1B original,
    # which is the property the test is about.
    function as_e5a(state; ifb_shift_s = 0.0, ggto = nothing)
        d = state.decoder.data
        ggto_fields = isnothing(ggto) ? (;) : (; A_0G = ggto, A_1G = 0.0, t_0G = 0, WN_0G = d.WN)
        e5a_data = GNSSDecoder.GalileoE5aData(; WN = d.WN, TOW = d.TOW, t_0e = d.t_0e,
            M_0 = d.M_0, e = d.e, sqrt_A = d.sqrt_A, Ω_0 = d.Ω_0, i_0 = d.i_0, ω = d.ω,
            i_dot = d.i_dot, Ω_dot = d.Ω_dot, Δn = d.Δn, C_uc = d.C_uc, C_us = d.C_us,
            C_rc = d.C_rc, C_rs = d.C_rs, C_ic = d.C_ic, C_is = d.C_is, t_0c = d.t_0c,
            a_f0 = d.a_f0, a_f1 = d.a_f1, a_f2 = d.a_f2,
            BGD_E1_E5a = d.BGD_E1_E5b /
                                           PositionVelocityTime.galileo_group_delay_scaling(
                GalileoE5aI(),
            ),
            E5a_SHS = GNSSDecoder.signal_ok,
            E5a_DVS = GNSSDecoder.navigation_data_valid, ggto_fields...)
        target = PositionVelocityTime.calc_uncorrected_time(state) + ifb_shift_s
        num_bits, code_phase = bits_and_code_phase(GalileoE5aI(), e5a_data.TOW, target)
        dec = GNSSDecoder.GNSSDecoderState(GNSSDecoder.GalileoE5aDecoderState(state.decoder.prn);
            data = e5a_data, raw_data = e5a_data,
            num_bits_after_valid_syncro_sequence = num_bits)
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
        target = PositionVelocityTime.calc_uncorrected_time(state) + ifb_shift_s
        num_bits, code_phase = bits_and_code_phase(GPSL2CM(), l2c_data.TOW, target)
        dec = GNSSDecoder.GNSSDecoderState(GNSSDecoder.GPSL2CMDecoderState(state.decoder.prn);
            data = l2c_data, raw_data = l2c_data,
            num_bits_after_valid_syncro_sequence = num_bits)
        SatelliteState(; decoder = dec, system = GPSL2CM(), code_phase = code_phase,
            carrier_doppler = 0.0Hz, carrier_phase = 0.0)
    end

    # GPS L5 (CNAV) copy of a GPS L1 C/A satellite — the L2C copy above on the L5 band:
    # the same CNAV message type, with the L5 group delay (ISC_L5I5 = 0) and L5 health
    # in place of the L2 ones. Only used to build a third band.
    function as_l5i(state)
        l2c = as_l2c(state)
        d = l2c.decoder.data
        l5_data = GNSSDecoder.GPSCNAVData(;
            (f => getfield(d, f) for f in fieldnames(GNSSDecoder.GPSCNAVData))...,
            ISC_L2C = nothing, l2_health = nothing, ISC_L5I5 = 0.0, l5_health = false)
        target = PositionVelocityTime.calc_uncorrected_time(state)
        num_bits, code_phase = bits_and_code_phase(GPSL5I(), l5_data.TOW, target)
        dec = GNSSDecoder.GNSSDecoderState(GNSSDecoder.GPSL5IDecoderState(state.decoder.prn);
            data = l5_data, raw_data = l5_data,
            num_bits_after_valid_syncro_sequence = num_bits)
        SatelliteState(; decoder = dec, system = GPSL5I(), code_phase = code_phase,
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
        ifb, extra, refs, ncomp = bcl([GPST(), GPST(), GST(), GST()], [:L5, :L5, :L1, :L1])
        @test ncomp == 2
        @test isempty(extra)
        @test isempty(refs)
        @test all(==(0), ifb)
        # Connected via a constellation spanning both bands ⇒ one component, one IFB.
        # L1 and L5 are equally populated ⇒ first-seen (L1) is the reference.
        _, extra2, refs2, ncomp2 = bcl([GPST(), GPST(), GST(), GST()], [:L1, :L5, :L1, :L5])
        @test ncomp2 == 1
        @test extra2 == [:L5]
        @test refs2 == [:L1]
        # Single constellation on two bands ⇒ connected by its shared clock ⇒ one IFB.
        _, extra3, refs3, ncomp3 = bcl([GPST(), GPST()], [:L1, :L5])
        @test ncomp3 == 1
        @test length(extra3) == 1
        @test refs3 == [:L1]
        # Disconnected yet one component carries an IFB: GPS on L1+L2 (connected) with
        # Galileo stranded on L5. The L2 bias is anchored to L1; L5 folds into its clock.
        ifb4, extra4, refs4, ncomp4 = bcl(
            [GPST(), GPST(), GPST(), GPST(), GST(), GST()],
            [:L1, :L1, :L2, :L2, :L5, :L5],
        )
        @test ncomp4 == 2
        @test extra4 == [:L2]
        @test refs4 == [:L1]
        @test ifb4 == [0, 0, 1, 1, 0, 0]
    end

    @testset "layout gate counts both measurements and distinct satellites" begin
        # GPS-only (no GGTO collapse possible): 3 position + 1 clock + num_ifb unknowns.
        # Only the classification fields of a measurement row are read here — the
        # broadcast time offsets are touched on the collapse path only, never reached
        # for a GPS-only constellation — so stand-in rows carrying those suffice.
        decide = PositionVelocityTime.decide_bias_layout
        gate(prns, bands) = decide([
            (; time_system = GPST(), band_id = band, prn = prn, time_offsets = ())
            for (prn, band) in zip(prns, bands)
        ])

        # Dual-band GPS ⇒ 1 IFB ⇒ needs 5 measurements; 4 is too few.
        @test gate(1:5, [:L1, :L5, :L1, :L5, :L1]) !== nothing
        @test gate(1:4, [:L1, :L5, :L1, :L5]) === nothing
        # Single-band GPS ⇒ no IFB ⇒ 4 satellites suffice.
        @test gate(1:4, fill(:L1, 4)) !== nothing

        # The measurement count alone is not enough. Only distinct satellites can
        # constrain the 3 position + 1 clock unknowns, because the extra bands of an
        # already-tracked satellite repeat its line of sight: two satellites on three
        # bands (6 ≥ 3 + 1 + 2) and three satellites on two bands (6 ≥ 3 + 1 + 1) both
        # clear the measurement count while remaining unsolvable.
        @test gate([1, 1, 1, 2, 2, 2], [:L1, :L2, :L5, :L1, :L2, :L5]) === nothing
        @test gate([1, 1, 2, 2, 3, 3], [:L1, :L5, :L1, :L5, :L1, :L5]) === nothing
        # Four distinct satellites on two bands clear both conditions, and so does the
        # leanest dual-frequency constellation: four satellites on L1, one of them also
        # tracked on L5 — that repeat is what makes the IFB observable.
        @test gate([1, 1, 2, 2, 3, 3, 4, 4], repeat([:L1, :L5], 4)) !== nothing
        @test gate([1, 2, 3, 4, 1], [:L1, :L1, :L1, :L1, :L5]) !== nothing
    end

    @testset "single-band fix reports no inter-frequency bias" begin
        # The L1-only GPS+Galileo fixtures share one band ⇒ no IFB unknown.
        pvt = calc_pvt(
            (signal_group(gps_l1_states(0.0Hz)), signal_group(galileo_e1b_states(0.0Hz)));
            kw...)
        @test isempty(pvt.inter_frequency_biases)
        @test length(pvt.sats) >= 4
    end

    # End-to-end through calc_pvt: a Galileo E1B (L1) + E5a (L5) two-band fix. The E5a
    # copies reproduce each E1B satellite's transmit time (same ephemeris/clock, BGD
    # matched, observables solved to the same epoch), so without an injected bias the
    # IFB is zero and the fix is unchanged; a uniform L5 delay is then recovered as the
    # inter-frequency bias rather than corrupting the position.
    @testset "calc_pvt estimates the IFB across the L1 and L5 bands" begin
        e1b = galileo_e1b_states(0.0Hz)
        ref = calc_pvt(signal_group(e1b); kw...)                          # L1-only Galileo fix
        @test length(ref.sats) >= 4

        # Consistent L5 copies ⇒ band grouping triggers, IFB ≈ 0, fix unchanged. The
        # copies reproduce each E1B transmit time exactly, so the residuals are the
        # L1-only fixtures' own ~m-level noise (the baseline), not zero.
        base = maximum(abs, [info.residual for info in values(ref.sats)])
        pvt0 = calc_pvt((signal_group(e1b), signal_group(map(as_e5a, e1b))); kw...)
        @test pvt0.reference_system == GST()
        @test haskey(pvt0.inter_frequency_biases, :L5)
        @test pvt0.inter_frequency_biases[:L5].reference == :L1
        @test length(pvt0.sats) == 2 * length(ref.sats)
        @test norm(pvt0.position - ref.position) < 1e-2
        @test abs(pvt0.inter_frequency_biases[:L5].value) < 1e-2m
        @test maximum(abs, [info.residual for info in values(pvt0.sats)]) ≈ base atol = 1e-2m

        # A uniform 12 m receiver L5 delay (signals appear farther ⇒ earlier transmit
        # time) is absorbed by the IFB, leaving the position and the residuals at the
        # baseline — without the IFB unknown the L5 satellites would carry ~12 m
        # residuals instead.
        δ = 12.0
        pvtδ = calc_pvt(
            (signal_group(e1b), signal_group(map(s -> as_e5a(s; ifb_shift_s = -δ / C), e1b)));
            kw...)
        @test pvtδ.inter_frequency_biases[:L5].value ≈ δ * m atol = 0.05m
        @test norm(pvtδ.position - ref.position) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvtδ.sats)]) ≈ base atol = 0.05m
    end

    # GPS two-band fix: L1 C/A + L2C. The L2C copies reproduce each L1 satellite's
    # transmit time exactly, so without an injected bias the L2 inter-frequency bias is
    # zero and the fix is unchanged; a uniform L2 delay is then recovered as the IFB
    # rather than corrupting the position. Exercises the new L2 band through calc_pvt.
    @testset "calc_pvt estimates the IFB across the GPS L1 and L2 bands" begin
        gps = gps_l1_states(0.0Hz)
        ref = calc_pvt(signal_group(gps); kw...)                          # L1-only GPS fix
        @test length(ref.sats) >= 4
        base = maximum(abs, [info.residual for info in values(ref.sats)])

        # Consistent L2C copies ⇒ band grouping triggers, IFB ≈ 0, fix unchanged.
        pvt0 = calc_pvt((signal_group(gps), signal_group(map(as_l2c, gps))); kw...)
        @test pvt0.reference_system == GPST()
        @test haskey(pvt0.inter_frequency_biases, :L2)
        @test pvt0.inter_frequency_biases[:L2].reference == :L1
        @test length(pvt0.sats) == 2 * length(ref.sats)
        @test norm(pvt0.position - ref.position) < 1e-2
        @test abs(pvt0.inter_frequency_biases[:L2].value) < 1e-2m
        @test maximum(abs, [info.residual for info in values(pvt0.sats)]) ≈ base atol = 1e-2m

        # A uniform 12 m receiver L2 delay is absorbed by the IFB, leaving the position
        # and residuals at the baseline.
        δ = 12.0
        pvtδ = calc_pvt(
            (signal_group(gps), signal_group(map(s -> as_l2c(s; ifb_shift_s = -δ / C), gps)));
            kw...)
        @test pvtδ.inter_frequency_biases[:L2].value ≈ δ * m atol = 0.05m
        @test norm(pvtδ.position - ref.position) < 1e-2
        @test maximum(abs, [info.residual for info in values(pvtδ.sats)]) ≈ base atol = 0.05m
    end

    # Regression for the disjoint-band bug: GPS on L1 only + Galileo on L5 only makes a
    # band's IFB column collinear with the stranded constellation's clock. The
    # component-aware layout must not return a degenerate (huge, canceling) solution.
    @testset "disjoint-band coverage is observable, not degenerate" begin
        gps_l1 = gps_l1_states(0.0Hz)                       # GPS on L1 only
        # Connected all-L1 reference fix, for the true position and inter-system offset.
        connected = calc_pvt(
            (signal_group(gps_l1), signal_group(galileo_e1b_states(0.0Hz))); kw...)

        # No GGTO to reconnect the split ⇒ fold: no IFB column. The L5 copies reproduce
        # the Galileo transmit times exactly, so the fix is the *same* as the connected
        # all-L1 reference — a finite, well-conditioned solution. (We compare to the
        # reference rather than bounding the magnitudes: these fixtures' true
        # inter-system offset is ~2.4e6 m, ~8 ms between the recordings. The degenerate
        # solution the bug produced instead had a ~1e7 GDOP and a huge IFB that does not
        # match the reference at all.)
        gal_l5 = map(as_e5a, galileo_e1b_states(0.0Hz))     # Galileo on L5 only
        pvt = calc_pvt((signal_group(gps_l1), signal_group(gal_l5)); kw...)
        @test isempty(pvt.inter_frequency_biases)            # L5 IFB folded into the Galileo clock
        @test length(pvt.sats) == length(gps_l1) + length(gal_l5)
        @test pvt.inter_system_biases[GST()] ≈ connected.inter_system_biases[GST()] rtol = 1e-6
        @test isfinite(pvt.dop.GDOP)
        @test pvt.dop.GDOP ≈ connected.dop.GDOP rtol = 1e-6
        @test norm(pvt.position - connected.position) < 1e-3
        @test all(isfinite, [info.residual for info in values(pvt.sats)])

        # With a (correct) broadcast GGTO, the collapse reconnects the bands: a clean
        # inter-frequency bias is recovered and the inter-system bias comes from the GGTO.
        true_isb = connected.inter_system_biases[GST()]
        gal_l5_ggto =
            map(s -> as_e5a(s; ggto = -ustrip(m, true_isb) / C), galileo_e1b_states(0.0Hz))
        pvt_ggto = calc_pvt((signal_group(gps_l1), signal_group(gal_l5_ggto)); kw...)
        @test pvt_ggto.reference_system == GPST()              # Galileo collapsed onto GPS
        @test haskey(pvt_ggto.inter_frequency_biases, :L5)   # reconnected ⇒ IFB observable
        @test pvt_ggto.inter_frequency_biases[:L5].reference == :L1
        @test pvt_ggto.inter_system_biases[GST()] ≈ true_isb atol = 5m
        @test abs(pvt_ggto.inter_frequency_biases[:L5].value) < 5m
        @test isfinite(pvt_ggto.dop.GDOP) && 0 < pvt_ggto.dop.GDOP < 1e4
        @test norm(pvt_ggto.position - connected.position) < 10
    end

    # Regression, end to end through calc_pvt: two real satellites tracked on three bands
    # clear the measurement count (6 ≥ 3 + 1 clock + 2 IFB) with two lines of sight, so
    # the design matrix is rank deficient and the epoch is unsolvable. calc_pvt must skip
    # it (return `prev_pvt`) rather than throw — while the layout gate still counted
    # measurements only, the epoch reached the solve and its cold-start step, which is
    # undamped and factorises HᵀH directly, threw SingularException out of LsqFit.
    @testset "unsolvable triple-band geometry is skipped, not thrown" begin
        gps = gps_l1_states(0.0Hz)
        two = gps[1:2]
        groups = (signal_group(two), signal_group(map(as_l2c, two)),
            signal_group(map(as_l5i, two)))
        rows = measurement_rows(groups)
        # All six measurements are usable and clear the measurement count, and it is the
        # distinct-satellite condition — two satellites for 3 + 1 unknowns — that rejects
        # the constellation as unsolvable.
        @test all(PositionVelocityTime.is_sat_healthy(state.decoder) for state in two)
        @test map(row -> row.band_id, rows) == [:L1, :L1, :L2, :L2, :L5, :L5]
        @test length(rows) >= 3 + 1 + 2
        @test PositionVelocityTime.decide_bias_layout(rows) === nothing

        ref = calc_pvt(signal_group(gps); kw...)               # L1-only GPS fix, as previous PVT
        @test calc_pvt(groups, ref; kw...) === ref             # warm start: previous fix kept
        cold = calc_pvt(groups; kw...)                         # cold start: nothing to fall back on
        @test isempty(cold.sats)
        @test iszero(cold.position)
    end
end

# Injects synthetic GGTO parameters into a Galileo SatelliteState's decoder.
function with_ggto(state; A_0G, A_1G = 0.0, t_0G = 0, WN_0G = 0)
    d = state.decoder
    new_data = GNSSDecoder.GalileoINAVData(d.data; A_0G, A_1G, t_0G, WN_0G)
    new_raw = GNSSDecoder.GalileoINAVData(d.raw_data; A_0G, A_1G, t_0G, WN_0G)
    new_decoder = GNSSDecoder.GNSSDecoderState(d; raw_data = new_raw, data = new_data)
    SatelliteState(;
        decoder = new_decoder,
        system = state.system,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end

@testset "PVT GPS + Galileo GGTO fallback" begin
    C = 299792458.0
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)

    # The full independent solution provides the reference position and the
    # inter-system time offset that the GGTO must encode.
    reference = calc_pvt((signal_group(gps), signal_group(gal)); approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gps_only = calc_pvt(signal_group(gps); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    # reference_system is GPS (most satellites), so the Galileo inter-system bias
    # is Galileo − GPS = −c·(GST − GPST). The physical broadcast GGTO that would
    # reproduce it is its negation, in seconds.
    Δ = -ustrip(m, reference.inter_system_biases[GST()]) / C

    # 3 GPS + 1 Galileo: the independent solve needs 3 + 2 = 5 satellites, so
    # without GGTO the constellation is under-determined and calc_pvt returns the
    # (origin) previous solution.
    subset = (signal_group(gps[1:3]), signal_group(gal[1:1]))
    @test calc_pvt(subset; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false).position == ECEF(0, 0, 0)

    # With the GGTO available the Galileo clock bias collapses onto GPS, so a
    # 4-satellite fix becomes possible and reproduces the full-constellation
    # position. A wrong offset sign would corrupt the single Galileo measurement
    # (here c·Δ ≈ 2.4e6 m), so reproducing the reference also pins the sign.
    subset_ggto = (signal_group(gps[1:3]), signal_group([with_ggto(gal[1]; A_0G = Δ)]))
    pvt = calc_pvt(subset_ggto; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt.position != ECEF(0, 0, 0)
    @test norm(pvt.position - reference.position) < 10
    @test norm(pvt.position - gps_only.position) < 10
    @test pvt.reference_system == GPST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GST()])

    # A satellite is identified as (time system, PRN) for the layout's
    # distinct-satellite count, not by PRN alone: a Galileo satellite sharing a GPS PRN
    # is still the fourth distinct satellite the collapsed layout needs, so relabelling
    # it leaves the fix untouched (PRN-only identity would count three and bail).
    collided = with_prn(with_ggto(gal[1]; A_0G = Δ), gps[1].decoder.prn)
    pvt_collided = calc_pvt((signal_group(gps[1:3]), signal_group([collided])); approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt_collided.position == pvt.position

    # A single Galileo satellite carrying the GGTO is enough to collapse the
    # whole Galileo set: here only the first of two Galileo satellites has it.
    @test !PositionVelocityTime.time_offset_available(gal[2].decoder, GPST())
    mixed = (signal_group(gps[1:2]), signal_group([with_ggto(gal[1]; A_0G = Δ), gal[2]]))
    pvt_mixed = calc_pvt(mixed; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt_mixed.position != ECEF(0, 0, 0)
    @test norm(pvt_mixed.position - reference.position) < 50

    # The GGTO offset is applied with the exact sign and magnitude: the Galileo
    # inter-system bias (relative to the GPS reference) equals −c·A_0G, because a
    # Galileo measurement is moved into the GPS frame by subtracting GST − GPST.
    A_big = 1.0e-6
    big = calc_pvt((signal_group(gps[1:3]), signal_group([with_ggto(gal[1]; A_0G = A_big)]));
        approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    @test big.inter_system_biases[GST()] ≈ -C * A_big * m rtol = 1e-6

    # A warm start seeds the shared clock column with the hub's own bias. The Galileo
    # row comes last and maps onto that column too, but its previous absolute bias —
    # the hub's plus its −c·A_0G ≈ −300 m inter-system bias — must not overwrite it.
    big_groups = (signal_group(gps[1:3]), signal_group([with_ggto(gal[1]; A_0G = A_big)]))
    workspace = PVTWorkspace()
    calc_pvt!(PVTSolution(), workspace, big_groups, big; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test big.dop.PDOP < PositionVelocityTime.MAX_TRUSTED_WARM_START_PDOP
    @test workspace.prev_ξ[4] == ustrip(m, big.time_correction)
    # The seed only sets the starting point: warm and cold starts reach the same fix.
    warm = calc_pvt(big_groups, big; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test norm(warm.position - big.position) < 1e-3
    @test warm.time ≈ big.time
    @test warm.inter_system_biases[GST()] ≈ big.inter_system_biases[GST()]

    # calc_steering_offset evaluates the OS SIS ICD word-type-10 polynomial, taking
    # the reference week difference modulo 64. `broadcast_time_offset` is the flattened
    # `BroadcastTimeOffset` the collection pass puts on every measurement row.
    g = with_ggto(gal[1]; A_0G = 5.0e-9, A_1G = 1.0e-15, t_0G = 100, WN_0G = 1134)
    ggto = PositionVelocityTime.broadcast_time_offset(g.decoder, GPST())
    @test PositionVelocityTime.calc_steering_offset(ggto, 132000.0) ≈
          5.0e-9 + 1.0e-15 * (132000.0 - 100 + 604800 * mod(1136 - 1134, 64))
    @test ggto.available
    @test PositionVelocityTime.time_offset_available(g.decoder, GPST())
    @test !PositionVelocityTime.time_offset_available(gal[1].decoder, GPST())
    @test !PositionVelocityTime.time_offset_available(gps[1].decoder, GPST())
    # A target no message broadcasts an offset toward is reported unavailable rather
    # than looked up and missed.
    @test !PositionVelocityTime.broadcast_time_offset(g.decoder, BDT()).available
    # An unavailable offset evaluates to NaN rather than a plausible zero — including
    # through the NaN-`t_0` branch `calc_steering_offset` takes for BeiDou D1/D2, which
    # a NaN `t_0` alone would silently turn into an all-zero polynomial.
    no_offset = PositionVelocityTime.broadcast_time_offset(gal[1].decoder, GPST())
    @test !no_offset.available
    @test isnan(PositionVelocityTime.calc_steering_offset(no_offset, 132000.0))

    # calc_hub_range_offsets turns the rows `decide_bias_layout` selected into the
    # per-satellite range corrections: −c·GGTO for the collapsed (Galileo) satellites at
    # their own transmit times, zero for the anchor system's. An independent layout
    # (no hub at all) needs no conversion at all.
    offset_rows = measurement_rows((
        signal_group(gps[1:2]),
        signal_group([
            with_ggto(gal[1]; A_0G = 5.0e-9, A_1G = 1.0e-15, t_0G = 100, WN_0G = 1134),
            gal[2],
        ]),
    ))
    # The measurement form of `time_offset_available` reads the flag the collection pass
    # flattened onto each row, and has to agree with the decoder form it was derived
    # from — row 3 is the Galileo satellite carrying the GGTO, row 4 the one without,
    # rows 1-2 are GPS, whose offset toward its own time system is not a broadcast
    # quantity.
    available(row, target) = PositionVelocityTime.time_offset_available(row, target)
    @test available(offset_rows[3], GPST())
    @test !available(offset_rows[4], GPST())
    @test !available(offset_rows[1], GPST())
    @test [available(row, GPST()) for row in offset_rows] ==
          [PositionVelocityTime.time_offset_available(state.decoder, GPST()) for state in
           [gps[1:2]; [with_ggto(gal[1]; A_0G = 5.0e-9, A_1G = 1.0e-15, t_0G = 100,
               WN_0G = 1134), gal[2]]]]
    # No message broadcasts an offset toward BDT, so every row reports it unavailable
    # rather than indexing something that is not there.
    @test !any(row -> available(row, BDT()), offset_rows)

    hub_rows = [offset_rows[3]]
    offsets = PositionVelocityTime.calc_hub_range_offsets(offset_rows, hub_rows, GPST())
    @test offsets[[1, 2]] == [0.0, 0.0]
    @test offsets[[3, 4]] ≈ [
        -C * PositionVelocityTime.calc_steering_offset(ggto, row.time) for
        row in offset_rows[3:4]
    ]
    # Both Galileo rows are converted by the one row that carries the offset, each at
    # its own transmit time — so the two corrections differ (by the polynomial's drift
    # term) rather than being one shared scalar.
    @test offsets[3] != offsets[4]
    # `inter_system_bias` is that correction for one row: the offset-carrying row's own,
    # or applied at another row's transmit time; zero for the hub's own system.
    @test PositionVelocityTime.inter_system_bias(offset_rows[3], GPST()) == offsets[3]
    @test PositionVelocityTime.inter_system_bias(
        offset_rows[3], GPST(), offset_rows[4].time) == offsets[4]
    @test PositionVelocityTime.inter_system_bias(offset_rows[1], GPST()) == 0.0
    @test isnan(PositionVelocityTime.inter_system_bias(offset_rows[4], GPST()))

    # The primary system: the hub whenever there is one, otherwise the system with the
    # most rows, a tie going to the one that appears first.
    select = PositionVelocityTime.select_primary_system
    @test select(offset_rows, nothing) == GPST()
    @test select(offset_rows[2:4], nothing) == GST()
    @test select(offset_rows[[3, 1, 2, 4]], nothing) == GST()
    @test select(offset_rows, GST()) == GST()
    @test select(offset_rows[1:0], GPST()) == GPST()
    @test_throws ArgumentError select(offset_rows[1:0], nothing)
    @test PositionVelocityTime.calc_hub_range_offsets(
        offset_rows, PositionVelocityTime.SatelliteMeasurement[], nothing) == zeros(4)
end

@testset "PVT primary system is the most-populated GNSS" begin
    # Galileo-majority independent mix (2 GPS + 5 Galileo): no GGTO collapse, so
    # the primary system — which sets the reported time scale and time_correction
    # — should be Galileo, the constellation with more satellites.
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)
    groups = (signal_group(gps[1:2]), signal_group(gal))
    pvt = calc_pvt(groups; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gal_only = calc_pvt(signal_group(gal); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)

    @test pvt.reference_system == GST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GPST()])
    # Independent clocks, one column per system in order of first appearance.
    layout = PositionVelocityTime.decide_bias_layout(measurement_rows(groups))
    @test layout.clock_systems == [GPST(), GST()]
    @test layout.bias_columns.clock_bias_indices ==
          [fill(1, length(gps[1:2])); fill(2, length(gal))]
    @test pvt.time ≈ gal_only.time
    @test norm(pvt.position - gal_only.position) < 10
end
