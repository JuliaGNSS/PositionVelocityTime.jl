@testset "PVT Galileo E1B with frequency offset of $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    galileo_e1b = GalileoE1B()
    states = galileo_e1b_states(freq_offset)

    # Fixture data was recorded on 2021-05-31. Pin `approximate_year` so
    # this test stays valid as the wall clock drifts past 2030 (when the
    # default `now()` anchor would start picking the wrong rollover cycle
    # for these archived fixtures).
    pvt = calc_pvt(
        states;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    expected_lla = LLA(; lat = 50.778851672464015, lon = 6.065568885758519, alt = 289.4069805158367)
    @test pvt.position ≈ ECEFfromLLA(wgs84)(expected_lla) rtol = 1e-8
    @test pvt.time ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.1183385390904732)
    @test pvt.velocity ≈ ECEF(0.0, 0.0, 0.0) atol = 9
    @test pvt.relative_clock_drift * get_center_frequency(galileo_e1b) ≈ -(1675.63Hz + freq_offset) atol = 0.01Hz

    warm_pvt = calc_pvt(
        states,
        pvt;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt)
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end
@testset "PVT GPS L1 with frequency offset of $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    gpsl1 = GPSL1CA()
    states = gps_l1_states(freq_offset)

    # Fixture data was recorded on 2021-05-31. See the note on the
    # Galileo testset above for why we pin `approximate_year`.
    pvt = calc_pvt(
        states;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    expected_lla = LLA(; lat = 50.77885249310784, lon = 6.0656199911189175, alt = 291.95658091689086)
    @test pvt.position ≈ ECEFfromLLA(wgs84)(expected_lla) rtol = 1e-8
    @test pvt.time ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.1491024351271335)
    @test pvt.velocity ≈ ECEF(0.0, 0.0, 0.0) atol = 2.5
    @test pvt.relative_clock_drift * get_center_frequency(gpsl1) ≈ -(1632.59Hz + freq_offset) atol = 0.01Hz

    warm_pvt = calc_pvt(
        states,
        pvt;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt)
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end

# Injects synthetic GGTO parameters into a Galileo SatelliteState's decoder.
function with_ggto(state; A_0G, A_1G = 0.0, t_0G = 0, WN_0G = 0)
    d = state.decoder
    new_data = GNSSDecoder.GalileoE1BData(d.data; A_0G, A_1G, t_0G, WN_0G)
    new_raw = GNSSDecoder.GalileoE1BData(d.raw_data; A_0G, A_1G, t_0G, WN_0G)
    new_decoder = GNSSDecoder.GNSSDecoderState(d; raw_data = new_raw, data = new_data)
    SatelliteState(;
        decoder = new_decoder,
        system = state.system,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end

# Relabels a SatelliteState's decoder PRN. The PRN is not used in the solve, only
# as the `sats`-dict key, so this isolates the (system, PRN) keying behaviour.
function with_prn(state, prn)
    d = state.decoder
    new_decoder = GNSSDecoder.GNSSDecoderState(;
        prn = prn,
        raw_data = d.raw_data,
        data = d.data,
        constants = d.constants,
        cache = d.cache,
        num_bits_after_valid_syncro_sequence = d.num_bits_after_valid_syncro_sequence,
        is_shifted_by_180_degrees = d.is_shifted_by_180_degrees,
    )
    SatelliteState(;
        decoder = new_decoder,
        system = state.system,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end

@testset "PVT GPS L1 + Galileo E1B combined, frequency offset $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    states = [gps_l1_states(freq_offset); galileo_e1b_states(freq_offset)]

    # Independent inter-system-bias solve: 8 GPS + 5 Galileo satellites, with one
    # clock bias per system (3 + 2 unknowns).
    pvt = calc_pvt(states; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    expected_pos = ECEF(4.0186793403460276e6, 427032.98199905816, 4.918251309950126e6)
    @test pvt.position ≈ expected_pos rtol = 1e-8
    @test pvt.velocity ≈ ECEF(-1.4405743822415678, 0.5393693783528187, -2.135825176671574) atol = 1e-3
    @test pvt.time ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.285)

    # One clock bias per GNSS; the reference (GPS, most satellites) bias is
    # reported as time_correction, with Galileo's offset relative to it.
    @test pvt.reference_system == GPST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GST()])

    # The combined fix agrees with each single-system fix and with the
    # GPS-anchored time scale.
    gps_only = calc_pvt(gps_l1_states(freq_offset); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gal_only = calc_pvt(galileo_e1b_states(freq_offset); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test norm(pvt.position - gps_only.position) < 5
    @test norm(pvt.position - gal_only.position) < 5
    @test pvt.time ≈ gps_only.time

    # A single receiver clock drift is estimated for the whole constellation
    # (one oscillator), so the combined drift agrees with each single-system drift
    # — to well within the inter-system spread, not split per system.
    @test pvt.relative_clock_drift ≈ gps_only.relative_clock_drift atol = 5e-8
    @test pvt.relative_clock_drift ≈ gal_only.relative_clock_drift atol = 5e-8

    # Satellites are keyed by (signal, PRN), so both constellations are represented
    # and nothing is merged or dropped across systems: the combined count equals the
    # sum of the single-system counts (PRNs here don't overlap).
    @test all(k -> k isa Tuple{Symbol,Int}, keys(pvt.sats))
    @test Set(first(k) for k in keys(pvt.sats)) == Set([:GPSL1CA, :GalileoE1B])
    @test length(pvt.sats) ==
          length(gps_only.sats) + length(gal_only.sats)
    let (sig, prn) = first(keys(pvt.sats))
        @test get_sat_info(pvt, sig, prn) === pvt.sats[(sig, prn)]
    end
    # A satellite not used in the fix returns `nothing` rather than throwing.
    @test get_sat_info(pvt, :GPSL1CA, 999) === nothing

    # Per-satellite post-fit residuals are exposed via SatInfo (modeled − measured
    # pseudorange, metres); finite and small for a good fix.
    resids = [info.residual for info in values(pvt.sats)]
    @test length(resids) == length(pvt.sats)
    @test all(isfinite, resids)
    @test maximum(abs, resids) < 10

    warm_pvt = calc_pvt(states, pvt; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt)
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end

@testset "PVT GPS + Galileo GGTO fallback" begin
    C = 299792458.0
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)

    # The full independent solution provides the reference position and the
    # inter-system time offset that the GGTO must encode.
    reference = calc_pvt([gps; gal]; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gps_only = calc_pvt(gps; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    # reference_system is GPS (most satellites), so the Galileo inter-system bias
    # is Galileo − GPS = −c·(GST − GPST). The physical broadcast GGTO that would
    # reproduce it is its negation, in seconds.
    Δ = -reference.inter_system_biases[GST()] / C

    # 3 GPS + 1 Galileo: the independent solve needs 3 + 2 = 5 satellites, so
    # without GGTO the constellation is under-determined and calc_pvt returns the
    # (origin) previous solution.
    subset = [gps[1:3]; gal[1:1]]
    @test calc_pvt(subset; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false).position == ECEF(0, 0, 0)

    # With the GGTO available the Galileo clock bias collapses onto GPS, so a
    # 4-satellite fix becomes possible and reproduces the full-constellation
    # position. A wrong offset sign would corrupt the single Galileo measurement
    # (here c·Δ ≈ 2.4e6 m), so reproducing the reference also pins the sign.
    subset_ggto = [gps[1:3]; [with_ggto(gal[1]; A_0G = Δ)]]
    pvt = calc_pvt(subset_ggto; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt.position != ECEF(0, 0, 0)
    @test norm(pvt.position - reference.position) < 10
    @test norm(pvt.position - gps_only.position) < 10
    @test pvt.reference_system == GPST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GST()])

    # A single Galileo satellite carrying the GGTO is enough to collapse the
    # whole Galileo set: here only the first of two Galileo satellites has it.
    @test !PositionVelocityTime.ggto_available(gal[2].decoder)
    mixed = [gps[1:2]; [with_ggto(gal[1]; A_0G = Δ), gal[2]]]
    pvt_mixed = calc_pvt(mixed; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt_mixed.position != ECEF(0, 0, 0)
    @test norm(pvt_mixed.position - reference.position) < 50

    # The GGTO offset is applied with the exact sign and magnitude: the Galileo
    # inter-system bias (relative to the GPS reference) equals −c·A_0G, because a
    # Galileo measurement is moved into the GPS frame by subtracting GST − GPST.
    A_big = 1.0e-6
    big = calc_pvt([gps[1:3]; [with_ggto(gal[1]; A_0G = A_big)]]; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test big.inter_system_biases[GST()] ≈ -C * A_big rtol = 1e-6

    # calc_ggto_offset evaluates the OS SIS ICD word-type-10 polynomial, taking
    # the reference week difference modulo 64.
    g = with_ggto(gal[1]; A_0G = 5.0e-9, A_1G = 1.0e-15, t_0G = 100, WN_0G = 1134)
    @test PositionVelocityTime.calc_ggto_offset(g.decoder, 132000.0) ≈
          5.0e-9 + 1.0e-15 * (132000.0 - 100 + 604800 * mod(1136 - 1134, 64))
    @test PositionVelocityTime.ggto_available(g.decoder)
    @test !PositionVelocityTime.ggto_available(gal[1].decoder)
    @test !PositionVelocityTime.ggto_available(gps[1].decoder)
end

@testset "PVT primary system is the most-populated GNSS" begin
    # Galileo-majority independent mix (2 GPS + 5 Galileo): no GGTO collapse, so
    # the primary system — which sets the reported time scale and time_correction
    # — should be Galileo, the constellation with more satellites.
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)
    states = [gps[1:2]; gal]
    pvt = calc_pvt(states; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gal_only = calc_pvt(gal; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)

    @test pvt.reference_system == GST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GPST()])
    @test pvt.time ≈ gal_only.time
    @test norm(pvt.position - gal_only.position) < 10
end

@testset "PVT sats keyed by (signal, PRN): same PRN across signals coexists" begin
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)
    healthy(s) = PositionVelocityTime.is_sat_healthy(s.decoder)

    # Pick a healthy satellite from each constellation and force the Galileo one to
    # share the GPS satellite's PRN, so the only thing telling them apart is the
    # signal tag in the key.
    gps_h = first(s for s in gps if healthy(s))
    gal_h_idx = findfirst(healthy, gal)
    shared_prn = gps_h.decoder.prn
    collided = with_prn(gal[gal_h_idx], shared_prn)
    states = [gps; collided; gal[setdiff(1:length(gal), gal_h_idx)]]

    baseline = calc_pvt([gps; gal]; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    pvt = calc_pvt(states; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)

    # Both the GPS and the Galileo satellite with the shared PRN survive as
    # distinct entries — with a bare-Int PRN key one would have overwritten the
    # other and the used-satellite count would have dropped.
    @test haskey(pvt.sats, (:GPSL1CA, shared_prn))
    @test haskey(pvt.sats, (:GalileoE1B, shared_prn))
    @test get_sat_info(pvt, :GPSL1CA, shared_prn) !==
          get_sat_info(pvt, :GalileoE1B, shared_prn)
    @test length(pvt.sats) == length(baseline.sats)
end
