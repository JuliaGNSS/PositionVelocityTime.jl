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
    expected_lla = LLA(; lat = 50.77885207231635, lon = 6.065568145321566, alt = 289.09688064471146)
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
    expected_lla = LLA(; lat = 50.778851781017025, lon = 6.065622611231713, alt = 291.96260731963366)
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

# Marks a GPS L1 C/A SatelliteState's decoder unhealthy: the six health bits are the
# subframe-1 SV health word, and a set MSB (nav data bad) is what `is_sat_healthy`
# rejects. GNSSDecoder 3 renamed the field `svhealth` -> `sv_health`, as in
# `fixtures.jl`'s `gps_l1ca_data`.
function with_health_bad(state)
    d = state.decoder
    unhealthy(data) =
        :sv_health in fieldnames(GNSSDecoder.GPSL1CAData) ?
        GNSSDecoder.GPSL1CAData(data; sv_health = "100000") :
        GNSSDecoder.GPSL1CAData(data; svhealth = "100000")
    new_decoder = GNSSDecoder.GNSSDecoderState(
        d; raw_data = unhealthy(d.raw_data), data = unhealthy(d.data))
    SatelliteState(;
        decoder = new_decoder,
        system = state.system,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end

# Relabels a SatelliteState's decoder PRN. Nothing in the least-squares solve depends on
# the PRN — it identifies a satellite, as the `sats`-dict key and for the layout's
# distinct-satellite count — so this isolates that identity from the measurements.
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
    expected_pos = ECEF(4.0186793226897363e6, 427033.09443239716, 4.918251247796992e6)
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
    @test maximum(abs, resids) < 10m

    # The rate-domain counterpart from the Doppler solve (modeled − measured range
    # rate, m/s): finite, and metre-per-second-level for a stationary receiver.
    rate_resids = [info.rate_residual for info in values(pvt.sats)]
    @test length(rate_resids) == length(pvt.sats)
    @test all(isfinite, rate_resids)
    @test maximum(abs, rate_resids) < 10.0m/s

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
    Δ = -ustrip(m, reference.inter_system_biases[GST()]) / C

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

    # A satellite is identified as (time system, PRN) for the layout's
    # distinct-satellite count, not by PRN alone: a Galileo satellite sharing a GPS PRN
    # is still the fourth distinct satellite the collapsed layout needs, so relabelling
    # it leaves the fix untouched (PRN-only identity would count three and bail).
    collided = with_prn(with_ggto(gal[1]; A_0G = Δ), gps[1].decoder.prn)
    pvt_collided = calc_pvt([gps[1:3]; [collided]]; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test pvt_collided.position == pvt.position

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
    @test big.inter_system_biases[GST()] ≈ -C * A_big * m rtol = 1e-6

    # calc_ggto_offset evaluates the OS SIS ICD word-type-10 polynomial, taking
    # the reference week difference modulo 64.
    g = with_ggto(gal[1]; A_0G = 5.0e-9, A_1G = 1.0e-15, t_0G = 100, WN_0G = 1134)
    @test PositionVelocityTime.calc_ggto_offset(g.decoder, 132000.0) ≈
          5.0e-9 + 1.0e-15 * (132000.0 - 100 + 604800 * mod(1136 - 1134, 64))
    @test PositionVelocityTime.ggto_available(g.decoder)
    @test !PositionVelocityTime.ggto_available(gal[1].decoder)
    @test !PositionVelocityTime.ggto_available(gps[1].decoder)

    # calc_ggto_range_offsets turns the decoder `decide_bias_layout` selected into the
    # per-satellite range corrections: −c·GGTO for the collapsed (Galileo) satellites at
    # their own transmit times, zero for the anchor system's. An independent layout
    # (`nothing`) needs no conversion at all.
    offset_systems = [GPST(), GST(), GPST(), GST()]
    offset_times = [100.0, 200.0, 300.0, 400.0]
    offsets = PositionVelocityTime.calc_ggto_range_offsets(
        g.decoder, offset_systems, offset_times)
    @test offsets[[1, 3]] == [0.0, 0.0]
    @test offsets[[2, 4]] ≈
          [-C * PositionVelocityTime.calc_ggto_offset(g.decoder, t) for t in (200.0, 400.0)]
    @test PositionVelocityTime.calc_ggto_range_offsets(
        nothing, offset_systems, offset_times) == zeros(4)
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

# Course over ground: azimuth of the velocity vector in the local ENU frame, degrees
# clockwise from North. The recorded fixtures are stationary, so exercise the
# computation directly with known ENU velocities rotated into ECEF at a fixed position.
@testset "calc_course_over_ground is the ENU velocity azimuth clockwise from North" begin
    user = ECEF(ECEFfromLLA(wgs84)(LLA(50.1, 8.7, 120.0)))
    ecef_from_enu = ECEFfromENU(user, wgs84)
    # ENU velocity → ECEF velocity vector: rotate the tip and drop the origin translation.
    ecef_velocity(e, n, u) = ECEF(ecef_from_enu(ENU(e, n, u))) - user
    cog(e, n, u) = PositionVelocityTime.calc_course_over_ground(user, ecef_velocity(e, n, u))

    @test cog(0.0, 10.0, 0.0) ≈ 0° atol = 1e-6°      # North
    @test cog(10.0, 0.0, 0.0) ≈ 90° atol = 1e-6°     # East
    @test cog(0.0, -10.0, 0.0) ≈ 180° atol = 1e-6°   # South
    @test cog(-10.0, 0.0, 0.0) ≈ 270° atol = 1e-6°   # West
    @test cog(10.0, 10.0, 0.0) ≈ 45° atol = 1e-6°    # North-East
    # A stationary receiver has no horizontal velocity ⇒ course defined as 0°.
    @test PositionVelocityTime.calc_course_over_ground(user, ECEF(0.0, 0.0, 0.0)) == 0.0°
    # The vertical component does not affect the horizontal course.
    @test cog(10.0, 0.0, 7.0) ≈ cog(10.0, 0.0, 0.0)
    # Always wrapped to [0, 360)°.
    @test 0° ≤ cog(-3.0, -4.0, 0.0) < 360°
end

# Every unsolvable epoch takes the same exit — `prev_pvt`, never a throw — whether the
# shortfall is in the raw satellite count, in how many of them are usable, or in the layout
# the constellation demands.
@testset "an unsolvable epoch returns prev_pvt instead of throwing" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    gps = gps_l1_states(0.0Hz)
    reference = calc_pvt(gps; kwargs...)
    @test reference.position != ECEF(0, 0, 0)

    # Too few states to hold 3 position + 1 clock unknown, down to none at all: with no
    # previous solution the default (origin) one comes back unchanged.
    for n in 0:3
        pvt = calc_pvt(gps[1:n]; kwargs...)
        @test pvt.position == ECEF(0, 0, 0)
        @test isnothing(pvt.time)
        @test isempty(pvt.sats)
        # Given a previous fix, that fix is what is carried forward.
        @test calc_pvt(gps[1:n], reference; kwargs...) === reference
    end

    # An unhealthy satellite is dropped before the count is taken, so 4 states carrying
    # only 3 usable ones exit exactly as 3 states do.
    unhealthy = [gps[1:3]; [with_health_bad(gps[4])]]
    @test length(unhealthy) == 4
    @test calc_pvt(unhealthy, reference; kwargs...) === reference
end

# Regression: a degenerate geometry must be reported, not thrown out of a least-squares
# solve as a SingularException. Observed while processing a dual-band recording, where an
# epoch's geometry left the velocity normal equations singular and the 4×4 solve threw
# before the epoch could be rejected.
@testset "degenerate geometry is rejected without throwing" begin
    @testset "positive_definite_cholesky rejects what is not positive definite" begin
        pdc = PositionVelocityTime.positive_definite_cholesky
        @test !isnothing(pdc(Symmetric([2.0 0.0; 0.0 3.0])))     # positive definite
        @test isnothing(pdc(Symmetric([1.0 1.0; 1.0 1.0])))      # semidefinite: zero pivot
        @test isnothing(pdc(Symmetric([1.0 0.0; 0.0 -1.0])))     # indefinite
        # Ill-conditioned to the point of being rank deficient in Float64: the pivot is
        # positive, so `issuccess` alone would accept it and return rounding noise.
        @test isnothing(pdc(Symmetric([1.0 0.0; 0.0 1e-14])))
        # A poor but genuinely solvable geometry is still accepted (cond(HᵀH) = 1e6).
        @test !isnothing(pdc(Symmetric([1.0 0.0; 0.0 1e-6])))
    end

    @testset "the geometry checks follow the design matrix rank" begin
        # `calc_DOP` and the velocity solve both decide rank through the normal-equations
        # matrix, so that is what these cases exercise.
        full_rank(H) = !isnothing(
            PositionVelocityTime.positive_definite_cholesky(Symmetric(H' * H)))
        @test full_rank([1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 1.0 0.0])
        @test !full_rank(repeat([1.0 0.0 1.0], 3))               # identical rows
        # A single-band, single-system position design over three satellites: four
        # columns constrained by three lines of sight is always rank deficient — the shape
        # the layout's distinct-satellite condition rejects before the solve.
        @test !full_rank([1.0 0.0 0.0 1.0; 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0])

        # The case no satellite count can see: four *distinct* satellites — plenty for the
        # 3 + 1 unknowns — whose lines of sight lie on a common cone (equal projection onto
        # z), which puts the clock column in their span.
        cosθ = 0.5
        sinθ = sqrt(1 - cosθ^2)
        cone = [sinθ 0.0 cosθ; 0.0 sinθ cosθ; -sinθ 0.0 cosθ; 0.0 -sinθ cosθ]
        @test !full_rank([cone ones(4)])
        # Tilting one satellite off the cone restores full rank.
        tilted = [cone[1:3, :]; normalize([0.0, -sinθ, 2cosθ])']
        @test full_rank([tilted ones(4)])
    end

    # The velocity solve assumes a full-rank position design (its caller establishes that
    # with `calc_DOP` first), so what is pinned here is the implication it relies on: a
    # degenerate line-of-sight set is degenerate for the velocity design too, and would
    # reach a singular solve if the DOP check were ever moved after it.
    @testset "a degenerate line-of-sight set is degenerate for the velocity design" begin
        states = gps_l1_states(0.0Hz)[1:4]
        times = map(PositionVelocityTime.calc_corrected_time, states)
        sat_pvs = map(
            (state, time) ->
                PositionVelocityTime.calc_satellite_position_and_velocity(state.decoder, time),
            states,
            times,
        )
        # Only H's first three columns — the line-of-sight unit vectors — are read; the
        # per-system clock columns are collapsed into a single drift column internally,
        # so one trailing column of ones stands in for them here.
        design(dirs) = reduce(vcat, [[normalize(d)' 1.0] for d in dirs])
        velocity_and_drift(dirs) = PositionVelocityTime.calc_user_velocity_and_clock_drift(
            sat_pvs, states, times, design(dirs))

        # Four well-spread directions ⇒ solvable: a finite [vx, vy, vz, ċ], plus one
        # finite post-fit range-rate residual per satellite.
        solution, rate_residuals = velocity_and_drift(
            [[1.0, 0.2, 0.9], [-0.5, 1.0, 0.7], [0.3, -1.0, 0.5], [0.0, 0.1, 1.0]])
        @test length(solution) == 4
        @test all(isfinite, solution)
        @test length(rate_residuals) == length(states)
        @test all(isfinite, rate_residuals)
        # A post-fit least-squares residual is orthogonal to every column of its design
        # — the line-of-sight columns and the common clock-drift column of ones. That is
        # what makes it a residual rather than an arbitrary modeled-minus-measured
        # difference, so it is pinned directly instead of only by magnitude.
        let A = design([[1.0, 0.2, 0.9], [-0.5, 1.0, 0.7], [0.3, -1.0, 0.5], [0.0, 0.1, 1.0]])
            @test norm(A' * rate_residuals) < 1e-6
        end

        # The degenerate cases are stated on the design matrix, since the solve itself no
        # longer tests rank. Two distinct directions (a satellite tracked on a second band
        # repeats its line of sight) leave the position columns rank deficient.
        solvable(dirs) = !isnothing(PositionVelocityTime.positive_definite_cholesky(
            Symmetric(design(dirs)' * design(dirs))))
        @test !solvable([[1.0, 0.2, 0.9], [1.0, 0.2, 0.9], [0.3, -1.0, 0.5], [0.3, -1.0, 0.5]])

        # Three independent directions on a common cone (equal projection onto z) put the
        # drift column in their span, leaving the fourth pivot exactly zero — the shape
        # that produced SingularException(4) from the 4×4 solve.
        cosθ = 0.5
        sinθ = sqrt(1 - cosθ^2)
        @test !solvable([
            [sinθ, 0.0, cosθ], [0.0, sinθ, cosθ], [-sinθ, 0.0, cosθ], [0.0, -sinθ, cosθ]])
    end
end
