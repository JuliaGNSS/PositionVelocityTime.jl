# BeiDou (BDS) support: the legacy D1/D2 message on B1I and B3I, and the BDS-3
# B-CNAV1/2/3 messages on B1C, B2a and B2b.
#
# Four things are genuinely new here rather than shared with GPS/Galileo, and the
# testsets below are organised around them:
#
#   1. the clock and ephemeris accessors, which have to read a BeiDou container as
#      readily as a GPS or Galileo one for the shared clock model to apply;
#   2. the B-CNAV quasi-Keplerian flavour — `A_REF` selected by `sat_type`, and an
#      outright `Ω̇` rather than GPS CNAV's delta off a reference rate;
#   3. the GEO reference frame, which no GPS or Galileo satellite exercises;
#   4. BDT as a fourth receiver-clock time system, with the BGTO constraint toward GPS.

@testset "BeiDou (BDS) support" begin
    C = 299792458.0
    kw = (;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )

    # ---- Fixture builders ---------------------------------------------------
    #
    # The Keplerian elements come from the Galileo E1B fixtures, which are real
    # ephemerides with geometry good enough to converge. Rebuilding them as BeiDou
    # messages is legitimate — the element set is the same shape — and it keeps the
    # BeiDou tests from depending on a second hand-written constellation. The
    # satellite positions do differ slightly from the Galileo ones, because BDCS's
    # Earth rotation rate is 7.2921150e-5 rad/s against GTRF's 7.2921151467e-5, so
    # nothing below asserts that the two agree.

    "A BeiDou D1 (B1I/B3I) decoder carrying `state`'s Galileo Keplerian ephemeris."
    function dnav_decoder(state; prn = state.decoder.prn, T_GD1 = -2.9e-9, extra...)
        d = state.decoder.data
        data = GNSSDecoder.BeiDouDNAVData(;
            SOW = Int64(d.TOW),
            WN = 800,
            SatH1 = false,
            AODC = Int64(1),
            URAI = Int64(2),
            AODE = Int64(3),
            t_0c = Int64(d.t_0c),
            a_f0 = d.a_f0,
            a_f1 = d.a_f1,
            a_f2 = d.a_f2,
            T_GD1 = T_GD1,
            T_GD2 = 3.7e-9,
            t_0e = Int64(d.t_0e),
            sqrt_A = d.sqrt_A,
            e = d.e,
            ω = d.ω,
            Δn = d.Δn,
            M_0 = d.M_0,
            Ω_0 = d.Ω_0,
            Ω_dot = d.Ω_dot,
            i_0 = d.i_0,
            i_dot = d.i_dot,
            C_uc = d.C_uc,
            C_us = d.C_us,
            C_rc = d.C_rc,
            C_rs = d.C_rs,
            C_ic = d.C_ic,
            C_is = d.C_is,
            extra...,
        )
        GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB1IDecoderState(prn);
            data = data,
            raw_data = data,
        )
    end

    "A BeiDou B2a (B-CNAV2) decoder carrying `state`'s ephemeris in quasi-Keplerian form."
    function b2a_decoder(state; sat_type = 3, extra...)
        d = state.decoder.data
        A_ref = sat_type == 3 ? 27_906_100.0 : 42_162_200.0
        data = GNSSDecoder.BeiDouB2aData(;
            SOW = Int64(d.TOW),
            WN = Int64(800),
            HS = Int64(0),
            IODE = Int64(0x2c),
            IODC = Int64(0x32c),           # low 8 bits must equal IODE (§7.4.3)
            t_0e = Int64(d.t_0e),
            sat_type = Int64(sat_type),
            ΔA = d.sqrt_A^2 - A_ref,
            A_dot = 0.0,
            Δn_0 = d.Δn,
            Δn_0_dot = 0.0,
            M_0 = d.M_0,
            e = d.e,
            ω = d.ω,
            SOW_mt10 = Int64(d.TOW),
            Ω_0 = d.Ω_0,
            i_0 = d.i_0,
            Ω_dot = d.Ω_dot,
            i_dot = d.i_dot,
            C_is = d.C_is,
            C_ic = d.C_ic,
            C_rs = d.C_rs,
            C_rc = d.C_rc,
            C_us = d.C_us,
            C_uc = d.C_uc,
            SOW_mt11 = Int64(d.TOW) + 3,   # adjacent frames pair MT10 with MT11
            t_0c = Int64(d.t_0c),
            a_f0 = d.a_f0,
            a_f1 = d.a_f1,
            a_f2 = d.a_f2,
            T_GD_B2ap = -4.1e-9,
            ISC_B2ad = 1.3e-9,
            T_GD_B1Cp = -2.2e-9,
            extra...,
        )
        GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB2aDecoderState(state.decoder.prn);
            data = data,
            raw_data = data,
        )
    end

    """
    A BeiDou B1C (B-CNAV1) decoder carrying `state`'s ephemeris in quasi-Keplerian
    form. B1C rides on L1, the same band GPS L1 C/A and Galileo E1 use, so it is what
    a mixed-constellation epoch needs to share a band rather than split into
    disconnected coverage components.
    """
    function b1c_decoder(state; sat_type = 3, extra...)
        d = state.decoder.data
        A_ref = sat_type == 3 ? 27_906_100.0 : 42_162_200.0
        # B-CNAV1 sends the hour of week and the 18 s seconds-of-hour count instead of
        # a seconds-of-week field, so the fixture's time of week has to be split.
        sow = Int64(d.TOW)
        data = GNSSDecoder.BeiDouB1CData(;
            HOW = sow ÷ 3600,
            soh = Int((sow % 3600) ÷ 18),
            WN = Int64(800),
            HS = Int64(0),
            IODE = Int64(0x2c),
            IODC = Int64(0x32c),
            t_0e = Int64(d.t_0e),
            sat_type = Int64(sat_type),
            ΔA = d.sqrt_A^2 - A_ref,
            A_dot = 0.0,
            Δn_0 = d.Δn,
            Δn_0_dot = 0.0,
            M_0 = d.M_0,
            e = d.e,
            ω = d.ω,
            Ω_0 = d.Ω_0,
            i_0 = d.i_0,
            Ω_dot = d.Ω_dot,
            i_dot = d.i_dot,
            C_is = d.C_is,
            C_ic = d.C_ic,
            C_rs = d.C_rs,
            C_rc = d.C_rc,
            C_us = d.C_us,
            C_uc = d.C_uc,
            t_0c = Int64(d.t_0c),
            a_f0 = d.a_f0,
            a_f1 = d.a_f1,
            a_f2 = d.a_f2,
            T_GD_B1Cp = -2.2e-9,
            ISC_B1Cd = 0.9e-9,
            T_GD_B2ap = -4.1e-9,
            extra...,
        )
        GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB1CDecoderState(state.decoder.prn);
            data = data,
            raw_data = data,
        )
    end

    galileo = galileo_e1b_states(0.0Hz)
    d1 = dnav_decoder(galileo[1])
    b2a = b2a_decoder(galileo[1])

    # ---- Classification -----------------------------------------------------

    @testset "classification: BDT, bands, absolute week, epoch" begin
        for signal in (BeiDouB1I(), BeiDouB3I(), BeiDouB1C_D(), BeiDouB2aI(), BeiDouB2bI())
            @test GNSSSignals.get_time_system(signal) == BDT()
        end
        # B1I and B3I sit in bands of their own; B1C shares L1 with GPS/Galileo and
        # B2a shares L5, which is what lets a mixed-constellation epoch share one
        # inter-frequency-bias column per band rather than one per constellation.
        @test GNSSSignals.get_band_id(BeiDouB1I()) == :B1I
        @test GNSSSignals.get_band_id(BeiDouB3I()) == :B3I
        @test GNSSSignals.get_band_id(BeiDouB1C_D()) == GNSSSignals.get_band_id(GPSL1CA())
        @test GNSSSignals.get_band_id(BeiDouB2aI()) == GNSSSignals.get_band_id(GPSL5I())
        @test GNSSSignals.get_band_id(BeiDouB2bI()) == GNSSSignals.get_band_id(GalileoE5bI())

        # The BDT epoch is 2006-01-01, 26 years after GPS's, so the same week/SOW pair
        # dates to a different instant — this is what keeps BeiDou's reported time
        # honest without any 14 s fudge in the solver.
        @test GNSSSignals.get_system_start_time(BDT()) == DateTime(2006, 1, 1)

        # 13-bit BDT week, broadcast absolute: no 1024-week rollover to resolve, so
        # `approximate_year` cannot change the answer.
        @test GNSSDecoder.get_week(d1; approximate_year = 2024) == 800
        @test GNSSDecoder.get_week(d1; approximate_year = 1999) == 800
        @test GNSSDecoder.get_week(b2a; approximate_year = 2024) == 800
    end

    @testset "time of week comes from the decoder, on every BeiDou signal" begin
        # GNSSDecoder answers this now, including B-CNAV1's reconstruction. Asserted
        # here because `calc_uncorrected_transmit_time` reads it for every signal, so a
        # `nothing` would propagate silently into a transmit time.
        @test GNSSDecoder.get_time_of_week(d1) == d1.data.SOW
        @test GNSSDecoder.get_time_of_week(b2a) == b2a.data.SOW

        # B-CNAV1 alone broadcasts no seconds of week: subframe 2 carries the hour of
        # week and subframe 1 the 18 s seconds-of-hour count
        # (BDS-SIS-ICD-B1C-1.0 §7.3), so SOW = HOW·3600 + SOH·18.
        b1c_data = GNSSDecoder.BeiDouB1CData(; HOW = Int64(36), soh = 25)
        b1c = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB1CDecoderState(20);
            data = b1c_data,
            raw_data = b1c_data,
        )
        @test GNSSDecoder.get_time_of_week(b1c) == 36 * 3600 + 25 * 18
    end

    # ---- Clock polynomial ----------------------------------------------------

    @testset "the clock model reads the BeiDou containers like any other" begin
        # The clock model and the propagator are written once for all eleven
        # signals, directly off the field names every GNSSDecoder container
        # spells the same way — pin the polynomial they build against the ICD
        # expression on a BeiDou container.
        t = 132769.0
        Δt = t - d1.data.t_0c
        @test PositionVelocityTime.calc_satellite_clock_drift(ephemeris(d1), t) ≈
              d1.data.a_f1 + 2 * d1.data.a_f2 * Δt
        @test PositionVelocityTime.correct_clock(ephemeris(d1), GNSSDecoder.get_group_delay(d1, BeiDouB3I()), t) ≈
              t - (
            d1.data.a_f0 + d1.data.a_f1 * Δt + d1.data.a_f2 * Δt^2 +
            PositionVelocityTime.calc_relativistic_correction(ephemeris(d1), t)
        )
    end

    # ---- B-CNAV quasi-Keplerian elements ------------------------------------

    @testset "B-CNAV elements: A_REF by sat_type, Ω̇ broadcast outright" begin
        # GNSSDecoder normalises the elements (`get_broadcast_ephemeris`); the
        # propagator's result depends on its doing so right.
        el = ephemeris(b2a)

        # MEO reference (27 906 100 m), not GPS CNAV's 26 559 710 m.
        @test el.A ≈ 27_906_100.0 + b2a.data.ΔA

        # An IGSO satellite with the same ΔA sits on the *other* reference — a
        # 14 000 km difference, so this cannot be a rounding-level distinction.
        igso = b2a_decoder(galileo[1]; sat_type = 2)
        @test ephemeris(igso).A ≈ 42_162_200.0 + igso.data.ΔA

        # Ω̇ is the broadcast value itself. GPS CNAV would add Ω̇_REF = −2.6e-9·π
        # here; doing that for BeiDou would be a ~8e-9 rad/s error, tens of metres of
        # along-track position within one fit interval.
        @test el.Ω_dot == b2a.data.Ω_dot
        @test !isapprox(el.Ω_dot, b2a.data.Ω_dot + -2.6e-9 * π; atol = 1e-13)

        # Ȧ and Δṅ_0 enter the same way as in GPS CNAV.
        accelerating = b2a_decoder(galileo[1]; A_dot = 0.5, Δn_0_dot = 1.0e-14)
        el_t = ephemeris(accelerating)
        @test el_t.A_dot == 0.5
        @test PositionVelocityTime.mean_motion(el_t, 100.0) ≈
              sqrt(el_t.μ / el_t.A^3) + accelerating.data.Δn_0 + 0.5 * 1.0e-14 * 100.0
    end

    # ---- Group delay --------------------------------------------------------
    # The group delay a range takes — B3I-referenced on every BeiDou message, the data
    # components' ISCs, B2a's ungated terms — is GNSSDecoder's `get_group_delay`, tested
    # there; `test/signals/galileo_e5b_e6b.jl` checks that its method table covers
    # exactly the pairings a solve can produce.

    # ---- The GEO reference frame --------------------------------------------

    @testset "GEO satellites use the tilted, rotated BDS frame" begin
        # GNSSDecoder's `get_orbit_class` decides the GEO branch: D1/D2
        # answers from the PRN partition (GEO = 1–5 and 59–63), B-CNAV from the
        # broadcast `sat_type`, GPS and Galileo from the constellation supertype.
        # B1I/B3I report `nothing` for a non-GEO satellite rather than claiming MEO —
        # the PRN partition separates GEO from non-GEO and no further — which is why
        # the predicate tests *for* GEO and not against MEO.
        @test GNSSDecoder.get_orbit_class(dnav_decoder(galileo[1]; prn = 3)) ===
              GNSSDecoder.geostationary_orbit
        @test isnothing(GNSSDecoder.get_orbit_class(dnav_decoder(galileo[1]; prn = 20)))
        @test GNSSDecoder.get_orbit_class(b2a) === GNSSDecoder.medium_earth_orbit
        @test ephemeris(dnav_decoder(galileo[1]; prn = 3)).is_geostationary
        @test ephemeris(dnav_decoder(galileo[1]; prn = 60)).is_geostationary
        @test !ephemeris(dnav_decoder(galileo[1]; prn = 20)).is_geostationary
        @test ephemeris(b2a_decoder(galileo[1]; sat_type = 1)).is_geostationary
        @test !ephemeris(b2a).is_geostationary
        # No GPS or Galileo satellite ever takes this branch.
        @test !ephemeris(galileo[1].decoder).is_geostationary

        # The GEO branch must be exactly "the standard projection, but in a frame
        # where Ω_k does not absorb the earth rotation, then R_z(ω̇_e·t_k)·R_x(−5°)".
        # That is checkable without re-deriving the propagator: feeding the *non-GEO*
        # path an `Ω̇` raised by ω̇_e cancels the `−ω̇_e·t_k` it folds into Ω_k, so it
        # produces precisely the auxiliary-frame vector X_GK the ICD's rotations act
        # on. The rotation matrices below are written out from the ICD definitions,
        # independently of `_rotate_beidou_geo`.
        ω_e = 7.2921150e-5
        elements = (; Ω_dot = -5.3e-9, t_0e = Int64(0))
        geo_frame = dnav_decoder(galileo[1]; prn = 3, elements...)
        auxiliary = dnav_decoder(
            galileo[1]; prn = 20, elements..., Ω_dot = elements.Ω_dot + ω_e)
        icd_R_z(φ) = [cos(φ) sin(φ) 0.0; -sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]
        icd_R_x(φ) = [1.0 0.0 0.0; 0.0 cos(φ) sin(φ); 0.0 -sin(φ) cos(φ)]

        for t in (0.0, 1800.0, 21600.0)
            pv_geo = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(geo_frame), t)
            pv_aux = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(auxiliary), t)
            rotation = icd_R_z(ω_e * t) * icd_R_x(-5 * π / 180)
            @test pv_geo.position ≈ rotation * pv_aux.position atol = 1e-6
            # Rotations preserve length, and the two frames differ by a real 5° tilt.
            @test norm(pv_geo.position) ≈ norm(pv_aux.position)
            @test norm(pv_geo.position - pv_aux.position) > 1e5

            # The velocity is where a GEO implementation most easily goes wrong: the
            # frame itself rotates, so the product rule contributes a term worth
            # kilometres per second. Check it against a central difference of the
            # position, which shares none of the derivative code.
            h = 1e-3
            forward = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(geo_frame), t + h).position
            backward = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(geo_frame), t - h).position
            @test pv_geo.velocity ≈ (forward - backward) / (2h) atol = 1e-4
        end

        # And the point of the whole construction: a genuinely geostationary satellite
        # comes out fixed in the earth-fixed frame. In the auxiliary frame such an
        # orbit is *not* equatorial — it is inclined by the same 5° the frame is
        # tilted by, with its ascending node at Ω = π so the two tilts cancel — which
        # is exactly why the naive "i = 0 means equatorial" reading of a BDS GEO
        # ephemeris is wrong. Mean motion is tuned to one revolution per sidereal day.
        μ = 3.986004418e14
        stationary = dnav_decoder(
            galileo[1];
            prn = 3,
            sqrt_A = (μ / ω_e^2)^(1 / 6),
            e = 0.0, ω = 0.0, M_0 = 0.0, Δn = 0.0,
            i_0 = 5 * π / 180, i_dot = 0.0,
            Ω_0 = π, Ω_dot = 0.0,
            C_uc = 0.0, C_us = 0.0, C_rc = 0.0, C_rs = 0.0, C_ic = 0.0, C_is = 0.0,
            t_0e = Int64(0),
        )
        anchor = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(stationary), 0.0)
        @test norm(anchor.position) ≈ (μ / ω_e^2)^(1 / 3)
        for t in (3600.0, 21600.0, 43082.0)
            later = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(stationary), t)
            # Fixed to well under a metre over half a day, and essentially at rest —
            # against the ~3 km/s the satellite moves at in the auxiliary frame.
            @test norm(later.position - anchor.position) < 1.0
            @test norm(later.velocity) < 1e-3
        end
    end

    # ---- Ionosphere ---------------------------------------------------------

    @testset "legacy Klobuchar set is BeiDou's own variant" begin
        # The two message families broadcast two different models, so each container
        # answers for exactly one of them: only the legacy D1/D2 one yields a
        # Klobuchar set, and only the B-CNAV ones yield BDGIM (see
        # `test/atmosphere/` for the models themselves).
        @test GNSSDecoder.get_klobuchar_parameters(d1) === nothing
        @test GNSSDecoder.get_klobuchar_parameters(b2a) === nothing
        @test GNSSDecoder.get_bdgim_parameters(d1) === nothing

        iono = dnav_decoder(
            galileo[1];
            α_0 = 1.0e-8, α_1 = 2.0e-8, α_2 = -3.0e-8, α_3 = 4.0e-8,
            β_0 = 90112.0, β_1 = 16384.0, β_2 = -196608.0, β_3 = 131072.0,
        )
        p = GNSSDecoder.get_klobuchar_parameters(iono)
        # A type of its own, not `KlobucharParameters`: BDS-SIS-ICD-B1I-3.0 §5.2.4.7
        # prescribes BeiDou's own algorithm (geographic pierce-point latitude,
        # exact geometry, a period clamped from above too), with the delay defined
        # at B1I rather than L1 (see `test/atmosphere/klobuchar.jl` for the model itself).
        @test p isa GNSSDecoder.BeiDouKlobucharParameters
        @test !(p isa GNSSDecoder.KlobucharParameters)
        @test (p.α_0, p.β_3) == (1.0e-8, 131072.0)

        geometry = (0.6, 1.2, Geodesy.LLA(50.8, 6.1, 180.0), bdt_epoch(132000.0))
        bds = PositionVelocityTime.ionospheric_delay(p, carrier_hz(BeiDouB1I()), geometry...)
        @test bds > 0
        # The same eight numbers fed through the GPS algorithm land elsewhere: the
        # model differs, and so does the reference carrier the 1/f² rescale starts
        # from (B1I is 14 MHz below L1, ~1.8 % on its own).
        as_gps_set = GNSSDecoder.KlobucharParameters(
            p.α_0, p.α_1, p.α_2, p.α_3, p.β_0, p.β_1, p.β_2, p.β_3)
        gps = PositionVelocityTime.ionospheric_delay(as_gps_set, carrier_hz(BeiDouB1I()), geometry...)
        @test !isapprox(bds, gps; rtol = 1e-3)
    end

    # ---- BGTO ---------------------------------------------------------------

    @testset "the BDT second-of-week counts 14 s behind a GPS time of week" begin
        # Structural, from the time scales' own definitions rather than any broadcast
        # value: GPST and GST are both TAI − 19 s and so count alike, BDT is TAI − 33 s.
        # Per satellite it is the seconds to add to reach the *reference* system's count,
        # applied whether or not a clock is constrained. Only each measurement's
        # `time_system` is read, so a stand-in carrying only that field is enough here.
        # GPS-referenced: BeiDou needs +14 s.
        stand_in(system) = (; time_system = system)
        offsets(systems, reference) = map(
            system -> PositionVelocityTime.calc_leap_second_correction(stand_in(system),
                reference),
            systems)
        @test offsets([GPST(), BDT(), GST(), BDT()], GPST()) == [0.0, 14.0, 0.0, 14.0]
        # BeiDou-referenced: the sign reverses, and a BeiDou-only epoch needs nothing —
        # which is why anchoring on GPST instead of the reference system would have
        # traded this bug for a 14 s error in the reported time of every BDT-referenced
        # fix, where `reference_time` is dated from the BDT epoch.
        @test offsets([BDT(), GPST(), BDT()], BDT()) == [0.0, -14.0, 0.0]
        @test all(iszero, offsets([BDT(), BDT(), BDT()], BDT()))
    end

    @testset "BDT–GPS time offset (BGTO) from all three message families" begin
        # PVT no longer parses these shapes — GNSSDecoder's `get_broadcast_time_offset` does —
        # but it still reads them, so each family is exercised through the two
        # functions the solve calls.
        #
        # GNSSDecoder reports `A_0` as broadcast, without the 14 s between the BDT and
        # GPST counts, which `calc_leap_second_correction` applies to the transmit
        # times instead; the steering offset is the polynomial alone.
        to_gpst(decoder, t) = steering(decoder, GPST(), t)

        @test !has_offset(d1, GPST())
        @test !has_offset(b2a, GPST())

        # Legacy D1: a bare two-term polynomial with no reference epoch, so the
        # argument is the time of week itself (BDS-SIS-ICD-B1I-3.0 §5.2.4.18).
        legacy = dnav_decoder(galileo[1]; A_0GPS = 8.0e-9, A_1GPS = 1.0e-15)
        @test has_offset(legacy, GPST())
        @test to_gpst(legacy, 132000.0) ≈ 8.0e-9 + 1.0e-15 * 132000.0

        # B-CNAV: three terms against an explicit (WN_0BGTO, t_0BGTO) epoch, and only
        # while GNSS_ID names GPS. GNSS ID 1 is GPS; 0 means "no parameters", which is
        # the trap — reading it as GPS would publish a zero offset as if measured.
        bgto = (; WN_0BGTO = Int64(799), t_0BGTO = Int64(100),
            A_0BGTO = 5.0e-9, A_1BGTO = 1.0e-15, A_2BGTO = 1.0e-21)
        gps_bgto = b2a_decoder(galileo[1]; GNSS_ID = Int64(1), bgto...)
        @test has_offset(gps_bgto, GPST())
        Δτ = 132000.0 - 100 + 604800 * (800 - 799)
        @test to_gpst(gps_bgto, 132000.0) ≈
              5.0e-9 + 1.0e-15 * Δτ + 1.0e-21 * Δτ^2
        for gnss_id in (Int64(0), Int64(2), Int64(3))   # unavailable / Galileo / GLONASS
            @test !has_offset(
                b2a_decoder(galileo[1]; GNSS_ID = gnss_id, bgto...), GPST())
        end

        # B-CNAV3 uses the same fields ...
        b2b_data = GNSSDecoder.BeiDouB2bData(;
            WN = Int64(800), GNSS_ID = Int64(1), WN_0BGTO = Int64(799),
            t_0BGTO = Int64(100), A_0BGTO = 5.0e-9, A_1BGTO = 0.0, A_2BGTO = 0.0)
        b2b = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB2bDecoderState(20); data = b2b_data, raw_data = b2b_data)
        @test has_offset(b2b, GPST())
        @test to_gpst(b2b, 132000.0) ≈ 5.0e-9

        # ... while B-CNAV1 pages one record per GNSS into a dictionary, so the GPS
        # one is looked up by key instead of gated on a field.
        b1c_data = GNSSDecoder.BeiDouB1CData(;
            WN = Int64(800),
            bgtos = Dictionary(
                [1, 2],
                [
                    GNSSDecoder.BeiDouB1CBGTO(; GNSS_ID = 1, WN_0BGTO = 799,
                        t_0BGTO = 100, A_0BGTO = 5.0e-9, A_1BGTO = 0.0, A_2BGTO = 0.0),
                    GNSSDecoder.BeiDouB1CBGTO(; GNSS_ID = 2, WN_0BGTO = 799,
                        t_0BGTO = 100, A_0BGTO = -7.0e-9, A_1BGTO = 0.0, A_2BGTO = 0.0),
                ],
            ),
        )
        b1c = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.BeiDouB1CDecoderState(20); data = b1c_data, raw_data = b1c_data)
        @test has_offset(b1c, GPST())
        # The GPS record, not the Galileo one that also sits in the dictionary.
        @test to_gpst(b1c, 132000.0) ≈ 5.0e-9
        # A satellite carrying only a Galileo offset has nothing to convert to GPS.
        galileo_only = GNSSDecoder.BeiDouB1CData(b1c_data;
            bgtos = Dictionary([2], [b1c_data.bgtos[2]]))
        @test !has_offset(
            GNSSDecoder.GNSSDecoderState(b1c; data = galileo_only, raw_data = galileo_only),
            GPST())
    end

    # ---- End to end ---------------------------------------------------------
    #
    # Synthesize observables from a known user position by solving the light-time
    # equation per satellite, then check `calc_pvt` recovers that position. This
    # exercises the whole BeiDou chain — SOW, the shimmed clock model, the
    # propagator, BDT as a time system — against an answer computed independently of
    # the solver.

    "Set a decoder's bit count and code phase so its corrected transmit time is `t_tx`."
    function observables_for(decoder, system, t_tx)
        # Undo the clock correction to recover the uncorrected time the decoder must
        # report. The correction changes by ~1e-12 s over its own size, so a couple of
        # fixed-point steps converge to well below numerical noise.
        t_uncorrected = t_tx
        for _ = 1:4
            t_uncorrected +=
                t_tx - PositionVelocityTime.correct_clock(ephemeris(decoder), GNSSDecoder.get_group_delay(decoder, system), t_uncorrected)
        end
        bit_rate = Float64(GNSSSignals.get_bit_rate(decoder) / Hz)
        chip_rate = Float64(GNSSSignals.get_chip_rate(system) / Hz)
        elapsed = t_uncorrected - GNSSDecoder.get_time_of_week(decoder)
        num_bits = floor(Int, elapsed * bit_rate)
        code_phase = (elapsed - num_bits / bit_rate) * chip_rate
        (num_bits, code_phase)
    end

    "A `SatelliteState` whose observables place the receiver at `user` at time `t_rx`."
    function ranging_state(decoder, system, user, t_rx)
        # Solve the light-time equation against the *same* range model the solver
        # uses, which rotates the satellite into the ECEF frame of the reception
        # instant (`predict_pseudoranges!`). Omitting that Sagnac term would leave a ±30 m
        # per-satellite inconsistency and the recovered position tens of metres out —
        # a synthesis error that looks exactly like a decoder bug.
        t_tx = t_rx - 0.075                      # ~one light-time as a seed
        for _ = 1:6
            position = PositionVelocityTime.calc_satellite_position_and_velocity(ephemeris(decoder), t_tx).position
            rotated = PositionVelocityTime.satellite_position_at_reception(user, position)
            t_tx = t_rx - norm(rotated - user) / C
        end
        num_bits, code_phase = observables_for(decoder, system, t_tx)
        SatelliteState(;
            decoder = GNSSDecoder.GNSSDecoderState(
                decoder; num_symbols_since_sync = num_bits),
            signal = system,
            code_phase = code_phase,
            carrier_doppler = 0.0Hz,
            carrier_phase = 0.0,
        )
    end

    user = ECEF(4.0e6, 0.6e6, 4.8e6)
    t_rx = 132769.0

    # Tolerance note: the pseudoranges are differences of times of week, ~1.3e5 s
    # apart in Float64, so one ULP of a transmit time is ~1.5e-11 s — about 4 mm of
    # range. A perfectly consistent synthesis therefore still lands a few millimetres
    # off, which is why the position assertions below are at the centimetre level, as
    # everywhere else in this suite.

    @testset "BeiDou-only fix recovers a synthesized position" begin
        b1i = [
            ranging_state(dnav_decoder(state; prn = 20 + i), BeiDouB1I(), user, t_rx)
            for (i, state) in enumerate(galileo)
        ]
        pvt = calc_pvt(signal_group(b1i); kw...)
        @test length(pvt.satellites) == length(galileo)
        @test pvt.reference_system == BDT()
        @test norm(pvt.position - user) < 1e-2
        # One constellation, one band, so no inter-system or inter-frequency biases.
        @test isempty(pvt.inter_system_biases)
        @test isempty(pvt.inter_frequency_biases)
        # The reported epoch is dated from BDT's own 2006 origin. Measuring the gap
        # back to that same origin only checks the arithmetic, so pin the *absolute*
        # instant against GPS instead: build a GPS constellation observing the very
        # same instant and require the two solutions to report one time.
        #
        # BDT week 0 second 0 is GPS week 1356, time of week 14 — the BDT epoch is
        # 2006-01-01 00:00:00 UTC, which is 14 s into GPS week 1356 because 14 leap
        # seconds had accrued since the GPS epoch. So the same instant is BDT
        # (800, t_rx) and GPST (2156, t_rx + 14); the broadcast GPS week is 10-bit, so
        # its field carries 2156 mod 1024 and `get_week` reconstructs the rest from
        # `approximate_year`. Agreement exercises what the circular check could not:
        # that `get_week_start_epoch` places both origins on the atomic scale correctly,
        # BDT's 33 s TAI offset against GPST's 19 s included.
        @test !isnothing(pvt.time)
        gps_matched = [
            ranging_state(
                GNSSDecoder.GNSSDecoderState(
                    state.decoder;
                    data = GNSSDecoder.GPSL1CAData(
                        state.decoder.data; WN = 2156 % 1024, TOW = Int64(t_rx) + 14),
                ),
                GPSL1CA(),
                user,
                t_rx + 14,
            ) for state in gps_l1_states(0.0Hz)
        ]
        pvt_gps = calc_pvt(signal_group(gps_matched); kw...)
        @test GNSSDecoder.get_week(gps_matched[1].decoder;
            approximate_year = 2021) == 2156
        @test norm(pvt_gps.position - user) < 1e-2
        # One instant, two constellations, two epochs, two TAI offsets: one answer.
        @test pvt_gps.time - pvt.time ≈ 0.0 atol = 1e-3

        # A B3I range off the same satellites is a second band, so it gains an IFB
        # column but must not move the position: B3I needs no group-delay correction
        # while B1I subtracts T_GD1, and the observables above account for each.
        b3i = [
            ranging_state(dnav_decoder(state; prn = 20 + i), BeiDouB3I(), user, t_rx)
            for (i, state) in enumerate(galileo)
        ]
        two_band = calc_pvt((signal_group(b1i), signal_group(b3i)); kw...)
        @test norm(two_band.position - user) < 1e-2
        @test collect(keys(two_band.inter_frequency_biases)) == [:B3I]
        @test abs(two_band.inter_frequency_biases[:B3I].value) < 1e-2m
    end

    @testset "a mixed fix across the GPS week wrap is dated in the right week" begin
        # BDT counts 14 s behind GPST, so GPS week 2157, 5 s is BDT week 800, 604791 s.
        # Seen by GPS and BeiDou together, the BeiDou transmit times on the GPS count
        # (SOW + 14) lie past 604800 s while the GPS ones have wrapped to a few seconds;
        # the reported time must still be GPS week 2157, 5 s.
        # GPS has more satellites, so it is the reference system.
        gps = [
            ranging_state(
                GNSSDecoder.GNSSDecoderState(
                    state.decoder;
                    data = GNSSDecoder.GPSL1CAData(
                        state.decoder.data; WN = 2157 % 1024, TOW = Int64(0)),
                ),
                GPSL1CA(),
                user,
                5.0,
            ) for state in gps_l1_states(0.0Hz)
        ]
        bds = [
            ranging_state(
                let decoder = dnav_decoder(state; prn = 20 + i)
                    data = GNSSDecoder.BeiDouDNAVData(decoder.data; SOW = Int64(604790))
                    GNSSDecoder.GNSSDecoderState(decoder; data, raw_data = data)
                end,
                BeiDouB1I(),
                user,
                604791.0,
            ) for (i, state) in enumerate(galileo[1:4])
        ]
        expected = get_week_start_epoch(GPST(), 2157) + 5.0
        gps_only = calc_pvt(signal_group(gps); kw...)
        @test gps_only.time - expected ≈ 0.0 atol = 1e-3
        mixed = calc_pvt((signal_group(gps), signal_group(bds)); kw...)
        @test mixed.reference_system == GPST()
        @test norm(mixed.position - user) < 1e-2
        @test mixed.time - expected ≈ 0.0 atol = 1e-3
    end

    @testset "ranging on the pilot recovers the same fix as ranging on the data" begin
        # This is the configuration a joint pilot/data tracking receiver actually
        # runs: both components tracked off one carrier/code loop, ranging on the
        # dataless pilot (no navigation-bit transitions to squash the discriminator)
        # and the message decoded from the data component. Such a receiver hands
        # `calc_pvt` the *pilot* as `system` and the *data* decoder separately, so
        # the pilot path is the common case rather than an exotic one.
        #
        # Two things differ from the data-ranged fix and both are worth a number
        # rather than a `hasmethod` check. `calc_uncorrected_transmit_time` divides the
        # pilot's code frequency by the *decoder's* data frequency to reduce the
        # shared code phase modulo one bit — for B1C that is 1.023 Mcps over
        # 100 Hz, 10230 chips per bit — and GNSSDecoder's `get_group_delay` gives the
        # pilot's own term, `T_GD_B1Cp` without the data component's `ISC_B1Cd`. The
        # observables are built against both, so the two fixes must agree.
        for (decoder_of, data_signal, pilot) in (
            (b1c_decoder, BeiDouB1C_D(), BeiDouB1C_P()),
            (b2a_decoder, BeiDouB2aI(), BeiDouB2aQ()),
        )
            on_data = [
                ranging_state(decoder_of(state), data_signal, user, t_rx)
                for state in galileo
            ]
            on_pilot = [
                ranging_state(decoder_of(state), pilot, user, t_rx) for state in galileo
            ]
            # The pilot really is a different measurement: its group delay differs
            # from the data component's by the ISC, so the observables are not copies.
            @test any(
                !isapprox(d.code_phase, p.code_phase) for (d, p) in zip(on_data, on_pilot)
            )

            pvt_data = calc_pvt(signal_group(on_data); kw...)
            pvt_pilot = calc_pvt(signal_group(on_pilot); kw...)
            @test length(pvt_pilot.satellites) == length(galileo)
            @test norm(pvt_pilot.position - user) < 1e-2
            @test norm(pvt_pilot.position - pvt_data.position) < 1e-2
            # Same satellite, same band — so a pilot-ranged epoch is one band and one
            # time system, with no bias columns of its own.
            @test isempty(pvt_pilot.inter_frequency_biases)
            @test pvt_pilot.reference_system == BDT()
        end
    end

    @testset "B-CNAV fix, and the BDS-3 quasi-Keplerian path end to end" begin
        b2a_states = [
            ranging_state(b2a_decoder(state), BeiDouB2aI(), user, t_rx)
            for state in galileo
        ]
        pvt = calc_pvt(signal_group(b2a_states); kw...)
        @test length(pvt.satellites) == length(galileo)
        @test norm(pvt.position - user) < 1e-2
    end

    @testset "GPS + BeiDou: separate clocks, and the BGTO constraint when scarce" begin
        # Synthesize both constellations against the same user position and the same
        # physical instant. Two things separate their reported times, and getting only
        # the small one right is how this test used to encode a bug rather than catch
        # it: it fed both constellations the same numeric reception time, which is only
        # correct if BDT and GPST count alike. They do not.
        #
        #   - Structural, 14 s: BDT is TAI − 33 s against GPST's TAI − 19 s, so BDT
        #     week 0 second 0 is GPS week 1356 time of week 14. A BeiDou
        #     second-of-week therefore reads 14 s *lower* than the GPS time-of-week
        #     for the same instant, which makes its raw pseudorange 4.2e9 m too long.
        #   - Broadcast steering, Δ: the tens-of-nanoseconds residual the BGTO
        #     publishes, which is what an inter-system bias should actually report.
        #
        # So the BeiDou reception time is `t_rx - 14 + Δ`. Without
        # `calc_leap_second_correction` the 14 s survives into the normal equations
        # as a parameter nine orders of magnitude above the others: the separate-clock
        # solve reports a nonsense inter-system bias, and the constrained solve — which
        # pins the BeiDou clock difference to the steering term — fails outright.
        Δ = 12.0e-9
        structural = -14.0
        gps = [
            ranging_state(state.decoder, GPSL1CA(), user, t_rx)
            for state in gps_l1_states(0.0Hz)
        ]
        beidou(prns; bgto = nothing) = [
            ranging_state(
                dnav_decoder(
                    galileo[i];
                    prn = prns[i],
                    (isnothing(bgto) ? (;) : (; A_0GPS = bgto, A_1GPS = 0.0))...,
                ),
                BeiDouB1I(),
                user,
                t_rx + structural + Δ,
            ) for i in eachindex(prns)
        ]

        # Enough satellites for two independent clocks. GPS lives only on L1 and
        # BeiDou only on B1I, so the band-coverage graph is disconnected and neither
        # band gets an inter-frequency-bias column — each constellation's chain delay
        # folds into its own clock. The inter-system bias then reads out as −c·Δ.
        plenty_beidou = beidou([21, 22, 23, 24])
        plenty = (signal_group(gps), signal_group(plenty_beidou))
        pvt = calc_pvt(plenty; kw...)
        @test length(pvt.satellites) == length(gps) + length(plenty_beidou)
        @test norm(pvt.position - user) < 1e-2
        @test haskey(pvt.inter_system_biases, BDT())
        # The steering term alone, ~3.6 m — not the 4.2e9 m the structural count
        # difference would contribute if it leaked into the readout.
        @test ustrip(m, pvt.inter_system_biases[BDT()]) ≈ -C * Δ atol = 1e-2
        @test abs(ustrip(m, pvt.inter_system_biases[BDT()])) < 100

        # Scarce, and on a shared band. B1C rides on L1 with GPS L1 C/A, so the
        # coverage graph is connected and no inter-frequency-bias column is needed —
        # which is what makes the BGTO constraint actually buy something. (On B1I it
        # would not: reconnecting the disjoint bands creates an IFB unknown that costs
        # back exactly the clock unknown the constraint pinned.)
        beidou_l1(prns; bgto = nothing) = [
            ranging_state(
                b1c_decoder(
                    galileo[i];
                    (isnothing(bgto) ? (;) :
                     (;
                        bgtos = Dictionary(
                            [1],
                            [
                                GNSSDecoder.BeiDouB1CBGTO(;
                                    GNSS_ID = 1, WN_0BGTO = 800, t_0BGTO = 0,
                                    A_0BGTO = bgto, A_1BGTO = 0.0, A_2BGTO = 0.0),
                            ],
                        )
                    ))...,
                ),
                BeiDouB1C_D(),
                user,
                t_rx + structural + Δ,
            ) for i in eachindex(prns)
        ]
        # Three GPS plus one BeiDou cannot pay for two independent clocks (3 + 2 = 5
        # unknowns against 4 measurements). With a broadcast BGTO the BeiDou clock is
        # tied to GPS and the epoch solves anyway — the BeiDou counterpart of the
        # Galileo GGTO constraint.
        scarce = (signal_group(gps[1:3]), signal_group(beidou_l1([21]; bgto = Δ)))
        # Pin the *shape*, not just the outcome: two time systems and four satellites,
        # so the independent layout is arithmetically impossible and the BGTO
        # constraint is the only route to a fix. That is what makes this the assertion
        # that catches a structural time offset — without the constraint, a 14 s error
        # is absorbed by the BeiDou clock and only shows in the reported bias; with the
        # clock difference pinned to the steering term there is nowhere for 4.2e9 m to
        # go and the position itself is wrong. Stated here because a suite whose every
        # case has enough satellites for the independent layout can miss a defect that
        # produces no fixes at all.
        @test length(unique(get_time_system(group.signal) for group in scarce)) == 2
        @test sum(length(group.satellites) for group in scarce) == 4
        scarce_layout = PositionVelocityTime.decide_bias_layout(measurements_of(scarce))
        @test scarce_layout.clock_systems == [GPST(), BDT()]
        @test length(scarce_layout.constraints) == 1

        constrained = calc_pvt(scarce; kw...)
        @test constrained.reference_system == GPST()
        @test length(constrained.satellites) == 4
        # A 4-satellite solve has no redundancy at all, so the millimetre-scale
        # time-of-week rounding noise described above is not averaged down here the
        # way it is in the 5-satellite fixes; a few centimetres is the honest bound.
        @test norm(constrained.position - user) < 5e-2
        # The constrained system's inter-system bias is the broadcast offset rather
        # than an estimate from the geometry: four satellites leave exactly the clock
        # difference unobservable, so the fit meets the constraint exactly.
        @test ustrip(m, constrained.inter_system_biases[BDT()]) ≈ -C * Δ atol = 1e-2
        # And the solution says so.
        bgto = only(constrained.applied_broadcast_inter_system_biases)
        @test (bgto.system, bgto.target) == (BDT(), GPST())
        @test ustrip(m, bgto.value) ≈ -C * Δ atol = 1e-2
        @test bgto.value ≈ constrained.inter_system_biases[BDT()] atol = 1e-6m

        # Without the BGTO the same four satellites are unsolvable — no constraint is
        # available and the independent layout cannot be paid for — and `calc_pvt`
        # says so by returning the previous solution rather than throwing.
        @test calc_pvt((signal_group(gps[1:3]), signal_group(beidou_l1([21]))); kw...).position ==
              PVTSolution().position

        # A Galileo satellite carrying a GGTO and a BeiDou one carrying a BGTO, in the
        # same epoch: only as many offsets are used as the satellites leave clock
        # differences unobservable.
        ggto_data = GNSSDecoder.GalileoINAVData(
            galileo[5].decoder.data;
            A_0G = 3.0e-9, A_1G = 0.0, t_0G = 0,
            WN_0G = galileo[5].decoder.data.WN,
        )
        galileo_ggto = SatelliteState(;
            decoder = GNSSDecoder.GNSSDecoderState(
                galileo[5].decoder; data = ggto_data, raw_data = ggto_data),
            signal = galileo[5].signal,
            code_phase = galileo[5].code_phase,
            carrier_doppler = galileo[5].carrier_doppler,
            carrier_phase = galileo[5].carrier_phase,
        )
        # Three GPS + one BeiDou + one Galileo, all on L1. Three independent clocks
        # would need six satellites, two clocks five — so one constraint is enough, and
        # a second would only let a broadcast offset override what the geometry
        # observes. BeiDou and Galileo tie on the satellite count; BeiDou appears first.
        mixed = (
            signal_group(gps[1:3]),
            signal_group(beidou_l1([21]; bgto = Δ)),
            signal_group([galileo_ggto]),
        )
        mixed_measurements = measurements_of(mixed)
        @test map(measurement -> measurement.time_system, mixed_measurements) ==
              [GPST(), GPST(), GPST(), BDT(), GST()]
        @test all(measurement -> measurement.band_id == :L1, mixed_measurements)
        layout = PositionVelocityTime.decide_bias_layout(mixed_measurements)
        @test !isnothing(layout)
        @test layout.clock_systems == [GPST(), BDT(), GST()]
        constraint = only(layout.constraints)
        @test (constraint.column, constraint.target_column) == (2, 1)
        @test constraint.value ≈ -C * Δ atol = 1e-2
        @test layout.num_clock_groups == 2
        mixed_pvt = calc_pvt(mixed; kw...)
        @test norm(mixed_pvt.position - user) < 5e-2
        @test ustrip(m, mixed_pvt.inter_system_biases[BDT()]) ≈ -C * Δ atol = 1e-2
        @test mixed_pvt.applied_broadcast_inter_system_biases ==
              [BroadcastInterSystemBias(BDT(), GPST(), constraint.value * m)]

        # The constraint goes to the system with the fewest satellites, wherever it
        # appears: with a second Galileo satellite, BeiDou's single one is the least
        # observed clock even though Galileo comes first. (A relabelled copy of a
        # measurement stands in for the second Galileo satellite; only the decision is under
        # test.)
        with_prn(measurement, prn) = PositionVelocityTime.SatelliteMeasurement(
            (f === :prn ? prn : getfield(measurement, f)
             for f in fieldnames(typeof(measurement)))...)
        gps_sats = mixed_measurements[1:3]
        beidou_sat, galileo_sat = mixed_measurements[4], mixed_measurements[5]
        fewest = PositionVelocityTime.decide_bias_layout([gps_sats[1:2];
            [galileo_sat, with_prn(galileo_sat, galileo_sat.prn + 1), beidou_sat]])
        @test fewest.clock_systems == [GPST(), GST(), BDT()]
        @test map(c -> (c.column, c.target_column), fewest.constraints) == [(3, 1)]

        # Two GPS + one BeiDou + one Galileo: four satellites against 3 + 3 unknowns,
        # so both offsets are needed, and both are used.
        both = PositionVelocityTime.decide_bias_layout(
            [gps_sats[1:2]; [beidou_sat, galileo_sat]])
        @test map(c -> (c.column, c.target_column), both.constraints) ==
              [(2, 1), (3, 1)]
        @test both.num_clock_groups == 1
        # One offset short of that — the Galileo satellite without its GGTO — the epoch
        # cannot be solved.
        plain_galileo_sat = only(measurements_of(signal_group([galileo[5]])))
        @test isnothing(PositionVelocityTime.decide_bias_layout(
            [gps_sats[1:2]; [beidou_sat, plain_galileo_sat]]))
        # It still leaves the offsets it has, not those of the epoch before: BeiDou tied
        # to GPS, Galileo on its own clock.
        reused = PositionVelocityTime.decide_bias_layout(
            [gps_sats[1:2]; [beidou_sat, galileo_sat]])
        @test !PositionVelocityTime.decide_bias_layout!(reused,
            [gps_sats[1:2]; [beidou_sat, plain_galileo_sat]])
        @test map(c -> (c.column, c.target_column), reused.constraints) == [(2, 1)]
        @test reused.clock_group_columns == [1, 1, 1, 2]
        @test reused.num_clock_groups == 2
        # With no offset to tie with, none.
        @test !PositionVelocityTime.decide_bias_layout!(reused,
            [gps_sats[1:2]; [plain_galileo_sat]])
        @test isempty(reused.constraints)
        @test reused.num_clock_groups == 2
        # One satellite per system, both offsets available: unsolvable even with every
        # clock tied, but the layout ties all three, through both offsets.
        @test !PositionVelocityTime.decide_bias_layout!(reused,
            [gps_sats[1:1]; [beidou_sat, galileo_sat]])
        @test sort(map(c -> (c.column, c.target_column), reused.constraints)) ==
              [(2, 1), (3, 1)]
        @test reused.clock_group_columns == [1, 1, 1]
        @test reused.num_clock_groups == 1

        # A constraint is not GPS-specific. A GPS-free scarce epoch ties BeiDou to
        # Galileo instead: the BGTO's GNSS_ID names Galileo (2) as well as GPS (1)
        # — BDS-SIS-ICD-B1C-1.0 §7.13.1. (The leap-second offset of BDT against GST is
        # the same 14 s as against GPST — both GPST and GST are TAI − 19 s — so the same
        # `structural` and steering values apply.)
        beidou_gal_bgto = ranging_state(
            b1c_decoder(
                galileo[4];
                bgtos = Dictionary(
                    [2],
                    [
                        GNSSDecoder.BeiDouB1CBGTO(;
                            GNSS_ID = 2, WN_0BGTO = 800, t_0BGTO = 0,
                            A_0BGTO = Δ, A_1BGTO = 0.0, A_2BGTO = 0.0),
                    ],
                ),
            ),
            BeiDouB1C_D(),
            user,
            t_rx + structural + Δ,
        )
        # Re-synthesized against this testset's own `user`/`t_rx`, like every other
        # satellite here; the top-level `galileo` states were built for another epoch.
        gal_states = [
            ranging_state(st.decoder, st.signal, user, t_rx) for st in galileo[1:3]
        ]
        gps_free = (signal_group(gal_states), signal_group([beidou_gal_bgto]))
        gst_layout = PositionVelocityTime.decide_bias_layout(measurements_of(gps_free))
        @test !isnothing(gst_layout)
        @test gst_layout.clock_systems == [GST(), BDT()]
        gst_constraint = only(gst_layout.constraints)
        @test (gst_constraint.column, gst_constraint.target_column) == (2, 1)

        gst_constrained = calc_pvt(gps_free; kw...)
        @test gst_constrained.reference_system == GST()
        @test length(gst_constrained.satellites) == 4
        @test norm(gst_constrained.position - user) < 5e-2
        @test ustrip(m, gst_constrained.inter_system_biases[BDT()]) ≈ -C * Δ atol = 1e-2
    end
end
