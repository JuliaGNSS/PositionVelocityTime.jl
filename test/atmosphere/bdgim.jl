@testset "BDGIM ionospheric model (BeiDou)" begin
    PVT = PositionVelocityTime
    # A plausible mid-latitude daytime coefficient set, in the ICD's TECu, and one
    # that could actually be broadcast: Table 7-10 makes α_1, α_3 and α_4 unsigned,
    # α_5 unsigned with a *negative* scale factor (so α_5 ≤ 0), and only α_2 and
    # α_6…α_9 two's complement. α_1 weights the degree-0 harmonic (P̃₀,₀ ≡ 1), so it
    # carries the bulk of the vertical TEC and the remaining eight shape it — the
    # predictive part A_0 spans only degrees 3-5 and is a near-zero-mean structure
    # term (±10-15 TECu globally), not a level.
    α = (35.0, -1.5, 2.25, 0.75, -0.5, 0.25, -0.125, 0.375, -0.625)
    lla = LLA(48.0, 11.0, 550.0)
    week, tow = 911, 302400.0   # BDT week 911, Wednesday noon → 2023-06-21T12:00
    mjd = PVT._bdgim_modified_julian_date(week, tow)

    # -- Transcription safety -------------------------------------------------
    #
    # BDS-SIS-ICD-B1C-1.0 Table 7-12 is 425 hand-copied numbers, and a single wrong
    # digit would shift TEC without failing anything else here. These checksums pin
    # the table as a whole: every entry contributes to its row sum, its column sum
    # and the absolute total, so no single-entry typo — including a dropped minus
    # sign, which the absolute total does not see but the signed sums do — survives
    # re-transcription. The expected values were computed from the ICD text layer
    # (cross-read against the rendered PDF page and against the same table restated
    # as BDS-SIS-ICD-B2b-1.0 Table 7-10), not from this package's own constants.
    @testset "Table 7-12 non-broadcast coefficients: shape and checksums" begin
        A0, A, B = PVT._BDGIM_A0, PVT._BDGIM_A, PVT._BDGIM_B
        @test length(A0) == 17
        @test length(A) == length(B) == 12          # k = 1…12
        @test all(==(17), length.(A)) && all(==(17), length.(B))   # j = 1…17
        @test length(PVT._BDGIM_PERIODS) == 12
        @test length(PVT._BDGIM_NM_PREDICTED) == 17
        @test length(PVT._BDGIM_NM) == 9

        # Constant row a_{0,j}
        @test sum(A0) ≈ -1.17 atol = 1e-9
        # Row sums over j, one per k (a_{k,j} then b_{k,j})
        row_sums_a = [
            -0.07, 0.27, -0.19, -0.01, 0.05, -0.25,
            0.17, 2.38, 0.36, 0.19, -0.04, 0.01,
        ]
        row_sums_b = [
            1.4, -0.2, -0.19, -0.03, 0.02, 0.09,
            -0.31, 1.39, 0.26, 0.13, 0.08, 0.05,
        ]
        @test collect(sum.(A)) ≈ row_sums_a atol = 1e-9
        @test collect(sum.(B)) ≈ row_sums_b atol = 1e-9
        # Column sums over k, one per j
        col_sums_a = [
            -0.41, -0.42, -1.86, 0.73, 1.06, -0.55, -0.13, 3.15, 0.24,
            -0.22, -0.38, 0.12, 0.22, 0.36, 0.97, 0.09, -0.1,
        ]
        col_sums_b = [
            0.94, -0.76, -0.76, 0.5, 0.21, -0.14, 0.13, 0.61, 0.57,
            0.57, 0.22, -0.26, -0.44, 0.1, 0.51, 0.19, 0.5,
        ]
        @test [sum(A[k][j] for k = 1:12) for j = 1:17] ≈ col_sums_a atol = 1e-9
        @test [sum(B[k][j] for k = 1:12) for j = 1:17] ≈ col_sums_b atol = 1e-9
        # Magnitude total and sparsity: catches a sign flip paired with a value typo
        # that happened to leave a signed sum intact.
        all_entries =
            vcat(collect(A0), reduce(vcat, collect.(A)), reduce(vcat, collect.(B)))
        @test length(all_entries) == 17 * 25
        @test sum(abs, all_entries) ≈ 44.79 atol = 1e-9
        @test count(!iszero, all_entries) == 318

        # Spot values at the table's corners and at a few of its largest entries.
        @test A0[1] == -0.61 && A0[17] == -0.12 && A0[16] == 0.0
        @test A[1][1] == -0.51 && A[12][17] == 0.0
        @test B[1][1] == 0.23 && B[12][8] == 0.07
        @test A[9][8] == 1.58          # k = 9 (T = 4028.71 d), j = 8 (n/m = 4/0)
        @test A[8][1] == 1.09          # k = 8 (T = 365.25 d), j = 1 (n/m = 3/0)
        @test A[9][3] == -1.22

        # Periods T_k [days] and the (n_j, m_j) header row, both from Table 7-12.
        @test PVT._BDGIM_PERIODS == (
            1.0, 0.5, 0.33, 14.6, 27.0, 121.6, 182.51, 365.25, 4028.71, 2014.35,
            1342.90, 1007.18,
        )
        @test PVT._BDGIM_NM_PREDICTED[1] == (3, 0)
        @test PVT._BDGIM_NM_PREDICTED[8] == (4, 0)
        @test PVT._BDGIM_NM_PREDICTED[17] == (5, -2)
        # Table 7-11, the nine broadcast harmonics.
        @test PVT._BDGIM_NM ==
              ((0, 0), (1, 0), (1, 1), (1, -1), (2, 0), (2, 1), (2, -1), (2, 2), (2, -2))
        # Suggested parameter values (§7.8.2, after Eq. 7-17). The pole longitude is
        # *negative* — the ICD text layer drops the minus sign, the PDF does not.
        @test PVT._BDGIM_RE == 6378.0
        @test PVT._BDGIM_HION == 400.0
        @test PVT._BDGIM_POLE_LAT ≈ deg2rad(80.27)
        @test PVT._BDGIM_POLE_LON ≈ deg2rad(-72.58)
    end

    @testset "Legendre functions (Eq. 7-12, 7-13)" begin
        # Against the closed forms, not against a reimplementation of the recursion.
        x = 0.3
        @test PVT._bdgim_legendre(0, 0, x) == 1.0
        @test PVT._bdgim_legendre(1, 0, x) ≈ x
        @test PVT._bdgim_legendre(1, 1, x) ≈ sqrt(1 - x^2)
        @test PVT._bdgim_legendre(2, 0, x) ≈ (3x^2 - 1) / 2
        @test PVT._bdgim_legendre(2, 1, x) ≈ 3x * sqrt(1 - x^2)
        @test PVT._bdgim_legendre(2, 2, x) ≈ 3 * (1 - x^2)
        @test PVT._bdgim_legendre(3, 0, x) ≈ (5x^3 - 3x) / 2
        @test PVT._bdgim_legendre(3, 3, -0.4) ≈ 15 * (1 - 0.16)^1.5
        @test PVT._bdgim_legendre(4, 2, 0.7) ≈ 7.5 * (7 * 0.49 - 1) * (1 - 0.49)
        # Normalisation N_{n,m} = √[(n−m)!(2n+1)(2−δ₀ₘ)/(n+m)!]
        @test PVT._bdgim_normalization(0, 0) == 1.0
        @test PVT._bdgim_normalization(1, 0) ≈ sqrt(3)
        @test PVT._bdgim_normalization(1, 1) ≈ sqrt(3)          # 0!·3·2/2!
        @test PVT._bdgim_normalization(2, 2) ≈ sqrt(10 / 24)
        @test PVT._bdgim_normalization(5, 0) ≈ sqrt(11)
    end

    @testset "reference frames (Eq. 7-8, 7-9, 7-10)" begin
        # At the zenith the Earth-central angle ψ vanishes, so the pierce point is
        # the user's own meridian and parallel.
        φ_g, λ_g = PVT._bdgim_pierce_point(deg2rad(48), deg2rad(11), π / 2, 0.0)
        @test φ_g ≈ deg2rad(48) atol = 1e-12
        @test λ_g ≈ deg2rad(11) atol = 1e-12
        # ψ at the horizon: π/2 − arcsin(Re/(Re+H_ion)) ≈ 0.3078 rad ≈ 17.6°, so the
        # pierce point due north of a user at 48°N sits at ~65.6°N.
        φ_north, λ_north = PVT._bdgim_pierce_point(deg2rad(48), deg2rad(11), 0.0, 0.0)
        ψ = π / 2 - asin(6378.0 / 6778.0)
        @test φ_north ≈ deg2rad(48) + ψ atol = 1e-12
        @test λ_north ≈ deg2rad(11) atol = 1e-12
        # The magnetic pole itself maps to geomagnetic latitude 90°.
        φ_m, _ = PVT._bdgim_solar_geomagnetic(
            PVT._BDGIM_POLE_LAT,
            PVT._BDGIM_POLE_LON,
            mjd,
        )
        @test φ_m ≈ π / 2 atol = 1e-6
        # S_lon = π(1 − 2·frac(MJD)) puts the sun's mean longitude at Greenwich at
        # 12:00 and at the date line at 00:00, so the solar-fixed longitude λ′ of a
        # point on the sub-solar meridian is 0 and of its antipode is ±π.
        @test last(PVT._bdgim_solar_geomagnetic(0.0, 0.0, 60116.5)) ≈ 0.0 atol = 1e-12
        @test abs(last(PVT._bdgim_solar_geomagnetic(0.0, π, 60116.5))) ≈ π atol = 1e-12
        @test last(PVT._bdgim_solar_geomagnetic(0.0, π, 60116.0)) ≈ 0.0 atol = 1e-12
    end

    @testset "BDT week/TOW → Modified Julian Date" begin
        # BDT week 0 starts 2006-01-01, MJD 53736.
        @test PVT._bdgim_modified_julian_date(0, 0.0) == 53736.0
        @test PVT._bdgim_modified_julian_date(0, 43200.0) == 53736.5
        @test PVT._bdgim_modified_julian_date(1, 0.0) == 53743.0
        # Cross-checked against Dates: MJD is JD − 2400000.5.
        @test mjd ≈ datetime2julian(DateTime(2023, 6, 21, 12)) - 2400000.5 atol = 1e-9
    end

    @testset "prediction reference epoch t_p (Eq. 7-15)" begin
        day = 60116.0
        # The nearest odd hour of the same day: 00:00–02:00 → 01:00, and so on.
        # Compared in hours past midnight, so the tolerance is not swamped by the
        # ~60000-day magnitude of the MJD itself.
        hour_of(mjd) = (PVT._bdgim_reference_epoch(mjd) - day) * 24
        @test hour_of(day) ≈ 1.0
        @test hour_of(day + 1.9 / 24) ≈ 1.0
        @test hour_of(day + 2.0 / 24) ≈ 3.0
        @test hour_of(day + 12.5 / 24) ≈ 13.0
        @test hour_of(day + 23.99 / 24) ≈ 23.0
        # Always on the same day, and always within an hour of the epoch.
        for frac in range(0, 1, length = 97)[1:(end-1)]
            t_p = PVT._bdgim_reference_epoch(day + frac)
            @test floor(t_p) == day
            @test abs(t_p - (day + frac)) ≤ 1 / 24 + 1e-12
        end
    end

    @testset "prediction part is continuous across its period boundaries" begin
        # β_j (Eq. 7-15) is a Fourier sum in the *absolute* MJD t_p, so nothing wraps:
        # crossing an MJD midnight, a whole T₁ = 1 day period, or the shortest
        # T₃ = 0.33 day period must not step. (This is what makes the two-hourly t_p
        # a quantisation of a smooth function rather than a model discontinuity.)
        for boundary in (60116.0, 60117.0, 60116.0 + 0.33, 60116.0 + 0.5)
            before = PVT._bdgim_prediction_amplitudes(boundary - 1e-7)
            after = PVT._bdgim_prediction_amplitudes(boundary + 1e-7)
            @test maximum(abs.(before .- after)) < 1e-5
        end
        # What *is* stepped is the argument: t_p is quantised to the nearest odd hour,
        # so the amplitudes are piecewise constant in the epoch, changing every two
        # hours and only by as much as the smooth series moved over those two hours.
        # (VTEC itself still varies continuously in between — the solar-fixed
        # longitude λ′ tracks the true epoch, not t_p.)
        β_of(t) = PVT._bdgim_prediction_amplitudes(PVT._bdgim_reference_epoch(t))
        block = [β_of(60116.0 + (0.1 + h) / 24) for h in (0.0, 0.5, 1.0, 1.8)]
        @test all(==(block[1]), block)                    # constant within a block
        jump = maximum(abs.(β_of(60116.0 + 2.05 / 24) .- β_of(60116.0 + 1.95 / 24)))
        # Two hours of the fastest term (T₃ = 0.33 d ≈ 7.9 h) is the scale here: a
        # fraction of a TECu of amplitude, well under the amplitudes themselves
        # (max |β_j| ≈ 2.7 TECu at this epoch).
        @test 0 < jump < 1.0
    end

    @testset "mapping function (Eq. 7-17)" begin
        @test PVT.bdgim_mapping_function(π / 2) == 1.0            # exactly 1 at zenith
        for el in deg2rad.((0, 5, 10, 30, 60, 89))
            @test PVT.bdgim_mapping_function(el) ≥ 1.0
        end
        # Monotonically increasing towards the horizon, to ~2.95 at 0° elevation.
        mfs = PVT.bdgim_mapping_function.(deg2rad.((90, 60, 30, 10, 5, 0)))
        @test issorted(mfs)
        @test mfs[end] ≈ 1 / sqrt(1 - (6378 / 6778)^2)
        @test mfs[end] ≈ 2.9547 atol = 1e-4
    end

    @testset "slant ≥ vertical, with equality at the zenith" begin
        # At the zenith the pierce point is the user's own position and M_F is 1, so
        # slant TEC and vertical TEC coincide exactly.
        vtec_zenith = PVT.bdgim_vtec(deg2rad(lla.lat), deg2rad(lla.lon), mjd, α)
        @test PVT.bdgim_stec(π / 2, 0.0, lla, mjd, α) ≈ vtec_zenith rtol = 1e-12
        for el in deg2rad.((5, 15, 30, 60, 85)), az in deg2rad.((0, 90, 200, 315))
            φ_g, λ_g = PVT._bdgim_pierce_point(
                deg2rad(lla.lat),
                deg2rad(lla.lon),
                el,
                az,
            )
            vtec = PVT.bdgim_vtec(φ_g, λ_g, mjd, α)
            stec = PVT.bdgim_stec(el, az, lla, mjd, α)
            @test stec ≥ vtec                                   # M_F ≥ 1
            @test stec ≈ PVT.bdgim_mapping_function(el) * vtec rtol = 1e-12
        end
    end

    @testset "non-positive VTEC is floored at zero" begin
        # The predictive part A_0 spans degrees 3-5 only, so it has no level of its
        # own and is negative over much of the globe (here ≈ −9.7 TECu). A broadcast
        # set too small to lift it — an all-zero α is a legal message, and a genuinely
        # quiet ionosphere broadcasts a small α₁ — therefore reaches Eq. 7-16 with a
        # negative sum. The ICD prescribes no handling of that; this package floors it
        # at zero, because a negative VTEC would enter as a group *advance*.
        zeros9 = ntuple(_ -> 0.0, 9)
        @test PVT.bdgim_vtec(deg2rad(48), deg2rad(11), mjd, zeros9) == 0.0
        @test PVT.bdgim_stec(deg2rad(30), 0.0, lla, mjd, zeros9) == 0.0
        @test PVT.ionospheric_delay(
            PVT.BDGIMParams(zeros9..., week),
            carrier_hz(BeiDouB1C_D()),
            deg2rad(30),
            0.0,
            lla,
            tow,
        ) == 0.0
        # Still zero, not negative, when a small positive α₁ does not quite lift it.
        small = (5.0, ntuple(_ -> 0.0, 8)...)
        @test PVT.bdgim_vtec(deg2rad(48), deg2rad(11), mjd, small) == 0.0
        # The clamp is a floor, not a blanket: an ordinary set is unaffected.
        @test PVT.bdgim_vtec(deg2rad(48), deg2rad(11), mjd, α) > 0.0
    end

    @testset "delay scales exactly as 1/f² across signals" begin
        p = PVT.BDGIMParams(α..., week)
        geometry = (deg2rad(35), deg2rad(120), lla, tow)
        d_b1c = PVT.ionospheric_delay(p, carrier_hz(BeiDouB1C_D()), geometry...)
        d_b2a = PVT.ionospheric_delay(p, carrier_hz(BeiDouB2aI()), geometry...)
        d_b2b = PVT.ionospheric_delay(p, carrier_hz(BeiDouB2bI()), geometry...)
        f(sig) = GNSSSignals.get_center_frequency(sig)
        # No reference frequency: each satellite's own carrier enters Eq. 7-6, so the
        # ratio of any two delays is exactly the inverse square of their frequencies.
        @test d_b2a / d_b1c ≈ (f(BeiDouB1C_D()) / f(BeiDouB2aI()))^2 rtol = 1e-12
        @test d_b2b / d_b1c ≈ (f(BeiDouB1C_D()) / f(BeiDouB2bI()))^2 rtol = 1e-12
        @test d_b2a > d_b1c   # B2a is the lower carrier, so the larger delay
        # The same set applied to a non-BeiDou satellite uses that satellite's carrier
        # too — the model is TEC, not a delay at a reference band.
        d_l5 = PVT.ionospheric_delay(p, carrier_hz(GPSL5I()), geometry...)
        @test d_l5 / d_b1c ≈ (f(BeiDouB1C_D()) / f(GPSL5I()))^2 rtol = 1e-12
    end

    @testset "magnitude sanity for a plausible daytime coefficient set" begin
        p = PVT.BDGIMParams(α..., week)
        # ~26 TECu of vertical TEC (α₁ = 35 minus the negative A_0 here) is ~4 m at
        # B1C; a low-elevation ray through the same ionosphere is a few metres more.
        # The family stays in the decimetres-to-metres band a single-frequency L-band
        # correction should, and never in the tens of metres.
        d_zenith = PVT.ionospheric_delay(p, carrier_hz(BeiDouB1C_D()), π / 2, 0.0, lla, tow)
        @test 1.0 < d_zenith < 8.0
        d_low = PVT.ionospheric_delay(p, carrier_hz(BeiDouB1C_D()), deg2rad(10), 0.0, lla, tow)
        @test d_low > d_zenith
        @test d_low < 25.0
        # A quieter ionosphere — the same shape with α₁ down to 12 TECu — lands in the
        # decimetres instead.
        quiet = PVT.BDGIMParams(12.0, α[2:9]..., week)
        @test 0.1 < PVT.ionospheric_delay(quiet, carrier_hz(BeiDouB1C_D()), π / 2, 0.0, lla, tow) < 1.5
    end

    @testset "regression pins" begin
        # Computed by this implementation, not ICD truth — BDS-SIS-ICD-B1C-1.0 ships
        # no worked example for BDGIM. They exist so a refactor that changes a number
        # has to say so out loud.
        vtec = PVT.bdgim_vtec(deg2rad(48), deg2rad(11), mjd, α)
        @test vtec ≈ 25.9776383593909 rtol = 1e-10
        stec = PVT.bdgim_stec(deg2rad(20), deg2rad(135), lla, mjd, α)
        @test stec ≈ 60.5472505715464 rtol = 1e-10
        @test PVT.ionospheric_delay(
            PVT.BDGIMParams(α..., week),
            carrier_hz(BeiDouB1C_D()),
            deg2rad(20),
            deg2rad(135),
            lla,
            tow,
        ) ≈ 9.8263262553667 rtol = 1e-10
        # The prediction path on its own: the Fourier amplitudes β_j of Eq. 7-15 at a
        # fixed t_p, which is where the whole of Table 7-12 lands. Pinned unfloored,
        # unlike the VTEC above.
        β = PVT._bdgim_prediction_amplitudes(60116.0 + 13 / 24)
        @test β[1] ≈ -1.5706740201187 rtol = 1e-10
        @test β[8] ≈ 2.2072433666267 rtol = 1e-10
        @test β[17] ≈ -0.5520394639058 rtol = 1e-10
    end

    # -- Decoder extraction and selection -------------------------------------

    bdgim_fields(α) = (;
        α_bdgim_1 = α[1], α_bdgim_2 = α[2], α_bdgim_3 = α[3], α_bdgim_4 = α[4],
        α_bdgim_5 = α[5], α_bdgim_6 = α[6], α_bdgim_7 = α[7], α_bdgim_8 = α[8],
        α_bdgim_9 = α[9],
    )
    function bdgim_decoder(kind, α, WN)
        state, data = if kind === :b1c
            GNSSDecoder.BeiDouB1CDecoderState(20),
            GNSSDecoder.BeiDouB1CData(; WN, bdgim_fields(α)...)
        elseif kind === :b2a
            GNSSDecoder.BeiDouB2aDecoderState(20),
            GNSSDecoder.BeiDouB2aData(; WN, bdgim_fields(α)...)
        else
            GNSSDecoder.BeiDouB2bDecoderState(20),
            GNSSDecoder.BeiDouB2bData(; WN, bdgim_fields(α)...)
        end
        GNSSDecoder.GNSSDecoderState(state; data, raw_data = data)
    end

    @testset "parameter extraction from B-CNAV decoders" begin
        for kind in (:b1c, :b2a, :b2b)
            p = PositionVelocityTime.bdgim_params(bdgim_decoder(kind, α, Int64(911)))
            @test p isa PositionVelocityTime.BDGIMParams
            @test (p.α_1, p.α_2, p.α_3, p.α_4, p.α_5, p.α_6, p.α_7, p.α_8, p.α_9) == α
            @test p.week_number == 911
            # A fresh decoder of the same kind has nothing yet.
            @test PositionVelocityTime.bdgim_params(
                bdgim_decoder(kind, ntuple(_ -> nothing, 9), nothing),
            ) === nothing
            # All nine coefficients *and* the week number are required: the week is
            # decoded in a different subframe/message type, so half-populated is the
            # ordinary transient state.
            @test PositionVelocityTime.bdgim_params(bdgim_decoder(kind, α, nothing)) ===
                  nothing
            partial = ntuple(i -> i == 5 ? nothing : α[i], 9)
            @test PositionVelocityTime.bdgim_params(
                bdgim_decoder(kind, partial, Int64(911)),
            ) === nothing
        end
        # Constellations and message families that broadcast no BDGIM set fall to the
        # generic method rather than dispatching on fields they do not have.
        @test PositionVelocityTime.bdgim_params(GNSSDecoderState(GPSL1CA(), 1)) === nothing
        @test PositionVelocityTime.bdgim_params(GNSSDecoderState(GalileoE1B(), 1)) ===
              nothing
        @test PositionVelocityTime.bdgim_params(GNSSDecoderState(BeiDouB1I(), 1)) ===
              nothing
        # ... and a B-CNAV decoder yields no Klobuchar or NTCM-G set.
        b1c = bdgim_decoder(:b1c, α, Int64(911))
        @test PositionVelocityTime.klobuchar_params(b1c) === nothing
        @test PositionVelocityTime.ntcm_g_params(b1c) === nothing
    end

    @testset "selection order: NTCM-G > BDGIM > GPS Klobuchar > BeiDou Klobuchar" begin
        mk(dec, sys) = SatelliteState(;
            decoder = dec,
            signal = sys,
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        klob_α = (4.6566129e-09, 1.4901161e-08, -5.96046e-08, -5.96046e-08)
        klob_β = (79872.0, 65536.0, -65536.0, -393216.0)
        gps = mk(
            GNSSDecoder.GNSSDecoderState(
                GNSSDecoderState(GPSL1CA(), 1);
                data = GNSSDecoder.GPSL1CAData(
                    GNSSDecoderState(GPSL1CA(), 1).data;
                    α_0 = klob_α[1], α_1 = klob_α[2], α_2 = klob_α[3], α_3 = klob_α[4],
                    β_0 = klob_β[1], β_1 = klob_β[2], β_2 = klob_β[3], β_3 = klob_β[4],
                ),
            ),
            GPSL1CA(),
        )
        gal = mk(
            GNSSDecoder.GNSSDecoderState(
                GNSSDecoderState(GalileoE1B(), 1);
                data = GNSSDecoder.GalileoINAVData(
                    GNSSDecoderState(GalileoE1B(), 1).data;
                    a_i0 = 121.13, a_i1 = 0.35, a_i2 = 0.013, WN = 1100,
                ),
            ),
            GalileoE1B(),
        )
        bds3 = mk(bdgim_decoder(:b2a, α, Int64(911)), BeiDouB2aI())
        dnav_data = GNSSDecoder.BeiDouDNAVData(;
            α_0 = 1.0e-8, α_1 = 2.0e-8, α_2 = -3.0e-8, α_3 = 4.0e-8,
            β_0 = 90112.0, β_1 = 16384.0, β_2 = -196608.0, β_3 = 131072.0,
        )
        bds2 = mk(
            GNSSDecoder.GNSSDecoderState(
                GNSSDecoder.BeiDouB1IDecoderState(6);
                data = dnav_data,
                raw_data = dnav_data,
            ),
            BeiDouB1I(),
        )
        # A second legacy satellite, needed below to place one *after* the BDS-3 one.
        # Its own PRN: a signal group is keyed by PRN, so the same satellite twice on
        # one signal is rejected structurally rather than silently counted twice.
        bds2_second = mk(
            GNSSDecoder.GNSSDecoderState(
                GNSSDecoder.BeiDouB1IDecoderState(7);
                data = dnav_data,
                raw_data = dnav_data,
            ),
            BeiDouB1I(),
        )

        select = PositionVelocityTime.select_ionospheric_parameters
        # Each of these satellites is on a signal of its own, so each is its own group;
        # the tuple's order is the order the scan sees them in, which is what the
        # preference rungs below are asserted against.
        one(state) = signal_group([state])
        # BDGIM alone corrects a BDS-3-only epoch, which used to get nothing at all.
        @test select(one(bds3)) isa PositionVelocityTime.BDGIMParams
        # It beats both Klobuchar sources, and loses to NTCM-G.
        @test select((one(bds3), one(gps))) isa PositionVelocityTime.BDGIMParams
        @test select((one(gps), one(bds3))) isa PositionVelocityTime.BDGIMParams
        @test select((one(bds3), one(bds2))) isa PositionVelocityTime.BDGIMParams
        @test select((one(bds2), one(bds3))) isa PositionVelocityTime.BDGIMParams
        @test select((one(bds3), one(gal))) isa PositionVelocityTime.NTCMGParams
        @test select((one(gal), one(bds3))) isa PositionVelocityTime.NTCMGParams
        @test select((one(gal), one(bds3), one(gps), one(bds2))) isa
              PositionVelocityTime.NTCMGParams
        # Without a BDS-3 satellite the previous order is untouched, and each
        # Klobuchar source keeps its own variant of the model.
        @test select((one(gps), one(bds2))) isa PositionVelocityTime.KlobucharParams
        @test select(one(bds2)) isa PositionVelocityTime.BeiDouKlobucharParams
        # The BeiDou branch asks both accessors: a legacy satellite alongside a BDS-3
        # one contributes its Klobuchar set without shadowing the BDGIM one. The two
        # legacy satellites are deliberately in separate groups, which is what puts one
        # of them *after* the BDS-3 satellite in scan order.
        @test select((one(bds2), one(bds3), one(bds2_second))) isa
              PositionVelocityTime.BDGIMParams
        # An epoch with no groups at all — the receiver tracking nothing yet.
        @test select(()) === nothing
    end
end
