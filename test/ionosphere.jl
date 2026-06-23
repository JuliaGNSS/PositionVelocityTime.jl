# Independent RTKLIB-style reference implementation (inputs in radians) used to
# cross-check the package's semicircle-based `klobuchar_group_delay`.
function klobuchar_reference(lat, lon, az, el, gps_time, α, β)
    ψ = 0.0137 / (el / π + 0.11) - 0.022
    φ = clamp(lat / π + ψ * cos(az), -0.416, 0.416)
    λ = lon / π + ψ * sin(az) / cos(φ * π)
    φ += 0.064 * cos((λ - 1.617) * π)
    tt = 43200.0 * λ + gps_time
    tt -= floor(tt / 86400.0) * 86400.0
    f = 1.0 + 16.0 * (0.53 - el / π)^3
    amp = max(α[1] + φ * (α[2] + φ * (α[3] + φ * α[4])), 0.0)
    per = max(β[1] + φ * (β[2] + φ * (β[3] + φ * β[4])), 72000.0)
    x = 2π * (tt - 50400.0) / per
    return f * (abs(x) < 1.57 ? 5e-9 + amp * (1 - x^2 / 2 + x^4 / 24) : 5e-9)
end

@testset "Klobuchar ionospheric model" begin
    # Exemplary Klobuchar coefficients (IS-GPS-200 SI units)
    α = (4.6566129e-09, 1.4901161e-08, -5.96046e-08, -5.96046e-08)
    β = (79872.0, 65536.0, -65536.0, -393216.0)

    @testset "klobuchar_group_delay matches independent reference" begin
        for lat in deg2rad.((10.0, 48.0, -33.0)),
            lon in deg2rad.((-120.0, 11.0, 150.0)),
            az in deg2rad.((0.0, 95.0, 270.0)),
            el in deg2rad.((5.0, 25.0, 89.0)),
            t in (0.0, 50400.0, 80000.0)

            got = PositionVelocityTime.klobuchar_group_delay(
                lat / π,
                lon / π,
                el / π,
                az / π,
                t,
                α,
                β,
            )
            ref = klobuchar_reference(lat, lon, az, el, t, α, β)
            @test got ≈ ref rtol = 1e-12
        end
    end

    @testset "physical behaviour" begin
        latsc, lonsc, t = 48 / 180, 11 / 180, 50400.0  # local ~14:00 → near diurnal peak
        d90 =
            PositionVelocityTime.klobuchar_group_delay(latsc, lonsc, 90 / 180, 0.0, t, α, β)
        d30 =
            PositionVelocityTime.klobuchar_group_delay(latsc, lonsc, 30 / 180, 0.0, t, α, β)
        d10 =
            PositionVelocityTime.klobuchar_group_delay(latsc, lonsc, 10 / 180, 0.0, t, α, β)
        @test d90 > 0
        @test d30 > d90            # obliquity increases slant delay at low elevation
        @test d10 > d30
        @test 299792458.0 * d90 < 15  # zenith L1 delay is a few metres, not tens
    end

    # --- helpers: build decoders carrying broadcast coefficients ---
    function gps_decoder_with(α, β)
        dec = GNSSDecoderState(GPSL1CA(), 1)
        GNSSDecoder.GNSSDecoderState(
            dec;
            data = GNSSDecoder.GPSL1CAData(
                dec.data;
                α_0 = α[1],
                α_1 = α[2],
                α_2 = α[3],
                α_3 = α[4],
                β_0 = β[1],
                β_1 = β[2],
                β_2 = β[3],
                β_3 = β[4],
            ),
        )
    end
    function galileo_decoder_with(a_i0, a_i1, a_i2, WN)
        dec = GNSSDecoderState(GalileoE1B(), 1)
        GNSSDecoder.GNSSDecoderState(
            dec;
            data = GNSSDecoder.GalileoE1BData(dec.data; a_i0, a_i1, a_i2, WN),
        )
    end
    mkstate(dec, sys) = SatelliteState(;
        decoder = dec,
        system = sys,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )

    @testset "parameter extraction from decoders" begin
        # Fresh decoders carry no coefficients
        @test PositionVelocityTime.klobuchar_params(GNSSDecoderState(GPSL1CA(), 1)) ===
              nothing
        @test PositionVelocityTime.ntcm_g_params(GNSSDecoderState(GalileoE1B(), 1)) ===
              nothing
        # Broadcast α/β are extracted into KlobucharParams
        kp = PositionVelocityTime.klobuchar_params(gps_decoder_with(α, β))
        @test kp isa PositionVelocityTime.KlobucharParams
        @test (kp.α_0, kp.α_1, kp.α_2, kp.α_3) == α
        @test (kp.β_0, kp.β_1, kp.β_2, kp.β_3) == β
        # Broadcast a_i0..a_i2 + WN are extracted into NTCMGParams
        np = PositionVelocityTime.ntcm_g_params(
            galileo_decoder_with(236.831641, -0.39362878, 0.00402826613, 1100),
        )
        @test np isa PositionVelocityTime.NTCMGParams
        @test np.a_i0 == 236.831641
        @test np.week_number == 1100
        # Wrong-system decoders never yield the other model's params
        @test PositionVelocityTime.ntcm_g_params(gps_decoder_with(α, β)) === nothing
        @test PositionVelocityTime.klobuchar_params(
            galileo_decoder_with(100.0, 0.0, 0.0, 1100),
        ) === nothing
    end

    @testset "constellation-wide model selection" begin
        gps_bare = mkstate(GNSSDecoderState(GPSL1CA(), 1), GPSL1CA())
        gps_klob = mkstate(gps_decoder_with(α, β), GPSL1CA())
        gal_bare = mkstate(GNSSDecoderState(GalileoE1B(), 1), GalileoE1B())
        gal_ntcm = mkstate(galileo_decoder_with(121.13, 0.35, 0.013, 1100), GalileoE1B())

        # Nothing decoded → no correction
        @test PositionVelocityTime.select_ionospheric_correction([gps_bare, gal_bare]) ===
              nothing
        # Only Klobuchar → Klobuchar
        @test PositionVelocityTime.select_ionospheric_correction([gps_klob, gps_bare]) isa
              PositionVelocityTime.KlobucharParams
        # Only Galileo → NTCM-G
        @test PositionVelocityTime.select_ionospheric_correction([gal_ntcm, gal_bare]) isa
              PositionVelocityTime.NTCMGParams
        # Both available → NTCM-G wins (more accurate)
        @test PositionVelocityTime.select_ionospheric_correction([gps_klob, gal_ntcm]) isa
              PositionVelocityTime.NTCMGParams
    end

    @testset "selected model applied to every satellite" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        sat = ECEF(user[1] + 1.0e7, user[2] + 1.0e7, user[3] + 2.0e7)
        klob = PositionVelocityTime.KlobucharParams(α..., β...)
        # Line-of-sight geometry (same user/sat for every call below)
        el, az = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat)
        lla = LLAfromECEF(wgs84)(user)

        # No correction → exactly zero for any system
        @test PositionVelocityTime.ionospheric_delay(
            nothing,
            GPSL1CA(),
            el,
            az,
            lla,
            50400.0,
        ) == 0.0
        @test PositionVelocityTime.ionospheric_delay(
            nothing,
            GalileoE1B(),
            el,
            az,
            lla,
            50400.0,
        ) == 0.0
        # Klobuchar applied to GPS *and* Galileo (E1 shares the L1 frequency, so
        # the delay is identical for both systems)
        d_gps =
            PositionVelocityTime.ionospheric_delay(klob, GPSL1CA(), el, az, lla, 50400.0)
        d_gal =
            PositionVelocityTime.ionospheric_delay(klob, GalileoE1B(), el, az, lla, 50400.0)
        @test d_gps > 0.0
        @test d_gps ≈ d_gal rtol = 1e-12

        # The delay scales as 1/f², so the same coefficients applied on a lower band
        # (here GPS L5, 1176.45 MHz) give the correct larger delay, not the L1 value.
        d_l5 = PositionVelocityTime.ionospheric_delay(klob, GPSL5I(), el, az, lla, 50400.0)
        # ratio of two Hz quantities is dimensionless
        f_ratio =
            GNSSSignals.get_center_frequency(GPSL1CA()) /
            GNSSSignals.get_center_frequency(GPSL5I())
        @test d_l5 ≈ d_gps * f_ratio^2 rtol = 1e-12
        @test d_l5 > d_gps
    end
end

@testset "NTCM-G ionospheric model (Galileo)" begin
    toecef(lat, lon, h) = ECEFfromLLA(wgs84)(LLA(lat, lon, h))
    HI = (236.831641, -0.39362878, 0.00402826613)        # high solar activity
    MED = (121.129893, 0.351254133, 0.0134635348)        # medium
    LOW = (2.580271, 0.127628236, 0.0252748384)          # low
    # Official Input/Output verification data (NTCM-G description, Annex D):
    # (a0, a1, a2, doy, ut, stn_lon, stn_lat, stn_h, sat_lon, sat_lat, sat_h, STEC)
    vectors = [
        (HI..., 105, 0, -62.34, 82.49, 78.11, 8.23, 54.29, 20281546.18, 33.7567),
        (HI..., 105, 12, -62.34, 82.49, 78.11, 81.09, 35.20, 20278071.03, 65.0500),
        (HI..., 105, 20, -52.81, 5.25, -25.76, 10.94, 44.72, 20450566.19, 252.0204),
        (HI..., 105, 16, -52.81, 5.25, -25.76, -70.26, 50.63, 20043030.72, 216.2278),
        (MED..., 105, 4, 40.19, -3.00, -23.32, 107.19, -10.65, 19943686.06, 108.8940),
        (MED..., 105, 20, 115.89, -31.80, 12.78, 131.65, -31.56, 20066111.12, 7.5508),
        (LOW..., 105, 0, 141.13, 39.14, 117.00, 165.14, -13.93, 20181976.50, 51.5270),
        (LOW..., 105, 20, -155.46, 19.80, 3754.69, -82.52, 20.64, 19937791.48, 67.4750),
    ]
    @testset "official Annex D test vectors" begin
        for v in vectors
            a0, a1, a2, doy, ut, slon, slat, sh, blon, blat, bh, expected = v
            u = toecef(slat, slon, sh)
            s = toecef(blat, blon, bh)
            el, az = PositionVelocityTime._elevation_azimuth(ENUfromECEF(u, wgs84), s)
            stec = PositionVelocityTime.ntcm_g_stec(
                el,
                az,
                LLAfromECEF(wgs84)(u),
                doy,
                ut,
                a0,
                a1,
                a2,
            )
            @test stec ≈ expected atol = 1e-3   # published values rounded to 4 dp
        end
    end

    @testset "Azpar (Eq. 2)" begin
        # |√(a0² + 1633.33·a1² + 4802000·a2² + 3266.67·a0·a2)|, in sfu
        @test PositionVelocityTime._azpar(HI...) ≈ 244.007 atol = 1e-2
        @test PositionVelocityTime._azpar(0.0, 0.0, 0.0) == 0.0
    end

    @testset "GST week/TOW → day-of-year and UT" begin
        # 1 week + 12 h past the GST epoch (1999-08-22 00:00 UTC) → 1999-08-29, 12:00
        doy, ut = PositionVelocityTime._galileo_doy_and_ut(1, 12 * 3600)
        @test doy == dayofyear(Date(1999, 8, 29))
        @test ut ≈ 12.0
    end

    @testset "selection and delay from a Galileo decoder" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        sat = ECEF(user[1] + 1.0e7, user[2] + 1.0e7, user[3] + 2.0e7)
        # Galileo decoder carrying broadcast a_i0..a_i2 and a week number
        gal = GNSSDecoderState(GalileoE1B(), 1)
        gal = GNSSDecoder.GNSSDecoderState(
            gal;
            data = GNSSDecoder.GalileoE1BData(
                gal.data;
                a_i0 = 236.831641,
                a_i1 = -0.39362878,
                a_i2 = 0.00402826613,
                WN = 1100,
            ),
        )
        state = SatelliteState(;
            decoder = gal,
            system = GalileoE1B(),
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        correction = PositionVelocityTime.select_ionospheric_correction([state])
        @test correction isa PositionVelocityTime.NTCMGParams
        el, az = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat)
        delay = PositionVelocityTime.ionospheric_delay(
            correction,
            GalileoE1B(),
            el,
            az,
            LLAfromECEF(wgs84)(user),
            200000.0,
        )
        @test delay > 0.0
        @test delay < 100.0   # a sane L1/E1 ionospheric delay magnitude (metres)
    end
end
