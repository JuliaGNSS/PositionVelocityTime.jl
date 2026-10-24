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

# Literal transcription of BDS-SIS-ICD-B1I-3.0 §5.2.4.7 (inputs in radians), used to
# cross-check `beidou_klobuchar_group_delay`. Powers are summed explicitly where the
# implementation uses Horner's scheme, and the clamps are spelled as the ICD's cases.
function bds_klobuchar_reference(lat, lon, az, el, t_E, α, β)
    R, h = 6378.0, 375.0                              # km, as printed in the ICD
    ψ = π / 2 - el - asin(R / (R + h) * cos(el))
    φ_M = asin(sin(lat) * cos(ψ) + cos(lat) * sin(ψ) * cos(az))
    λ_M = lon + asin(sin(ψ) * sin(az) / cos(φ_M))
    t = mod(t_E + λ_M * 43200 / π, 86400)
    A_2 = sum(α[n+1] * (φ_M / π)^n for n = 0:3)
    A_2 = A_2 < 0 ? 0.0 : A_2
    A_4 = sum(β[n+1] * (φ_M / π)^n for n = 0:3)
    A_4 = A_4 >= 172800 ? 172800.0 : (A_4 < 72000 ? 72000.0 : A_4)
    I_z = abs(t - 50400) < A_4 / 4 ? 5e-9 + A_2 * cos(2π * (t - 50400) / A_4) : 5e-9
    return I_z / sqrt(1 - (R / (R + h) * cos(el))^2)
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
    function galileo_decoder_with(a_i0, a_i1, a_i2, WN; prn = 1)
        dec = GNSSDecoderState(GalileoE1B(), prn)
        GNSSDecoder.GNSSDecoderState(
            dec;
            data = GNSSDecoder.GalileoINAVData(dec.data; a_i0, a_i1, a_i2, WN),
        )
    end
    mkstate(dec, sys) = SatelliteState(;
        decoder = dec,
        signal = sys,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )

    @testset "parameter extraction from decoders" begin
        # Fresh decoders carry no coefficients
        @test PositionVelocityTime.get_klobuchar_parameters(GNSSDecoderState(GPSL1CA(), 1)) ===
              nothing
        @test PositionVelocityTime.get_ntcm_g_parameters(GNSSDecoderState(GalileoE1B(), 1)) ===
              nothing
        # Broadcast α/β are extracted into KlobucharParameters
        kp = PositionVelocityTime.get_klobuchar_parameters(gps_decoder_with(α, β))
        @test kp isa PositionVelocityTime.KlobucharParameters
        @test (kp.α_0, kp.α_1, kp.α_2, kp.α_3) == α
        @test (kp.β_0, kp.β_1, kp.β_2, kp.β_3) == β
        # Broadcast a_i0..a_i2 + WN are extracted into NTCMGParameters
        np = PositionVelocityTime.get_ntcm_g_parameters(
            galileo_decoder_with(236.831641, -0.39362878, 0.00402826613, 1100),
        )
        @test np isa PositionVelocityTime.NTCMGParameters
        @test np.a_i0 == 236.831641
        @test np.week == 1100
        # Wrong-system decoders never yield the other model's params
        @test PositionVelocityTime.get_ntcm_g_parameters(gps_decoder_with(α, β)) === nothing
        @test PositionVelocityTime.get_klobuchar_parameters(
            galileo_decoder_with(100.0, 0.0, 0.0, 1100),
        ) === nothing
    end

    @testset "constellation-wide model selection" begin
        # Distinct PRNs within a constellation: a signal group is keyed by PRN, so the
        # bare and the coefficient-carrying satellite of one signal have to be two
        # satellites rather than the same one twice.
        gps_bare = mkstate(GNSSDecoderState(GPSL1CA(), 2), GPSL1CA())
        gps_klob = mkstate(gps_decoder_with(α, β), GPSL1CA())
        gal_bare = mkstate(GNSSDecoderState(GalileoE1B(), 2), GalileoE1B())
        gal_ntcm = mkstate(galileo_decoder_with(121.13, 0.35, 0.013, 1100), GalileoE1B())

        select = PositionVelocityTime.select_ionospheric_parameters
        # Nothing decoded → no correction
        @test select((signal_group([gps_bare]), signal_group([gal_bare]))) === nothing
        # Only Klobuchar → Klobuchar. Both satellites are on one signal, so they are
        # one group.
        @test select(signal_group([gps_klob, gps_bare])) isa
              PositionVelocityTime.KlobucharParameters
        # Only Galileo → NTCM-G
        @test select(signal_group([gal_ntcm, gal_bare])) isa
              PositionVelocityTime.NTCMGParameters
        # Both available → NTCM-G wins (more accurate)
        @test select((signal_group([gps_klob]), signal_group([gal_ntcm]))) isa
              PositionVelocityTime.NTCMGParameters
    end

    @testset "selected model applied to every satellite" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        sat = ECEF(user[1] + 1.0e7, user[2] + 1.0e7, user[3] + 2.0e7)
        klob = PositionVelocityTime.KlobucharParameters(α..., β...)
        # Line-of-sight geometry (same user/sat for every call below)
        el, az = calc_elevation_and_azimuth(ENUfromECEF(user, wgs84), sat)
        lla = LLAfromECEF(wgs84)(user)

        # No correction → exactly zero for any system
        @test PositionVelocityTime.ionospheric_delay(
            nothing,
            carrier_hz(GPSL1CA()),
            el,
            az,
            lla,
            50400.0,
        ) == 0.0
        @test PositionVelocityTime.ionospheric_delay(
            nothing,
            carrier_hz(GalileoE1B()),
            el,
            az,
            lla,
            50400.0,
        ) == 0.0
        # Klobuchar applied to GPS *and* Galileo (E1 shares the L1 frequency, so
        # the delay is identical for both systems)
        d_gps =
            PositionVelocityTime.ionospheric_delay(klob, carrier_hz(GPSL1CA()), el, az, lla, 50400.0)
        d_gal =
            PositionVelocityTime.ionospheric_delay(klob, carrier_hz(GalileoE1B()), el, az, lla, 50400.0)
        @test d_gps > 0.0
        @test d_gps ≈ d_gal rtol = 1e-12

        # The delay scales as 1/f², so the same coefficients applied on a lower band
        # (here GPS L5, 1176.45 MHz) give the correct larger delay, not the L1 value.
        d_l5 = PositionVelocityTime.ionospheric_delay(klob, carrier_hz(GPSL5I()), el, az, lla, 50400.0)
        # ratio of two Hz quantities is dimensionless
        f_ratio =
            GNSSSignals.get_center_frequency(GPSL1CA()) /
            GNSSSignals.get_center_frequency(GPSL5I())
        @test d_l5 ≈ d_gps * f_ratio^2 rtol = 1e-12
        @test d_l5 > d_gps
    end
end

@testset "BDS Klobuchar variant (BeiDou B1I)" begin
    # A plausible broadcast set. The BeiDou coefficients come in the same s·π⁻ⁿ
    # units as the GPS ones, so the numbers are interchangeable between variants —
    # which is exactly why sharing them must not mean sharing the algorithm.
    α = (1.0e-8, 2.0e-8, -3.0e-8, 4.0e-8)
    β = (90112.0, 16384.0, -196608.0, 131072.0)
    PVT = PositionVelocityTime

    @testset "beidou_klobuchar_group_delay matches the ICD transcription" begin
        for lat in deg2rad.((10.0, 48.0, -33.0)),
            lon in deg2rad.((-120.0, 11.0, 150.0)),
            az in deg2rad.((0.0, 95.0, 270.0)),
            el in deg2rad.((5.0, 25.0, 89.0)),
            t in (0.0, 43200.0, 50400.0, 80000.0)

            got = PVT.beidou_klobuchar_group_delay(lat, lon, el, az, t, α, β)
            @test got ≈ bds_klobuchar_reference(lat, lon, az, el, t, α, β) rtol = 1e-12
        end
    end

    @testset "zenith at 14:00 local time reduces to 5 ns + A₂" begin
        # At E = π/2 the pierce point is the user (ψ = 0), the obliquity is exactly
        # 1, and local 14:00 sits at the cosine's peak, so the ICD's whole chain
        # collapses to 5·10⁻⁹ + Σ αₙ(φ_u/π)ⁿ — a closed form the implementation
        # must hit whatever the azimuth.
        lat = deg2rad(30.0)
        φ = lat / π
        expected = 5.0e-9 + α[1] + α[2] * φ + α[3] * φ^2 + α[4] * φ^3
        for az in (0.0, 2.1)
            got = PVT.beidou_klobuchar_group_delay(lat, 0.0, π / 2, az, 50400.0, α, β)
            @test got ≈ expected rtol = 1e-14
        end
    end

    @testset "obliquity, floors and clamps" begin
        lat, lon, az = deg2rad(48.0), deg2rad(11.0), deg2rad(120.0)
        # With A₂ = 0 the vertical delay is the 5 ns night floor at any time of
        # day, so the slant result isolates the ICD's exact obliquity factor.
        no_amp = (0.0, 0.0, 0.0, 0.0)
        for el in deg2rad.((5.0, 30.0, 60.0, 90.0))
            got = PVT.beidou_klobuchar_group_delay(lat, lon, el, az, 12345.0, no_amp, β)
            @test got ≈ 5.0e-9 / sqrt(1 - (6378.0 / 6753.0 * cos(el))^2) rtol = 1e-14
        end
        # A₂ floor: a negative amplitude polynomial is floored at zero, i.e. the
        # day-time value collapses onto the night floor.
        neg = PVT.beidou_klobuchar_group_delay(
            lat, lon, 0.6, az, 50400.0, (-1.0e-8, 0.0, 0.0, 0.0), β)
        flat = PVT.beidou_klobuchar_group_delay(lat, lon, 0.6, az, 50400.0, no_amp, β)
        @test neg == flat
        # A₄ clamps, from both sides: a period polynomial beyond a bound behaves
        # exactly as one pinned at that bound. 47000 s BDT puts the pierce point's
        # local time inside the day window for either clamped period.
        delay_for(β) = PVT.beidou_klobuchar_group_delay(lat, lon, 0.6, az, 47000.0, α, β)
        @test delay_for((1.0e9, 0.0, 0.0, 0.0)) == delay_for((172800.0, 0.0, 0.0, 0.0))
        @test delay_for((-1.0e9, 0.0, 0.0, 0.0)) == delay_for((72000.0, 0.0, 0.0, 0.0))
        @test delay_for((1.0e9, 0.0, 0.0, 0.0)) != delay_for((-1.0e9, 0.0, 0.0, 0.0))
    end

    @testset "not the GPS algorithm" begin
        # Same coefficients, same line of sight: the two variants must disagree —
        # geographic vs geomagnetic latitude alone moves the polynomials' argument,
        # and the exact obliquity differs from the GPS fit by ~1 % at mid elevation.
        # This is the regression guard for the bug where the BeiDou broadcast set
        # was fed through the IS-GPS-200 algorithm.
        lat, lon = deg2rad(48.0), deg2rad(11.0)
        el, az, t = deg2rad(30.0), deg2rad(200.0), 50000.0
        bds = PVT.beidou_klobuchar_group_delay(lat, lon, el, az, t, α, β)
        gps = PVT.klobuchar_group_delay(lat / π, lon / π, el / π, az / π, t, α, β)
        @test !isapprox(bds, gps; rtol = 1e-3)
    end

    @testset "the delay is defined at B1I and rescales by 1/f²" begin
        p = PVT.BeiDouKlobucharParameters(α..., β...)
        lla = LLA(48.0, 11.0, 550.0)
        el, az, t = deg2rad(35.0), deg2rad(120.0), 50000.0
        # At B1I itself the metre conversion is exactly c times the ICD's I_B1I.
        d_b1i = PVT.ionospheric_delay(p, carrier_hz(BeiDouB1I()), el, az, lla, t)
        seconds = PVT.beidou_klobuchar_group_delay(
            deg2rad(lla.lat), deg2rad(lla.lon), el, az, t, α, β)
        @test d_b1i ≈ 299792458.0 * seconds rtol = 1e-14
        # Every other carrier gets that delay rescaled by 1/f² from B1I — including
        # non-BeiDou satellites, since one model corrects the whole solve.
        f(s) = GNSSSignals.get_center_frequency(s)
        d_b3i = PVT.ionospheric_delay(p, carrier_hz(BeiDouB3I()), el, az, lla, t)
        d_l1 = PVT.ionospheric_delay(p, carrier_hz(GPSL1CA()), el, az, lla, t)
        @test d_b3i / d_b1i ≈ (f(BeiDouB1I()) / f(BeiDouB3I()))^2 rtol = 1e-12
        @test d_l1 / d_b1i ≈ (f(BeiDouB1I()) / f(GPSL1CA()))^2 rtol = 1e-12
        @test d_l1 < d_b1i < d_b3i    # 1575.42 > 1561.098 > 1268.52 MHz
    end
end
