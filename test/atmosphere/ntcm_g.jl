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
            el, az = calc_elevation_and_azimuth(ENUfromECEF(u, wgs84), s)
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
        @test PositionVelocityTime._ntcm_g_azpar(HI...) ≈ 244.007 atol = 1e-2
        @test PositionVelocityTime._ntcm_g_azpar(0.0, 0.0, 0.0) == 0.0
    end

    @testset "the day of year and UT off the epoch's TAI calendar label" begin
        # 1 week + 12 h past the GST epoch (1999-08-22T00:00:19 TAI) → 1999-08-29
        # at 12:00:19 on the TAI calendar label the date is read off.
        lla = LLA(48.0, 11.0, 550.0)
        el, az, f = deg2rad(30), deg2rad(60), carrier_hz(GalileoE1B())
        delay = PositionVelocityTime.ionospheric_delay(GNSSDecoder.NeQuickGParameters(HI...),
            f, el, az, lla, get_week_start_epoch(GST(), 1) + 12 * 3600)
        stec = PositionVelocityTime.ntcm_g_stec(el, az, lla,
            dayofyear(Date(1999, 8, 29)), 12 + 19 / 3600, HI...)
        @test delay ≈ 40.3 / f^2 * stec * 1.0e16
    end

    @testset "selection and delay from a Galileo decoder" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        sat = ECEF(user[1] + 1.0e7, user[2] + 1.0e7, user[3] + 2.0e7)
        # Galileo decoder carrying broadcast a_i0..a_i2 and a week number
        gal = GNSSDecoderState(GalileoE1B(), 1)
        gal = GNSSDecoder.GNSSDecoderState(
            gal;
            data = GNSSDecoder.GalileoINAVData(
                gal.data;
                a_i0 = 236.831641,
                a_i1 = -0.39362878,
                a_i2 = 0.00402826613,
                WN = 1100,
            ),
        )
        state = SatelliteState(;
            decoder = gal,
            signal = GalileoE1B(),
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        parameters = PositionVelocityTime.select_ionospheric_parameters(signal_group([state]))
        @test parameters isa GNSSDecoder.NeQuickGParameters
        el, az = calc_elevation_and_azimuth(ENUfromECEF(user, wgs84), sat)
        delay = PositionVelocityTime.ionospheric_delay(
            parameters,
            carrier_hz(GalileoE1B()),
            el,
            az,
            LLAfromECEF(wgs84)(user),
            week_start(GST(), 1100) + 200000.0,
        )
        @test delay > 0.0
        @test delay < 100.0   # a sane L1/E1 ionospheric delay magnitude (metres)
    end
end
