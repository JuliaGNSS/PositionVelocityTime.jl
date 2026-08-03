# Independent reimplementation of the Saastamoinen model (RTKLIB `tropmodel`),
# written with sin(elevation) instead of cos(z) and explicit intermediate terms,
# to cross-check the package's `saastamoinen_delay`.
function saastamoinen_reference(lat, height, elevation, humidity)
    (height < -100.0 || height > 1.0e4) && return 0.0
    h = max(height, 0.0)
    P = 1013.25 * (1 - 2.2557e-5 * h)^5.2568
    T = 288.16 - 6.5e-3 * h
    e = 6.108 * humidity * exp((17.15 * T - 4684.0) / (T - 38.45))
    el = max(elevation, 0.05)                         # Orekit's low elevation bound
    mapping = 1 / sin(el)                             # = 1/cos(z), z = π/2 - el
    zhd = 0.0022768 * P / (1 - 0.00266 * cos(2lat) - 0.00028 * h / 1000)
    zwd = 0.002277 * (1255.0 / T + 0.05) * e
    return (zhd + zwd) * mapping
end

@testset "Saastamoinen tropospheric model" begin
    @testset "matches independent reference" begin
        for lat in deg2rad.((0.0, 48.0, -67.0)),
            h in (0.0, 550.0, 3000.0, 9000.0),
            el in deg2rad.((-10.0, 0.0, 1.0, 5.0, 15.0, 45.0, 90.0)),
            humi in (0.0, 0.5, 0.7, 1.0)

            got = PositionVelocityTime.saastamoinen_delay(lat, h, el, humi)
            ref = saastamoinen_reference(lat, h, el, humi)
            @test got ≈ ref rtol = 1e-12
        end
    end

    @testset "physical behaviour" begin
        lat, humi = deg2rad(48.0), 0.7
        zenith = PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(90.0), humi)
        # Zenith total delay at sea level is a couple of metres
        @test 2.2 < zenith < 2.6
        # Obliquity: lower elevation → larger slant delay (≈ 1/sin(el))
        d30 = PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(30.0), humi)
        d10 = PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(10.0), humi)
        @test d30 > zenith
        @test d10 > d30
        @test d30 ≈ zenith / sin(deg2rad(30.0)) rtol = 1e-12
        # Delay decreases with user height (thinner atmosphere above)
        @test PositionVelocityTime.saastamoinen_delay(lat, 3000.0, deg2rad(90.0), humi) <
              zenith
        # Wet component grows with humidity; dry-only is the floor
        @test PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(90.0), 1.0) >
              PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(90.0), 0.0)
    end

    @testset "low elevations are bounded" begin
        lat, humi = deg2rad(48.0), 0.7
        bound = PositionVelocityTime._LOW_ELEVATION_THRESHOLD
        @test bound == 0.05  # Orekit's DEFAULT_LOW_ELEVATION_THRESHOLD
        at_bound = PositionVelocityTime.saastamoinen_delay(lat, 0.0, bound, humi)
        # Everything at or below the bound — grazing, at the horizon, below it — reuses
        # the delay computed at the bound, so the correction stays finite and continuous
        # across the horizon instead of diverging or dropping to zero.
        for el in (bound, 1e-3, 0.0, deg2rad(-5.0), -π / 2)
            @test PositionVelocityTime.saastamoinen_delay(lat, 0.0, el, humi) == at_bound
        end
        # The cap is just the 1/sin(el) mapping evaluated at the bound — ≈ 20 × zenith,
        # a few tens of metres instead of the hundreds an unbounded mapping gives. It is
        # a bound, not an accurate low-elevation delay (see #62).
        zenith = PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(90.0), humi)
        @test at_bound ≈ zenith / sin(bound) rtol = 1e-12
        @test at_bound < 50.0
        # Above the bound the model is untouched and still monotonic in elevation
        just_above = PositionVelocityTime.saastamoinen_delay(lat, 0.0, 0.06, humi)
        @test just_above < at_bound
        @test PositionVelocityTime.saastamoinen_delay(lat, 0.0, deg2rad(5.0), humi) <
              just_above
    end

    @testset "out-of-range heights return zero" begin
        lat, humi = deg2rad(48.0), 0.7
        # User height outside the model's valid range (−100 m … 10 km) returns exactly
        # zero: too high, then too low.
        @test PositionVelocityTime.saastamoinen_delay(lat, 1.5e4, deg2rad(45.0), humi) == 0.0
        @test PositionVelocityTime.saastamoinen_delay(lat, -200.0, deg2rad(45.0), humi) == 0.0
    end

    @testset "tropospheric_delay from line-of-sight geometry" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        lla = LLAfromECEF(wgs84)(user)
        # Satellite roughly overhead → near-zenith, small slant delay (~2 m)
        sat_up = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 2.0e7))
        el_up, _ = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat_up)
        d_up = PositionVelocityTime.tropospheric_delay(el_up, lla)
        @test 2.0 < d_up < 2.6
        # A low-elevation satellite has a larger slant delay than the near-zenith one
        sat_low = ECEF(user[1] + 2.0e7, user[2] + 2.0e6, user[3] + 2.0e6)
        el_low, _ = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat_low)
        @test PositionVelocityTime.tropospheric_delay(el_low, lla) > d_up
        # Frequency independence is implicit: tropospheric_delay takes no system/freq.
    end
end
