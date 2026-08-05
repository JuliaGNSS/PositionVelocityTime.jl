# Independent reimplementation of the tropospheric model — Saastamoinen zenith
# delays (RTKLIB `tropmodel` at zenith) mapped by the Niell (1996) mapping
# functions (RTKLIB `tropmapf`/`nmf`) — to cross-check the package implementation.
# Written differently on purpose: latitudes in degrees, table rows as vectors with
# `searchsortedlast` interpolation, and the seasonal term as a phase shift.

const NMF_LATITUDE_NODES = 15.0:15.0:75.0

# Niell (1996) tables 3 and 4 (corrected paper), rows (a, b, c) per column latitude
const NMF_HYDROSTATIC_AVG = (
    [1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3],
    [2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3],
    [62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3],
)
const NMF_HYDROSTATIC_AMP = (
    [0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5],
    [0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5],
    [0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5],
)
const NMF_WET_REF = (
    [5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4],
    [1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3],
    [4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2],
)
const NMF_HEIGHT_REF = (2.53e-5, 5.49e-3, 1.14e-3)

function nmf_interpolate(row, abs_lat_deg)
    abs_lat_deg <= first(NMF_LATITUDE_NODES) && return first(row)
    abs_lat_deg >= last(NMF_LATITUDE_NODES) && return last(row)
    i = searchsortedlast(NMF_LATITUDE_NODES, abs_lat_deg)
    f = (abs_lat_deg - NMF_LATITUDE_NODES[i]) / step(NMF_LATITUDE_NODES)
    return (1 - f) * row[i] + f * row[i+1]
end

# Normalized continued fraction: 1 at zenith, finite at the horizon
function continued_fraction(sinel, (a, b, c))
    return (1 + a / (1 + b / (1 + c))) / (sinel + a / (sinel + b / (sinel + c)))
end

function niell_reference(lat_deg, height, elevation, doy)
    el = max(elevation, 0.0)
    # Seasonal phase: day-of-year cosine anchored at doy 28, southern hemisphere
    # shifted half a year (expressed here as a phase offset of π)
    phase = 2π * (doy - 28) / 365.25 + (lat_deg < 0 ? π : 0.0)
    alat = abs(lat_deg)
    hydrostatic_abc = ntuple(3) do i
        nmf_interpolate(NMF_HYDROSTATIC_AVG[i], alat) -
        nmf_interpolate(NMF_HYDROSTATIC_AMP[i], alat) * cos(phase)
    end
    wet_abc = ntuple(i -> nmf_interpolate(NMF_WET_REF[i], alat), 3)
    # Height correction, frozen below 0.05 rad where its 1/sin(el) diverges
    el_ht = max(el, 0.05)
    dm = (1 / sin(el_ht) - continued_fraction(sin(el_ht), NMF_HEIGHT_REF)) * height / 1000
    mh = continued_fraction(sin(el), hydrostatic_abc) + dm
    mw = continued_fraction(sin(el), wet_abc)
    return mh, mw
end

function slant_delay_reference(lat_deg, height, elevation, humidity, doy)
    (height < -100.0 || height > 1.0e4) && return 0.0
    h = max(height, 0.0)
    P = 1013.25 * (1 - 2.2557e-5 * h)^5.2568
    T = 288.16 - 6.5e-3 * h
    e = 6.108 * humidity * exp((17.15 * T - 4684.0) / (T - 38.45))
    zhd = 0.0022768 * P / (1 - 0.00266 * cosd(2 * lat_deg) - 0.00028 * h / 1000)
    zwd = 0.002277 * (1255.0 / T + 0.05) * e
    mh, mw = niell_reference(lat_deg, h, elevation, doy)
    return zhd * mh + zwd * mw
end

@testset "Tropospheric model (Saastamoinen zenith × Niell mapping)" begin
    @testset "matches independent reference" begin
        for lat in (0.0, 48.0, -67.0, 80.0),
            h in (0.0, 550.0, 3000.0, 9000.0),
            el in deg2rad.((-10.0, 0.0, 1.0, 3.0, 5.0, 15.0, 45.0, 90.0)),
            humi in (0.0, 0.7, 1.0),
            doy in (1, 28, 119, 211, 365)

            got = PositionVelocityTime.tropospheric_delay(
                el,
                LLA(lat, 11.0, h),
                doy;
                humidity = humi,
            )
            ref = slant_delay_reference(lat, h, el, humi, doy)
            @test got ≈ ref rtol = 1e-12

            mh, mw = PositionVelocityTime.niell_mapping_functions(
                deg2rad(lat),
                max(h, 0.0),
                el,
                doy,
            )
            mh_ref, mw_ref = niell_reference(lat, max(h, 0.0), el, doy)
            @test mh ≈ mh_ref rtol = 1e-12
            @test mw ≈ mw_ref rtol = 1e-12
        end
    end

    @testset "zenith delay unchanged (regression)" begin
        lat, humi = 48.0, 0.7
        # Both mapping functions are exactly 1 at zenith (normalized continued
        # fraction) and the height correction vanishes there, so the zenith slant
        # delay is exactly the sum of the zenith delays — as before the remap.
        mh, mw = PositionVelocityTime.niell_mapping_functions(
            deg2rad(lat),
            550.0,
            π / 2,
            100,
        )
        @test mh ≈ 1.0 rtol = 1e-12
        @test mw ≈ 1.0 rtol = 1e-12
        zhd, zwd = PositionVelocityTime.saastamoinen_zenith_delays(deg2rad(lat), 0.0, humi)
        ztd = PositionVelocityTime.tropospheric_delay(π / 2, LLA(lat, 11.0, 0.0), 100)
        @test ztd ≈ zhd + zwd rtol = 1e-12
        # Standard atmosphere at sea level, 70 % RH: ZHD ≈ 2.31 m, ZWD ≈ 0.12 m (#62)
        @test ztd ≈ 2.43 atol = 0.05
        @test 2.2 < ztd < 2.6
        # Delay decreases with user height (thinner atmosphere above)
        @test PositionVelocityTime.tropospheric_delay(π / 2, LLA(lat, 11.0, 3000.0), 100) <
              ztd
        # Wet component grows with humidity; dry-only is the floor
        @test PositionVelocityTime.tropospheric_delay(
            π / 2,
            LLA(lat, 11.0, 0.0),
            100;
            humidity = 1.0,
        ) > PositionVelocityTime.tropospheric_delay(
            π / 2,
            LLA(lat, 11.0, 0.0),
            100;
            humidity = 0.0,
        )
    end

    @testset "low-elevation values match the Niell-mapped references of #62" begin
        # Issue #62's comparison table: lat 48°, sea level, 70 % RH, average (no
        # seasonal) coefficients — the doy is chosen so the seasonal cosine is ≈ 0.
        lla = LLA(48.0, 11.0, 0.0)
        doy = 119  # 28 + 365.25/4 → cos ≈ 0
        for (el_deg, mh_expected, slant_expected) in (
            (15.0, 3.79, 9.2),
            (10.0, 5.54, 13.5),
            (5.0, 10.09, 24.6),
            (3.0, 14.54, 35.5),
            (0.0, 36.46, 90.9),
        )
            mh, _ = PositionVelocityTime.niell_mapping_functions(
                deg2rad(48.0),
                0.0,
                deg2rad(el_deg),
                doy,
            )
            @test mh ≈ mh_expected rtol = 0.03
            slant = PositionVelocityTime.tropospheric_delay(deg2rad(el_deg), lla, doy)
            @test slant ≈ slant_expected rtol = 0.03
        end
        # The flat-slab 1/sin(el) mapping the Niell functions replaced over-predicted
        # by ≈ 14 % at 5° — make sure that bias is gone, not just reduced.
        ztd = PositionVelocityTime.tropospheric_delay(π / 2, lla, doy)
        flat_slab_5deg = ztd / sind(5.0)
        niell_5deg = PositionVelocityTime.tropospheric_delay(deg2rad(5.0), lla, doy)
        @test niell_5deg < 0.90 * flat_slab_5deg
    end

    @testset "finite and monotone down to the horizon" begin
        for h in (0.0, 550.0, 3000.0)
            lla = LLA(48.0, 11.0, h)
            delays = map(range(0.0, π / 2; length = 200)) do el
                PositionVelocityTime.tropospheric_delay(el, lla, 100)
            end
            @test all(isfinite, delays)
            # Strictly decreasing in elevation, no clamp needed above the horizon
            @test all(diff(delays) .< 0)
            # A grazing ray at sea level is ≈ 90 m — large, but finite and physical,
            # unlike the diverging 1/sin(el) mapping (#62)
            @test delays[1] < 120.0
            # Below the horizon the horizon value is reused: continuous, finite
            at_horizon = delays[1]
            for el in (-1e-6, -0.05, deg2rad(-10.0), -π / 2)
                @test PositionVelocityTime.tropospheric_delay(el, lla, 100) == at_horizon
            end
        end
    end

    @testset "seasonal term and hemispheres" begin
        el = deg2rad(5.0)
        # The hydrostatic coefficients swing over the year at mid/high latitude...
        mh_winter, mw_winter =
            PositionVelocityTime.niell_mapping_functions(deg2rad(48.0), 0.0, el, 28)
        mh_summer, mw_summer =
            PositionVelocityTime.niell_mapping_functions(deg2rad(48.0), 0.0, el, 211)
        @test mh_winter != mh_summer
        # ...the wet coefficients are latitude-only (Niell 1996)
        @test mw_winter == mw_summer
        # Southern hemisphere is the northern one shifted half a year
        mh_south, _ =
            PositionVelocityTime.niell_mapping_functions(deg2rad(-48.0), 0.0, el, 28)
        @test mh_south ≈ mh_summer rtol = 1e-6
        # At the equator seasons cancel out of the tables (amplitude column is
        # constant below 15°): no day-of-year dependence at all
        mh_eq_1, _ = PositionVelocityTime.niell_mapping_functions(0.0, 0.0, el, 28)
        mh_eq_2, _ = PositionVelocityTime.niell_mapping_functions(0.0, 0.0, el, 211)
        @test mh_eq_1 == mh_eq_2
    end

    @testset "out-of-range heights return zero" begin
        # User height outside the model's valid range (−100 m … 10 km) returns exactly
        # zero: too high, then too low.
        @test PositionVelocityTime.tropospheric_delay(
            deg2rad(45.0),
            LLA(48.0, 11.0, 1.5e4),
            100,
        ) == 0.0
        @test PositionVelocityTime.tropospheric_delay(
            deg2rad(45.0),
            LLA(48.0, 11.0, -200.0),
            100,
        ) == 0.0
    end

    @testset "tropospheric_delay from line-of-sight geometry" begin
        user = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 550.0))
        lla = LLAfromECEF(wgs84)(user)
        # Satellite roughly overhead → near-zenith, small slant delay (~2 m)
        sat_up = ECEFfromLLA(wgs84)(LLA(48.0, 11.0, 2.0e7))
        el_up, _ = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat_up)
        d_up = PositionVelocityTime.tropospheric_delay(el_up, lla, 100)
        @test 2.0 < d_up < 2.6
        # A low-elevation satellite has a larger slant delay than the near-zenith one
        sat_low = ECEF(user[1] + 2.0e7, user[2] + 2.0e6, user[3] + 2.0e6)
        el_low, _ = PositionVelocityTime._elevation_azimuth(ENUfromECEF(user, wgs84), sat_low)
        @test PositionVelocityTime.tropospheric_delay(el_low, lla, 100) > d_up
        # Frequency independence is implicit: tropospheric_delay takes no system/freq.
    end

    @testset "day of year from system week and time of week" begin
        # GPS week 0 starts 1980-01-06 → doy 6; half a year of weeks later lands
        # mid-year. The ≤ 18 s GPS−UTC offset is irrelevant at day resolution.
        @test PositionVelocityTime._day_of_year(GPST(), 0, 0.0) == 6
        @test PositionVelocityTime._day_of_year(GPST(), 0, 86400.0 * 10) == 16
        # GST week 0 starts 1999-08-21 23:59:47 UTC (doy 233)
        @test PositionVelocityTime._day_of_year(GST(), 0, 14.0) == 234
    end
end
