# ===========================================================================
#  NTCM-G — Galileo single-frequency ionospheric model
#  (European GNSS (Galileo) NTCM-G Ionospheric Model Description, Issue 1.0,
#   May 2022). Driven by the broadcast Effective Ionisation Level coefficients
#   a_i0, a_i1, a_i2. Returns slant TEC (TECU); equation numbers below refer to
#   that document.
# ===========================================================================

# NTCM-G model coefficients k1..k12 (Table 3)
const _NTCM_G_K = (
    0.92519,
    0.16951,
    0.00443,
    0.06626,
    0.00899,
    0.21289,
    -0.15414,
    -0.38439,
    1.14023,
    1.20556,
    1.41808,
    0.13985,
)
const _NTCM_G_RE = 6371.0               # Earth mean radius [km] (Table 2)
const _NTCM_G_HI = 450.0                # ionospheric pierce point height [km]
const _NTCM_G_GNP_LAT = deg2rad(79.74)  # geomagnetic North pole latitude
const _NTCM_G_GNP_LON = deg2rad(-71.78) # geomagnetic North pole longitude

# Effective Ionisation Level Azpar [sfu] from the broadcast coefficients (Eq. 2).
function _ntcm_g_azpar(a_i0, a_i1, a_i2)
    radicand = a_i0^2 + 1633.33 * a_i1^2 + 4802000.0 * a_i2^2 + 3266.67 * a_i0 * a_i2
    return sqrt(max(radicand, 0.0))
end

# Ionospheric pierce point geographic latitude/longitude [rad] (Eq. 24-26),
# given user geodetic lat/lon and satellite elevation/azimuth [rad].
function _ntcm_g_pierce_point(φ_u, λ_u, elevation, azimuth)
    ψ = π / 2 - elevation -
        asin(_NTCM_G_RE / (_NTCM_G_RE + _NTCM_G_HI) * cos(elevation))
    # Clamp to asin's domain: the φ_pp argument is a unit dot-product that can
    # overshoot ±1 by a rounding ulp, and the λ_pp argument genuinely diverges as
    # the pierce point approaches a pole (cos(φ_pp) → 0).
    φ_pp = asin(clamp(sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth), -1.0, 1.0))
    λ_pp = λ_u + asin(clamp(sin(ψ) * sin(azimuth) / cos(φ_pp), -1.0, 1.0))
    return φ_pp, λ_pp
end

# Sun's declination [rad] for the day of year (Eq. 28).
_ntcm_g_sun_declination(doy) = deg2rad(23.44) * sin(deg2rad(0.9856 * (doy - 80.7)))

# Modified Single Layer Model mapping function (Eq. 32-33).
function _ntcm_g_mapping_function(elevation)
    sinz = _NTCM_G_RE / (_NTCM_G_RE + _NTCM_G_HI) *
           sin(0.9782 * (π / 2 - elevation))
    return 1.0 / sqrt(1.0 - sinz^2)
end

"""
    ntcm_g_vtec(φ_pp, λ_pp, day_of_year, ut, azpar) -> Float64

Vertical TEC in TECU at the ionospheric pierce point (geographic latitude/longitude
`φ_pp`/`λ_pp` in radians) on `day_of_year`, at universal time `ut` in hours, and
Effective Ionisation Level `azpar` in solar flux units. Implements the NTCM-G model
`VTEC = F1·F2·F3·F4·F5` (Eq. 3-15).
"""
function ntcm_g_vtec(φ_pp, λ_pp, day_of_year, ut, azpar)
    k = _NTCM_G_K
    lt = ut + rad2deg(λ_pp) / 15                      # local time [h] (Eq. 27)
    δ = _ntcm_g_sun_declination(day_of_year)

    # Solar zenith angle dependence (Eq. 29-30)
    cosχ3 = cos(φ_pp - δ) + 0.4                        # cosχ*** (PF1 = 0.4)
    cosχ2 = cos(φ_pp - δ) - (2 / π) * φ_pp * sin(δ)    # cosχ**

    # F1 — local-time dependency (Eq. 4-7)
    V_D = 2π * (lt - 14) / 24
    V_SD = 2π * lt / 12
    V_TD = 2π * lt / 8
    F1 =
        cosχ3 +
        cosχ2 * (
            k[1] * cos(V_D) +
            k[2] * cos(V_SD) +
            k[3] * sin(V_SD) +
            k[4] * cos(V_TD) +
            k[5] * sin(V_TD)
        )

    # F2 — seasonal dependency (Eq. 8-10)
    V_A = 2π * (day_of_year - 18) / 365.25
    V_SA = 4π * (day_of_year - 6) / 365.25
    F2 = 1 + k[6] * cos(V_A) + k[7] * cos(V_SA)

    # Geomagnetic latitude of the pierce point (Eq. 31)
    φ_m = asin(
        sin(φ_pp) * sin(_NTCM_G_GNP_LAT) +
        cos(φ_pp) * cos(_NTCM_G_GNP_LAT) * cos(λ_pp - _NTCM_G_GNP_LON),
    )

    # F3 — geomagnetic field dependency (Eq. 11), φ_m in radians
    F3 = 1 + k[8] * cos(φ_m)

    # F4 — equatorial (Appleton) anomaly dependency (Eq. 12-14), φ_m in degrees
    φ_m_deg = rad2deg(φ_m)
    EC1 = -(φ_m_deg - 16.0)^2 / (2 * 12.0^2)
    EC2 = -(φ_m_deg + 10.0)^2 / (2 * 13.0^2)
    F4 = 1 + k[9] * exp(EC1) + k[10] * exp(EC2)

    # F5 — solar activity dependency (Eq. 15)
    F5 = k[11] + k[12] * azpar

    return max(F1 * F2 * F3 * F4 * F5, 0.0)
end

"""
    ntcm_g_stec(elevation, azimuth, lla, day_of_year, ut, a_i0, a_i1, a_i2) -> Float64

Slant TEC in TECU along the user→satellite line of sight using NTCM-G, for the
satellite `elevation`/`azimuth` (radians) seen from the user geodetic position
`lla` (a `Geodesy.LLA`), on `day_of_year`, at universal time `ut` (hours), and the
broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`.
"""
function ntcm_g_stec(elevation, azimuth, lla, day_of_year, ut, a_i0, a_i1, a_i2)
    φ_pp, λ_pp =
        _ntcm_g_pierce_point(deg2rad(lla.lat), deg2rad(lla.lon), elevation, azimuth)
    vtec = ntcm_g_vtec(φ_pp, λ_pp, day_of_year, ut, _ntcm_g_azpar(a_i0, a_i1, a_i2))
    return _ntcm_g_mapping_function(elevation) * vtec
end

# Day of year and universal time (hours) from a Galileo System Time (GST) week
# number and time of week [s], read off the TAI calendar label of the epoch, dated from
# `week_start_epoch` as `day_of_year` dates the tropospheric one. The label runs ahead
# of UTC by the accrued leap seconds (37 s), 0.01 h of UT, negligible for the day-of-year
# and local-time inputs of NTCM-G.
function _galileo_doy_and_ut(week, time_of_week)
    t = DateTime(week_start_epoch(GST(), week) + time_of_week)
    return dayofyear(t), hour(t) + minute(t) / 60 + second(t) / 3600 + millisecond(t) / 3.6e6
end
