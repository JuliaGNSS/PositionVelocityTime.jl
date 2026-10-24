# BDGIM, the BeiDou Global Ionospheric delay correction Model (BDS-SIS-ICD-B1C-1.0 §7.8;
# the same in BDS-SIS-ICD-B2a-1.0 §7.8 and, with equation and table numbers one lower,
# BDS-SIS-ICD-B2b-1.0 §7.7). Equation numbers are the B1C ones. Nine broadcast α_i weight
# normalised harmonics of the pierce point's solar-fixed geomagnetic coordinates, on top
# of a predictive part A_0 of 17 harmonics whose amplitudes are Fourier series fixed in
# the ICD. Returns TECU.

# Suggested parameter values, BDS-SIS-ICD-B1C-1.0 §7.8.2 (the paragraph after
# Eq. 7-17). Re and H_ion in km — only their ratio enters, so the unit is free.
const _BDGIM_RE = 6378.0                    # mean radius of the Earth [km]
const _BDGIM_HION = 400.0                   # ionospheric single-layer shell altitude [km]
const _BDGIM_POLE_LAT = 80.27 / 180 * π     # north magnetic pole geographic latitude [rad]
const _BDGIM_POLE_LON = -72.58 / 180 * π    # north magnetic pole geographic longitude [rad]

# Degree/order (n_i, m_i) of the nine broadcast harmonics, Table 7-11. A negative
# m selects the sine rather than the cosine term (Eq. 7-11).
const _BDGIM_NM =
    ((0, 0), (1, 0), (1, 1), (1, -1), (2, 0), (2, 1), (2, -1), (2, 2), (2, -2))

# Degree/order (n_j, m_j) of the seventeen predictive harmonics, from Table 7-12's header
# row (Eq. 7-14's "Table 7-11" is an ICD slip: that table lists only the broadcast pairs).
const _BDGIM_NM_PREDICTED = (
    (3, 0), (3, 1), (3, -1), (3, 2), (3, -2), (3, 3), (3, -3),   # degree 3
    (4, 0), (4, 1), (4, -1), (4, 2), (4, -2),                    # degree 4
    (5, 0), (5, 1), (5, -1), (5, 2), (5, -2),                    # degree 5
)

# Prediction periods T_k [days], k = 1…12, from the rightmost column of Table 7-12.
const _BDGIM_PERIODS =
    (1.0, 0.5, 0.33, 14.6, 27.0, 121.6, 182.51, 365.25, 4028.71, 2014.35, 1342.90, 1007.18)

# Table 7-12, the non-broadcast coefficients [TECU]: the constant row a_{0,j} and the
# twelve (a_{k,j}, b_{k,j}) pairs, 17 entries each.
#
# Transcribed three ways, which agreed entry for entry: the text layer of
# BDS-SIS-ICD-B1C-1.0 Table 7-12, the rendered page image of it (printed page 50), and the
# text layer of the independently typeset BDS-SIS-ICD-B2b-1.0 (2020) Table 7-10.
# `test/atmosphere/bdgim.jl` pins row, column and absolute-sum checksums against later typos.
const _BDGIM_A0 = (
    -0.61, -1.31, -2.0, -0.03, 0.15, -0.48, -0.4, 2.28, -0.16, -0.21, -0.1, -0.13, 0.21,
    0.68, 1.06, 0.0, -0.12,
)

const _BDGIM_A = (
    (-0.51, -0.43, 0.34, -0.01, 0.17, 0.02, -0.06, 0.3, 0.44, -0.28, -0.31, -0.17, 0.04,
        0.39, -0.12, 0.12, 0.0),
    (-0.06, -0.05, 0.06, 0.17, 0.15, 0.0, 0.11, -0.05, -0.16, 0.02, 0.11, 0.04, 0.12, 0.07,
        0.02, -0.14, -0.14),
    (0.01, -0.03, 0.01, -0.01, 0.05, -0.03, 0.05, -0.03, -0.01, 0.0, -0.08, -0.04, 0.0,
        -0.02, -0.03, 0.0, -0.03),
    (-0.01, 0.0, 0.01, 0.0, 0.01, 0.0, -0.01, -0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (0.0, 0.0, 0.03, 0.01, 0.02, 0.01, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (-0.19, -0.02, 0.12, -0.1, 0.06, 0.0, -0.02, -0.08, -0.02, -0.07, 0.01, 0.03, 0.15,
        0.06, -0.05, -0.03, -0.1),
    (-0.18, 0.06, -0.55, -0.02, 0.09, -0.08, 0.0, 0.86, -0.18, -0.05, -0.07, 0.04, 0.14,
        -0.03, 0.37, -0.11, -0.12),
    (1.09, -0.14, -0.21, 0.52, 0.27, 0.0, 0.11, 0.17, 0.23, 0.35, -0.05, 0.02, -0.6, 0.02,
        0.01, 0.27, 0.32),
    (-0.34, -0.09, -1.22, 0.05, 0.15, -0.29, -0.17, 1.58, -0.06, -0.15, 0.0, 0.13, 0.28,
        -0.08, 0.62, -0.01, -0.04),
    (-0.13, 0.07, -0.37, 0.05, 0.06, -0.11, -0.07, 0.46, 0.0, -0.04, 0.01, 0.07, 0.09,
        -0.05, 0.15, -0.01, 0.01),
    (-0.06, 0.13, -0.07, 0.03, 0.02, -0.05, -0.05, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
    (-0.03, 0.08, -0.01, 0.04, 0.01, -0.02, -0.02, -0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
)

const _BDGIM_B = (
    (0.23, -0.2, -0.31, 0.16, -0.03, 0.02, 0.04, 0.18, 0.34, 0.45, 0.19, -0.25, -0.12, 0.18,
        0.4, -0.09, 0.21),
    (0.02, -0.08, -0.06, -0.11, 0.15, -0.14, 0.01, 0.01, 0.04, -0.14, -0.05, 0.08, 0.08,
        -0.01, 0.01, 0.11, -0.12),
    (0.0, -0.02, -0.03, -0.05, -0.01, -0.07, -0.03, -0.01, 0.02, -0.01, 0.03, -0.1, 0.01,
        0.05, -0.01, 0.04, 0.0),
    (0.0, -0.02, 0.01, 0.0, -0.01, 0.01, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (0.01, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (-0.09, 0.07, 0.03, 0.06, 0.09, 0.01, 0.02, 0.0, -0.04, -0.02, -0.01, 0.01, -0.1, 0.0,
        -0.01, 0.02, 0.05),
    (0.15, -0.31, 0.13, 0.05, -0.09, -0.03, 0.06, -0.36, 0.08, 0.05, 0.06, -0.02, -0.05,
        0.06, -0.2, 0.04, 0.07),
    (0.5, -0.08, -0.38, 0.36, 0.14, 0.04, 0.0, 0.25, 0.17, 0.27, -0.03, -0.03, -0.32, -0.1,
        0.2, 0.1, 0.3),
    (0.0, -0.11, -0.22, 0.01, 0.02, -0.03, -0.01, 0.49, -0.03, -0.02, 0.01, 0.02, 0.04,
        -0.04, 0.16, -0.02, -0.01),
    (0.05, 0.03, 0.07, 0.02, -0.01, 0.03, 0.02, -0.04, -0.01, -0.01, 0.02, 0.03, 0.02,
        -0.04, -0.04, -0.01, 0.0),
    (0.03, -0.02, 0.04, -0.01, -0.03, 0.02, 0.01, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
    (0.04, -0.02, -0.04, 0.0, -0.01, 0.0, 0.01, 0.07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
)

"""
    _bdgim_legendre(n, m, x) -> Float64

Un-normalised associated Legendre function `P_{n,m}(x)` (Eq. 7-13), by the ICD's
recursion: `P_{n,n} = (2n−1)!!·(1−x²)^(n/2)`, `P_{m+1,m} = x·(2m+1)·P_{m,m}`, and
`P_{n,m} = [(2n−1)·x·P_{n−1,m} − (n+m−1)·P_{n−2,m}] / (n−m)`, for `5 ≥ n ≥ m ≥ 0`.
"""
function _bdgim_legendre(n::Int, m::Int, x)
    # (2m−1)!!, an empty product at m = 0.
    double_factorial = 1.0
    for i in 1:m
        double_factorial *= 2i - 1
    end
    p_mm = double_factorial * (1 - x^2)^(m / 2)
    n == m && return p_mm
    p_prev = p_mm                       # P_{m,m}
    p_cur = x * (2m + 1) * p_mm         # P_{m+1,m}
    for k in (m+2):n
        p_prev, p_cur = p_cur, ((2k - 1) * x * p_cur - (k + m - 1) * p_prev) / (k - m)
    end
    return p_cur
end

# Normalisation N_{n,m} of Eq. 7-12: √[(n−m)!·(2n+1)·(2 − δ_{0,m}) / (n+m)!]. Degree ≤ 5
# keeps (n+m)! exact in `Int`.
_bdgim_normalisation(n::Int, m::Int) =
    sqrt(factorial(n - m) * (2n + 1) * (m == 0 ? 1 : 2) / factorial(n + m))

# One normalised spherical-harmonic term of Eq. 7-11 / Eq. 7-14, at the solar-fixed
# geomagnetic latitude/longitude φ′/λ′ [rad]: P̃_{|n|,|m|}(sin φ′) times cos(m·λ′)
# for m ≥ 0 and sin(−m·λ′) for m < 0.
function _bdgim_harmonic(n::Int, m::Int, sinφ′, λ′)
    p = _bdgim_normalisation(abs(n), abs(m)) * _bdgim_legendre(abs(n), abs(m), sinφ′)
    return m ≥ 0 ? p * cos(m * λ′) : p * sin(-m * λ′)
end

"""
    _bdgim_pierce_point(φ_u, λ_u, elevation, azimuth) -> (φ_g, λ_g)

Geographic latitude/longitude [rad] of the pierce point's Earth projection (Eq. 7-7 and
7-8), for a user at geodetic `φ_u`/`λ_u` seeing a satellite at `elevation`/`azimuth`,
all in radians.
"""
function _bdgim_pierce_point(φ_u, λ_u, elevation, azimuth)
    ψ = π / 2 - elevation - asin(_BDGIM_RE / (_BDGIM_RE + _BDGIM_HION) * cos(elevation))
    # Clamped against a rounding ulp past ±1.
    φ_g = asin(
        clamp(
            sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth),
            -1.0,
            1.0,
        ),
    )
    # Eq. 7-8 prints a one-argument arctan of a quotient, which folds the offset into
    # (−π/2, π/2) and is wrong when the pierce point is beyond the pole. The quotient's
    # terms are the great-circle pair (cos φ_g·sin Δλ, cos φ_g·cos Δλ) up to a positive
    # factor, so the two-argument `atan` of them gives the full circle.
    λ_g = λ_u + atan(
        sin(ψ) * sin(azimuth) * cos(φ_u),
        cos(ψ) - sin(φ_u) * sin(φ_g),
    )
    return φ_g, λ_g
end

"""
    _bdgim_solar_geomagnetic(φ_g, λ_g, mjd) -> (φ′, λ′)

Latitude/longitude [rad] of the pierce point in the solar-fixed geomagnetic frame (Eq. 7-9
and 7-10), from its geographic `φ_g`/`λ_g` [rad] and the Modified Julian Date `mjd`. The
pole is the ICD's (80.27°N, 72.58°W); the sun's mean longitude is
`S_lon = π·(1 − 2·(t − ⌊t⌋))`, 180°E at MJD midnight. Only the longitude is referred to
the sun: `φ′ = φ_m`.
"""
function _bdgim_solar_geomagnetic(φ_g, λ_g, mjd)
    Δλ = λ_g - _BDGIM_POLE_LON
    φ_m = asin(
        clamp(
            sin(_BDGIM_POLE_LAT) * sin(φ_g) +
            cos(_BDGIM_POLE_LAT) * cos(φ_g) * cos(Δλ),
            -1.0,
            1.0,
        ),
    )
    # The two-argument arctan, as in `_bdgim_pierce_point`.
    λ_m = atan(
        cos(φ_g) * sin(Δλ) * cos(_BDGIM_POLE_LAT),
        sin(_BDGIM_POLE_LAT) * sin(φ_m) - sin(φ_g),
    )
    s_lon = π * (1 - 2 * (mjd - floor(mjd)))
    # Minus the Eq. 7-9 longitude of the sun itself (φ_g = 0, λ_g = S_lon).
    λ′ = λ_m - atan(
        sin(s_lon - _BDGIM_POLE_LON),
        sin(_BDGIM_POLE_LAT) * cos(s_lon - _BDGIM_POLE_LON),
    )
    return φ_m, λ′
end

"""
    _bdgim_reference_epoch(mjd) -> Float64

The prediction reference epoch `t_p` (Eq. 7-15), in days (MJD): the odd hour of the
day nearest `mjd`, which is always on the same day. Only the amplitudes β_j step every two
hours; VTEC stays continuous, as λ′ follows the true epoch.
"""
function _bdgim_reference_epoch(mjd)
    day = floor(mjd)
    hour_of_day = (mjd - day) * 24
    return day + (2 * floor(hour_of_day / 2) + 1) / 24
end

"""
    _bdgim_prediction_amplitudes(t_p) -> NTuple{17,Float64}

The seventeen prediction amplitudes `β_j` (Eq. 7-15) at `t_p` (days, MJD):
`β_j = a_{0,j} + Σ_{k=1}^{12} [a_{k,j}·cos(ω_k·t_p) + b_{k,j}·sin(ω_k·t_p)]` with
`ω_k = 2π/T_k` and the coefficients of Table 7-12. The ~4028.71-day term carries the solar
cycle.
"""
function _bdgim_prediction_amplitudes(t_p)
    # The twelve (sin, cos) pairs are shared by all seventeen amplitudes.
    trig = ntuple(k -> sincos(2π / _BDGIM_PERIODS[k] * t_p), Val(12))
    return ntuple(Val(17)) do j
        β = _BDGIM_A0[j]
        for k in eachindex(trig)
            s, c = trig[k]
            β += _BDGIM_A[k][j] * c + _BDGIM_B[k][j] * s
        end
        β
    end
end

"""
    bdgim_vtec(φ_g, λ_g, mjd, α) -> Float64

Vertical TEC in TECU at the pierce point `φ_g`/`λ_g` (geographic, radians) at Modified
Julian Date `mjd`, from the nine broadcast coefficients `α` (TECU):
`VTEC = A_0 + Σ α_i·A_i` (Eq. 7-16, with Eq. 7-11 and 7-14/7-15).

Floored at zero, as `ntcm_g_vtec` is, although the ICD prescribes no floor: a negative
VTEC would enter [`ionospheric_delay`](@ref) as an advance. It is reachable: A_0 (degrees
3–5 only) is a ±10–15 TECU term with no level of its own, negative over much of the globe,
and a small broadcast set (an all-zero α is legal) does not lift it.
"""
function bdgim_vtec(φ_g, λ_g, mjd, α)
    φ′, λ′ = _bdgim_solar_geomagnetic(φ_g, λ_g, mjd)
    sinφ′ = sin(φ′)
    # A_0, the predictive part (Eq. 7-14).
    β = _bdgim_prediction_amplitudes(_bdgim_reference_epoch(mjd))
    a_0 = sum(
        β[j] * _bdgim_harmonic(_BDGIM_NM_PREDICTED[j]..., sinφ′, λ′) for
        j in eachindex(β)
    )
    # The broadcast part (Eq. 7-11).
    broadcast_part = sum(
        α[i] * _bdgim_harmonic(_BDGIM_NM[i]..., sinφ′, λ′) for i in eachindex(α)
    )
    return max(a_0 + broadcast_part, 0.0)
end

"""
    bdgim_mapping_function(elevation) -> Float64

BDGIM's single-layer obliquity factor `M_F` (Eq. 7-17) at `elevation` (radians):
`1/√(1 − (Re/(Re+H_ion)·cos E)²)`, from 1 at the zenith to ~2.95 at the horizon.
"""
bdgim_mapping_function(elevation) =
    1 / sqrt(1 - (_BDGIM_RE / (_BDGIM_RE + _BDGIM_HION) * cos(elevation))^2)

"""
    bdgim_stec(elevation, azimuth, lla, mjd, α) -> Float64

Slant TEC in TECU, `M_F · VTEC`, along the line of sight at `elevation`/`azimuth`
(radians) from `lla` (a `Geodesy.LLA`) at Modified Julian Date `mjd`, from the broadcast
`α`: steps (1)–(5) of §7.8.2, step (6) being [`ionospheric_delay`](@ref).
"""
function bdgim_stec(elevation, azimuth, lla, mjd, α)
    φ_g, λ_g = _bdgim_pierce_point(deg2rad(lla.lat), deg2rad(lla.lon), elevation, azimuth)
    return bdgim_mapping_function(elevation) * bdgim_vtec(φ_g, λ_g, mjd, α)
end

# Modified Julian Date of a BDT week and time of week [s], read off the TAI calendar label
# as `day_of_year` is; its 37 s lead on UTC is negligible for a date.
const _MJD_EPOCH = TAITime(DateTime(1858, 11, 17))   # MJD 0, as a TAI calendar label
_beidou_modified_julian_date(week, time_of_week) =
    (week_start_epoch(BDT(), week) + time_of_week - _MJD_EPOCH) / 86_400
