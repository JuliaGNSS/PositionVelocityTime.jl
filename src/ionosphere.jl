# ===========================================================================
#  Ionospheric correction — constellation-wide model selection
#
#  The broadcast ionospheric coefficients are global to a GNSS (every GPS
#  satellite broadcasts the same Klobuchar α/β; every Galileo satellite the same
#  NTCM-G Effective Ionisation Level a_i0..a_i2). So rather than choosing a model
#  per satellite, a single model is chosen for the whole solve from whatever the
#  healthy decoders have delivered, and applied to *all* satellites:
#
#    - No coefficients decoded            → no correction (0 m for every sat).
#    - Only Klobuchar (GPS) decoded       → Klobuchar for every sat.
#    - Only NTCM-G (Galileo) decoded      → NTCM-G for every sat.
#    - Both decoded                       → NTCM-G for every sat (more accurate).
#
#  Only data actually decoded from the navigation message is used; there are no
#  user-supplied fallback coefficients. The ionospheric delay scales as 1/f², so a
#  model is evaluated for the constellation it came from and then rescaled to each
#  satellite's actual carrier frequency. This makes a single set of coefficients
#  applicable across constellations and frequency bands (L1/E1, L2, L5/E5a, E6, …),
#  not just the band it was broadcast on.
# ===========================================================================

# The frequency the Klobuchar group delay refers to (GPS L1), taken from GNSSSignals
# so it stays consistent with the per-satellite carrier frequencies used below. Kept
# as a Hz quantity: only the (f_L1/f) ratio is used, which is dimensionless.
const _GPS_L1_FREQUENCY = GNSSSignals.get_center_frequency(GPSL1CA())

"""
    KlobucharParams(α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3)

The eight Klobuchar ionospheric coefficients decoded from a GPS L1 navigation
message, in IS-GPS-200 SI units (seconds and seconds·semicircle⁻ⁿ). The field
names mirror `GNSSDecoder.GPSL1CAData` (`α_0…α_3`, `β_0…β_3`).
"""
struct KlobucharParams
    α_0::Float64
    α_1::Float64
    α_2::Float64
    α_3::Float64
    β_0::Float64
    β_1::Float64
    β_2::Float64
    β_3::Float64
end

"""
    NTCMGParams(a_i0, a_i1, a_i2, week_number::Integer)

The broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`
(sfu, sfu/deg, sfu/deg²) decoded from an E1B navigation message, together with the
GST week number needed to derive the day of year and universal time for NTCM-G.
"""
struct NTCMGParams
    a_i0::Float64
    a_i1::Float64
    a_i2::Float64
    week_number::Int
end

"""
    klobuchar_params(decoder) -> Union{KlobucharParams,Nothing}

Klobuchar α/β decoded from a GPS navigation message (LNAV `GPSL1CAData`, CNAV
`GPSCNAVData` on L5/L2C, or CNAV-2 `GPSL1C_DData`), or `nothing` if they have not
been broadcast yet or the decoder is not a GPS signal. The same single-frequency
Klobuchar model is broadcast on all GPS civil signals.
"""
klobuchar_params(decoder) = nothing
function klobuchar_params(decoder::GNSSDecoder.GNSSDecoderState{<:GPSNavData})
    d = decoder.data
    # All eight coefficients must be present: they are decoded together (subframe 4,
    # page 18), but guard each so a partially-populated decoder returns `nothing`
    # rather than throwing when a `nothing` hits a `Float64` field.
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return KlobucharParams(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

"""
    ntcm_g_params(decoder) -> Union{NTCMGParams,Nothing}

NTCM-G Effective Ionisation Level coefficients decoded from a Galileo navigation
message (I/NAV on E1B or F/NAV on E5a — both broadcast `a_i0…a_i2`), or `nothing`
if they (or the week number) have not been decoded yet or the decoder is not Galileo.
"""
ntcm_g_params(decoder) = nothing
function ntcm_g_params(decoder::GNSSDecoder.GNSSDecoderState{<:GalileoNavData})
    d = decoder.data
    (isnothing(d.a_i0) || isnothing(d.a_i1) || isnothing(d.a_i2) || isnothing(d.WN)) &&
        return nothing
    return NTCMGParams(d.a_i0, d.a_i1, d.a_i2, d.WN)
end

"""
    select_ionospheric_correction(states) -> Union{KlobucharParams,NTCMGParams,Nothing}

Scan all (healthy) satellite decoders and pick the single ionospheric correction
to apply to the whole solve. NTCM-G is preferred whenever Galileo coefficients are
available (it is the more accurate model); otherwise Klobuchar is used if GPS
coefficients are available; otherwise `nothing` (no correction). The coefficients
are global to a constellation, so the first decoder that carries each set is used.
"""
function select_ionospheric_correction(states)
    klobuchar = nothing
    ntcm_g = nothing
    for state in states
        klobuchar === nothing && (klobuchar = klobuchar_params(state.decoder))
        ntcm_g === nothing && (ntcm_g = ntcm_g_params(state.decoder))
    end
    ntcm_g !== nothing && return ntcm_g   # most accurate model when available
    return klobuchar                      # KlobucharParams, or nothing if neither
end

"""
    ionospheric_delay(correction, system, elevation, azimuth, lla, gps_time) -> Float64

Slant ionospheric group delay in metres for one satellite (`system` is the
satellite's GNSS, used for its carrier frequency), using the constellation-wide
`correction` returned by [`select_ionospheric_correction`](@ref):

- `::Nothing` → `0.0` (no coefficients were decoded).
- [`KlobucharParams`](@ref) → Klobuchar model (IS-GPS-200, Fig. 20-4).
- [`NTCMGParams`](@ref) → NTCM-G model.

The line of sight is given by the satellite `elevation`/`azimuth` (radians) and
the user geodetic position `lla` (a `Geodesy.LLA`); `gps_time` is the system time
of week of the measurement in seconds. The geometry is taken precomputed — and
shared across satellites and with [`tropospheric_delay`](@ref) — so a whole-epoch
correction does the user geodetic conversion only once. Derive the geometry from
ECEF with `LLAfromECEF(wgs84)(user)` and
[`_elevation_azimuth`](@ref)`(ENUfromECEF(user, wgs84), sat)`.
"""
ionospheric_delay(::Nothing, system, elevation, azimuth, lla, gps_time) = 0.0

function ionospheric_delay(p::KlobucharParams, system, elevation, azimuth, lla, gps_time)
    # IS-GPS-200 works in semicircles: lat/lon in deg/180, elevation/azimuth in rad/π.
    l1_seconds = klobuchar_group_delay(
        lla.lat / 180,
        lla.lon / 180,
        elevation / π,
        azimuth / π,
        gps_time,
        (p.α_0, p.α_1, p.α_2, p.α_3),
        (p.β_0, p.β_1, p.β_2, p.β_3),
    )
    # The Klobuchar broadcast coefficients define the group delay at the GPS L1
    # frequency (IS-GPS-200). The ionospheric delay scales as 1/f², so rescale it
    # to this satellite's actual carrier frequency — yielding the correct (larger)
    # delay on lower bands such as L2, L5/E5a and E6. For an L1/E1 signal the
    # factor is exactly 1. The ratio of two Hz quantities is dimensionless.
    f = GNSSSignals.get_center_frequency(system)
    return SPEEDOFLIGHT * l1_seconds * (_GPS_L1_FREQUENCY / f)^2
end

function ionospheric_delay(p::NTCMGParams, system, elevation, azimuth, lla, gps_time)
    doy, ut = _galileo_doy_and_ut(p.week_number, gps_time)
    stec = ntcm_g_stec(elevation, azimuth, lla, doy, ut, p.a_i0, p.a_i1, p.a_i2) # TECU
    f = ustrip(Hz, GNSSSignals.get_center_frequency(system))
    # Eq. 1: group delay [m] = 40.3 / f² · STEC, with STEC in electrons/m² (1 TECU = 1e16).
    return 40.3 / f^2 * stec * 1.0e16
end

"""
    klobuchar_group_delay(φ_u, λ_u, E, A, gps_time, α, β)

Klobuchar single-frequency ionospheric group delay for GPS L1 (IS-GPS-200N,
Fig. 20-4), returned in **seconds**. All angles are in **semicircles**:
`φ_u`/`λ_u` are the user geodetic latitude/longitude, `E`/`A` the satellite
elevation/azimuth. `gps_time` is GPS system time in seconds. `α`/`β` are the
4-element Klobuchar coefficient tuples (SI units).
"""
function klobuchar_group_delay(φ_u, λ_u, E, A, gps_time, α, β)
    # Earth-centred angle between user and ionospheric pierce point (semicircles)
    ψ = 0.0137 / (E + 0.11) - 0.022
    # Geodetic latitude of the ionospheric pierce point (IPP), clamped per ICD
    φ_i = clamp(φ_u + ψ * cos(A * π), -0.416, 0.416)
    # Geodetic longitude of the IPP
    λ_i = λ_u + ψ * sin(A * π) / cos(φ_i * π)
    # Geomagnetic latitude of the IPP
    φ_m = φ_i + 0.064 * cos((λ_i - 1.617) * π)
    # Local time at the IPP (seconds), wrapped to [0, 86400)
    t = mod(4.32e4 * λ_i + gps_time, 86400.0)
    # Obliquity / slant factor
    F = 1.0 + 16.0 * (0.53 - E)^3
    # Amplitude (s) and period (s) of the cosine model, with ICD floors
    AMP = max(α[1] + φ_m * (α[2] + φ_m * (α[3] + φ_m * α[4])), 0.0)
    PER = max(β[1] + φ_m * (β[2] + φ_m * (β[3] + φ_m * β[4])), 72000.0)
    x = 2π * (t - 50400.0) / PER
    return abs(x) < 1.57 ? F * (5.0e-9 + AMP * (1 - x^2 / 2 + x^4 / 24)) : F * 5.0e-9
end

"""
    _elevation_azimuth(enu_from_ecef::ENUfromECEF, sat_position) -> (elevation, azimuth)

Elevation and azimuth (radians) of `sat_position` (ECEF) in the local East-North-Up
frame defined by `enu_from_ecef = ENUfromECEF(user_position, wgs84)`. Azimuth is
measured clockwise from North. The transform is taken precomputed so it can be built
once per user position and reused across satellites.
"""
function _elevation_azimuth(enu_from_ecef::ENUfromECEF, sat_position)
    sat_enu = get_sat_enu(enu_from_ecef, ECEF(sat_position))
    elevation = sat_enu.ϕ
    # `SphericalFromCartesian` measures θ counter-clockwise from East (the ENU +x
    # axis); the ionospheric models use azimuth measured clockwise from North, i.e.
    # π/2 − θ. Only cos/sin of the azimuth are used downstream, so the wrap is moot.
    azimuth = π / 2 - sat_enu.θ
    return elevation, azimuth
end

# ===========================================================================
#  NTCM-G — Galileo single-frequency ionospheric model
#  (European GNSS (Galileo) NTCM-G Ionospheric Model Description, Issue 1.0,
#   May 2022). Driven by the broadcast Effective Ionisation Level coefficients
#   a_i0, a_i1, a_i2. Returns slant TEC (TECU); equation numbers below refer to
#   that document.
# ===========================================================================

# NTCM-G model coefficients k1..k12 (Table 3)
const _NTCM_K = (
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
const _NTCM_RE = 6371.0                 # Earth mean radius [km] (Table 2)
const _NTCM_HI = 450.0                  # ionospheric pierce point height [km]
const _NTCM_GNP_LAT = deg2rad(79.74)    # geomagnetic North pole latitude
const _NTCM_GNP_LON = deg2rad(-71.78)   # geomagnetic North pole longitude

# Effective Ionisation Level Azpar [sfu] from the broadcast coefficients (Eq. 2).
function _azpar(a_i0, a_i1, a_i2)
    radicand = a_i0^2 + 1633.33 * a_i1^2 + 4802000.0 * a_i2^2 + 3266.67 * a_i0 * a_i2
    return sqrt(max(radicand, 0.0))
end

# Ionospheric pierce point geographic latitude/longitude [rad] (Eq. 24-26),
# given user geodetic lat/lon and satellite elevation/azimuth [rad].
function _pierce_point(φ_u, λ_u, elevation, azimuth)
    ψ = π / 2 - elevation - asin(_NTCM_RE / (_NTCM_RE + _NTCM_HI) * cos(elevation))
    # Clamp to asin's domain: the φ_pp argument is a unit dot-product that can
    # overshoot ±1 by a rounding ulp, and the λ_pp argument genuinely diverges as
    # the pierce point approaches a pole (cos(φ_pp) → 0).
    φ_pp = asin(clamp(sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth), -1.0, 1.0))
    λ_pp = λ_u + asin(clamp(sin(ψ) * sin(azimuth) / cos(φ_pp), -1.0, 1.0))
    return φ_pp, λ_pp
end

# Sun's declination [rad] for the day of year (Eq. 28).
_sun_declination(doy) = deg2rad(23.44) * sin(deg2rad(0.9856 * (doy - 80.7)))

# Modified Single Layer Model mapping function (Eq. 32-33).
function _mslm_mapping_function(elevation)
    sinz = _NTCM_RE / (_NTCM_RE + _NTCM_HI) * sin(0.9782 * (π / 2 - elevation))
    return 1.0 / sqrt(1.0 - sinz^2)
end

"""
    ntcm_g_vtec(φ_pp, λ_pp, doy, ut, azpar) -> Float64

Vertical TEC in TECU at the ionospheric pierce point (geographic latitude/longitude
`φ_pp`/`λ_pp` in radians) for day of year `doy`, universal time `ut` in hours, and
Effective Ionisation Level `azpar` in solar flux units. Implements the NTCM-G model
`VTEC = F1·F2·F3·F4·F5` (Eq. 3-15).
"""
function ntcm_g_vtec(φ_pp, λ_pp, doy, ut, azpar)
    k = _NTCM_K
    lt = ut + rad2deg(λ_pp) / 15                      # local time [h] (Eq. 27)
    δ = _sun_declination(doy)

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
    V_A = 2π * (doy - 18) / 365.25
    V_SA = 4π * (doy - 6) / 365.25
    F2 = 1 + k[6] * cos(V_A) + k[7] * cos(V_SA)

    # Geomagnetic latitude of the pierce point (Eq. 31)
    φ_m = asin(
        sin(φ_pp) * sin(_NTCM_GNP_LAT) +
        cos(φ_pp) * cos(_NTCM_GNP_LAT) * cos(λ_pp - _NTCM_GNP_LON),
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
    ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2) -> Float64

Slant TEC in TECU along the user→satellite line of sight using NTCM-G, for the
satellite `elevation`/`azimuth` (radians) seen from the user geodetic position
`lla` (a `Geodesy.LLA`), day of year `doy`, universal time `ut` (hours), and the
broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`.
"""
function ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2)
    φ_pp, λ_pp = _pierce_point(deg2rad(lla.lat), deg2rad(lla.lon), elevation, azimuth)
    vtec = ntcm_g_vtec(φ_pp, λ_pp, doy, ut, _azpar(a_i0, a_i1, a_i2))
    return _mslm_mapping_function(elevation) * vtec
end

# Day of year and universal time (hours) from a Galileo System Time (GST) week
# number and time of week [s]. GST epoch is 1999-08-22 00:00:00 UTC; GST is
# continuous (offset from UTC by leap seconds, ~18 s), which is negligible for
# the day-of-year / UT inputs of NTCM-G.
function _galileo_doy_and_ut(week_number, time_of_week)
    total_seconds = week_number * 604800 + time_of_week
    days = floor(Int, total_seconds / 86400)
    ut = (total_seconds - days * 86400) / 3600
    return dayofyear(Date(1999, 8, 22) + Day(days)), ut
end
