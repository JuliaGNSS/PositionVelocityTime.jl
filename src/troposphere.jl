# ===========================================================================
#  Tropospheric correction — Saastamoinen zenith delays × Niell mapping functions
#
#  The tropospheric delay is non-dispersive — the same on every GNSS frequency.
#  Unlike the ionosphere there are no broadcast coefficients: the delay is a
#  "blind" function of the user position and a standard atmosphere, mapped to the
#  line of sight by the satellite elevation. The model has two independent parts:
#
#   1. Zenith delays: the Saastamoinen hydrostatic and wet zenith delays driven by
#      a standard atmosphere, exactly as in RTKLIB's `tropmodel` (which GNSS-SDR
#      and PocketSDR reuse for single-point positioning), with the relative
#      humidity defaulting to 70 %.
#   2. Obliquity mapping: the Niell mapping functions (Niell 1996, as in RTKLIB's
#      `tropmapf`/`nmf` used for PPP, and Orekit's `NiellMappingFunctionModel`) —
#      a continued fraction in sin(elevation) with separate hydrostatic and wet
#      coefficients interpolated over latitude, a seasonal term, and a height
#      correction. Unlike the flat-slab 1/cos(z) mapping RTKLIB's `tropmodel`
#      itself uses, it accounts for Earth curvature and ray bending, so it stays
#      finite at the horizon and does not over-predict at low elevation
#      (+14 % ≈ +3 m at 5° for 1/cos(z); see #62).
#
#  Saastamoinen's own −B·tan²(z) + δR slant correction (carried by Orekit's
#  `ModifiedSaastamoinenModel`) is deliberately NOT included: it is the 1972
#  curvature/bending fix *to the 1/cos(z) mapping*, superseded by the mapping
#  function — adding both would double-count, and the tan² term inverts the sign
#  of the delay below ≈ 1.9° elevation, forcing Orekit's low-elevation clamp.
#
#  Reference: A. E. Niell (1996), "Global mapping functions for the atmosphere
#  delay at radio wavelengths", J. Geophys. Res. 101(B2), 3227–3246. The original
#  paper has typos in eqs. (4)–(5); the coefficients below are the corrected
#  Table 3/4 values, identical in RTKLIB (`rtkcmn.c`) and Orekit.
# ===========================================================================

const _DEFAULT_RELATIVE_HUMIDITY = 0.7  # GNSS-SDR / RTKLIB default

# Niell hydrostatic coefficients: averages and seasonal amplitudes at latitudes
# 15°, 30°, 45°, 60°, 75° (linearly interpolated in |latitude|, constant outside).
# The coefficient at day-of-year `doy` is `average − amplitude·cos(2π(doy−28)/365.25)`,
# phase-shifted half a year in the southern hemisphere.
const _NMF_HYDROSTATIC_AVERAGE = (
    a = (1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3),
    b = (2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3),
    c = (62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3),
)
const _NMF_HYDROSTATIC_AMPLITUDE = (
    a = (0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5),
    b = (0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5),
    c = (0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5),
)
# Niell wet coefficients: latitude-dependent only — the wet delay's short
# correlation time makes a seasonal average meaningless (Niell 1996, §4).
const _NMF_WET = (
    a = (5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4),
    b = (1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3),
    c = (4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2),
)
# Height-correction coefficients (Niell 1996 eq. 6): the hydrostatic mapping
# refers to sea level, and the correction `(1/sin E − m_ht(E)) · height[km]`
# accounts for the atmosphere column below a raised site.
const _NMF_HEIGHT = (a = 2.53e-5, b = 5.49e-3, c = 1.14e-3)

# Lowest elevation the *height-correction term* is evaluated at. Its leading
# 1/sin(E) is the same flat-slab form the continued fractions replace, so the term
# diverges at the horizon (and is 0·∞ = NaN there at height 0) even though the
# mapping functions themselves stay finite; below this bound the term is frozen at
# its value here. The freeze is harmless — the term is ≈ 7 cm per km of site height
# at the bound — and the value coincides with the ≈ 3° lower validity limit of the
# mapping-function literature (Niell 1996, VMF3, IERS Conventions 2010 ch. 9) and
# with Orekit's `DEFAULT_LOW_ELEVATION_THRESHOLD`.
const _HEIGHT_CORRECTION_MIN_ELEVATION = 0.05  # rad ≈ 2.87°

# Piecewise-linear interpolation of a Niell coefficient over |latitude| in degrees,
# with nodes at 15°, 30°, 45°, 60°, 75° and held constant outside (RTKLIB `interpc`).
function _interpolate_nmf_coefficient(values::NTuple{5,Float64}, abs_latitude_deg)
    x = abs_latitude_deg / 15.0
    i = floor(Int, x)
    i < 1 && return values[1]
    i >= 5 && return values[5]
    w = x - i
    return values[i] * (1.0 - w) + values[i+1] * w
end

# Normalized Marini continued fraction (Herring normalization): exactly 1 at
# zenith, finite at the horizon (sin E = 0), strictly decreasing in elevation.
function _marini_mapping(sin_el, a, b, c)
    return (1.0 + a / (1.0 + b / (1.0 + c))) /
           (sin_el + a / (sin_el + b / (sin_el + c)))
end

"""
    niell_mapping_functions(latitude, height, elevation, doy) -> (hydrostatic, wet)

The Niell (1996) hydrostatic and wet mapping functions: dimensionless factors that
map the zenith tropospheric delays to a satellite at `elevation` (radians).
`latitude` is in radians, `height` is the site's ellipsoidal height in metres (it
enters the hydrostatic height correction; RTKLIB likewise substitutes the
ellipsoidal height for the orthometric one), and `doy` is the day of year driving
the hydrostatic seasonal term (phase-shifted half a year in the southern
hemisphere). Elevations below the horizon are evaluated at the horizon, where both
factors are finite (≈ 37 hydrostatic, ≈ 57 wet); both decrease strictly to exactly
1 at zenith.
"""
function niell_mapping_functions(latitude, height, elevation, doy)
    el = max(elevation, 0.0)
    sin_el = sin(el)
    abs_lat = abs(rad2deg(latitude))
    # Seasonal term, southern hemisphere shifted by half a year
    cosy = cos(2π * ((doy - 28.0) / 365.25 + (latitude < 0.0 ? 0.5 : 0.0)))
    ah = _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AVERAGE.a, abs_lat) -
         _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AMPLITUDE.a, abs_lat) * cosy
    bh = _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AVERAGE.b, abs_lat) -
         _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AMPLITUDE.b, abs_lat) * cosy
    ch = _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AVERAGE.c, abs_lat) -
         _interpolate_nmf_coefficient(_NMF_HYDROSTATIC_AMPLITUDE.c, abs_lat) * cosy
    aw = _interpolate_nmf_coefficient(_NMF_WET.a, abs_lat)
    bw = _interpolate_nmf_coefficient(_NMF_WET.b, abs_lat)
    cw = _interpolate_nmf_coefficient(_NMF_WET.c, abs_lat)
    el_ht = max(el, _HEIGHT_CORRECTION_MIN_ELEVATION)
    dm =
        (1.0 / sin(el_ht) -
         _marini_mapping(sin(el_ht), _NMF_HEIGHT.a, _NMF_HEIGHT.b, _NMF_HEIGHT.c)) *
        height / 1.0e3
    hydrostatic = _marini_mapping(sin_el, ah, bh, ch) + dm
    wet = _marini_mapping(sin_el, aw, bw, cw)
    return hydrostatic, wet
end

"""
    saastamoinen_zenith_delays(latitude, height, humidity) -> (hydrostatic, wet)

Saastamoinen zenith hydrostatic and wet tropospheric delays in metres, from a
standard atmosphere (15 °C at sea level) evaluated at the geodetic `latitude`
(radians) and `height` (metres, clamped below at 0). `humidity` is the relative
humidity (0…1) feeding the water-vapour partial pressure of the wet term. The
formulas are those of RTKLIB's `tropmodel` at zenith.
"""
function saastamoinen_zenith_delays(latitude, height, humidity)
    hgt = max(height, 0.0)
    # Standard atmosphere
    pressure = 1013.25 * (1.0 - 2.2557e-5 * hgt)^5.2568          # hPa
    temperature = 15.0 - 6.5e-3 * hgt + 273.16                  # K (15 °C at sea level)
    e = 6.108 * humidity * exp((17.15 * temperature - 4684.0) / (temperature - 38.45))  # hPa
    hydrostatic =
        0.0022768 * pressure /
        (1.0 - 0.00266 * cos(2.0 * latitude) - 0.00028 * hgt / 1.0e3)
    wet = 0.002277 * (1255.0 / temperature + 0.05) * e
    return hydrostatic, wet
end

"""
    tropospheric_delay(elevation, lla, doy; humidity = 0.7) -> Float64

Slant tropospheric group delay in metres for a satellite at `elevation` (radians)
seen from the user geodetic position `lla` (a `Geodesy.LLA`): the Saastamoinen
zenith delays (see [`saastamoinen_zenith_delays`](@ref)) mapped to the line of
sight by the Niell mapping functions (see [`niell_mapping_functions`](@ref)),
whose seasonal term is driven by the day of year `doy`. The delay is
non-dispersive (frequency independent). `humidity` is the relative humidity (0…1)
used for the wet component. The delay is finite and monotone in elevation all the
way down to the horizon (≈ 90 m for a grazing ray at sea level); elevations below
the horizon reuse the horizon value. Returns `0.0` when the user height is outside
the model's valid range (−100 m … 10 km).

The geometry is taken precomputed so a whole-epoch correction shares the user
geodetic conversion and the elevation with [`ionospheric_delay`](@ref).
"""
function tropospheric_delay(elevation, lla, doy; humidity = _DEFAULT_RELATIVE_HUMIDITY)
    (lla.alt < -100.0 || lla.alt > 1.0e4) && return 0.0
    latitude = deg2rad(lla.lat)
    height = max(lla.alt, 0.0)
    zenith_hydrostatic, zenith_wet = saastamoinen_zenith_delays(latitude, height, humidity)
    hydrostatic_mapping, wet_mapping =
        niell_mapping_functions(latitude, height, elevation, doy)
    return zenith_hydrostatic * hydrostatic_mapping + zenith_wet * wet_mapping
end

# Day of year (1–366) of the GNSS system time given by the absolute `week` and
# `time_of_week` [s] of `system`. The system time scale's offset from UTC (leap
# seconds, ≤ ~18 s) is negligible for a seasonal argument with a one-year period.
function _day_of_year(system, week, time_of_week)
    epoch = get_system_start_time(system)
    t = epoch + Millisecond(round(Int, (week * 604800 + time_of_week) * 1000))
    return dayofyear(t)
end
