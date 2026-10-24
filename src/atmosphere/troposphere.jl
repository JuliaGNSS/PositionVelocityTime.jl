# Tropospheric correction: non-dispersive, and "blind" — a function of the user position
# and a standard atmosphere, with no broadcast coefficients. Two parts:
#
#  1. Zenith delays: Saastamoinen's hydrostatic and wet delays from a standard
#     atmosphere, as in RTKLIB's `tropmodel` (reused by GNSS-SDR and PocketSDR), at 70 %
#     relative humidity by default.
#  2. Mapping: the Niell functions (as in RTKLIB's `tropmapf`/`nmf` and Orekit's
#     `NiellMappingFunctionModel`). Unlike `tropmodel`'s flat-slab 1/cos(z), they account
#     for Earth curvature and ray bending, stay finite at the horizon and do not
#     over-predict low (1/cos(z) is +14 % ≈ +3 m at 5°; see #62).
#
# Saastamoinen's −B·tan²(z) + δR slant term (Orekit's `ModifiedSaastamoinenModel`) is
# left out: it corrects the 1/cos(z) mapping the Niell functions replace, so it would
# double-count, and it inverts the delay's sign below ≈ 1.9°.
#
# Reference: A. E. Niell (1996), "Global mapping functions for the atmosphere delay at
# radio wavelengths", J. Geophys. Res. 101(B2), 3227–3246. The coefficients are the
# corrected Table 3/4 values (the paper's eqs. (4)–(5) have typos), as in RTKLIB
# (`rtkcmn.c`) and Orekit.

const _DEFAULT_RELATIVE_HUMIDITY = 0.7  # GNSS-SDR / RTKLIB default

# Niell hydrostatic coefficients: averages and seasonal amplitudes at latitudes 15°, 30°,
# 45°, 60°, 75°. A coefficient is `average − amplitude·cos(2π(day_of_year−28)/365.25)`,
# half a year later in the southern hemisphere.
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
# Niell wet coefficients: latitude-dependent only (Niell 1996, §4).
const _NMF_WET = (
    a = (5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4),
    b = (1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3),
    c = (4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2),
)
# Height-correction coefficients (Niell 1996 eq. 6): the hydrostatic mapping refers to sea
# level, corrected by `(1/sin E − m_ht(E)) · height[km]`.
const _NMF_HEIGHT = (a = 2.53e-5, b = 5.49e-3, c = 1.14e-3)

# Lowest elevation the height correction is evaluated at: its 1/sin(E) diverges at the
# horizon (NaN at height 0), so below this it is frozen, ≈ 7 cm per km of height. The
# bound is the ≈ 3° validity limit of the mapping-function literature (Niell 1996, VMF3,
# IERS Conventions 2010 ch. 9) and Orekit's `DEFAULT_LOW_ELEVATION_THRESHOLD`.
const _HEIGHT_CORRECTION_MIN_ELEVATION = 0.05  # rad ≈ 2.87°

# Linear interpolation of a Niell coefficient over |latitude| (degrees), constant outside
# 15°…75° (RTKLIB `interpc`).
function _interpolate_nmf_coefficient(values::NTuple{5,Float64}, abs_latitude_deg)
    x = abs_latitude_deg / 15.0
    i = floor(Int, x)
    i < 1 && return values[1]
    i >= 5 && return values[5]
    w = x - i
    return values[i] * (1.0 - w) + values[i+1] * w
end

# Normalised Marini continued fraction (Herring): 1 at zenith, finite at the horizon.
function _marini_mapping(sin_el, a, b, c)
    return (1.0 + a / (1.0 + b / (1.0 + c))) /
           (sin_el + a / (sin_el + b / (sin_el + c)))
end

"""
    niell_mapping_functions(latitude, height, elevation, day_of_year) -> (hydrostatic, wet)

The Niell (1996) hydrostatic and wet mapping functions at `elevation` (radians), for
`latitude` (radians), the site's ellipsoidal `height` (metres; used for the orthometric
one, as in RTKLIB) and `day_of_year`. Below the horizon they are evaluated at it, where
both are finite (≈ 37 hydrostatic, ≈ 57 wet); both fall to 1 at zenith.
"""
function niell_mapping_functions(latitude, height, elevation, day_of_year)
    el = max(elevation, 0.0)
    sin_el = sin(el)
    abs_lat = abs(rad2deg(latitude))
    # Seasonal term, southern hemisphere shifted by half a year
    cosy = cos(2π * ((day_of_year - 28.0) / 365.25 + (latitude < 0.0 ? 0.5 : 0.0)))
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
    saastamoinen_zenith_delays(latitude, height, relative_humidity) -> (hydrostatic, wet)

Saastamoinen zenith hydrostatic and wet delays in metres, from a standard atmosphere
(15 °C at sea level) at geodetic `latitude` (radians) and `height` (metres, clamped at 0),
with `relative_humidity` (0…1) for the wet term: RTKLIB's `tropmodel` at zenith.
"""
function saastamoinen_zenith_delays(latitude, height, relative_humidity)
    hgt = max(height, 0.0)
    # Standard atmosphere
    pressure = 1013.25 * (1.0 - 2.2557e-5 * hgt)^5.2568          # hPa
    temperature = 15.0 - 6.5e-3 * hgt + 273.16                  # K (15 °C at sea level)
    e = 6.108 * relative_humidity *
        exp((17.15 * temperature - 4684.0) / (temperature - 38.45))  # hPa
    hydrostatic =
        0.0022768 * pressure /
        (1.0 - 0.00266 * cos(2.0 * latitude) - 0.00028 * hgt / 1.0e3)
    wet = 0.002277 * (1255.0 / temperature + 0.05) * e
    return hydrostatic, wet
end

"""
    tropospheric_delay(elevation, lla, day_of_year; relative_humidity = 0.7) -> Float64

Slant tropospheric delay in metres at `elevation` (radians) from `lla` (a
`Geodesy.LLA`), the geometry [`ionospheric_delay`](@ref) takes too: the
[`saastamoinen_zenith_delays`](@ref) mapped by the [`niell_mapping_functions`](@ref).
Finite and monotone down to the horizon (≈ 90 m at sea level), and the horizon value
below it; `0.0` outside the valid heights, −100 m … 10 km.
"""
function tropospheric_delay(elevation, lla, day_of_year;
    relative_humidity = _DEFAULT_RELATIVE_HUMIDITY)
    (lla.alt < -100.0 || lla.alt > 1.0e4) && return 0.0
    latitude = deg2rad(lla.lat)
    height = max(lla.alt, 0.0)
    zenith_hydrostatic, zenith_wet = saastamoinen_zenith_delays(latitude, height,
        relative_humidity)
    hydrostatic_mapping, wet_mapping =
        niell_mapping_functions(latitude, height, elevation, day_of_year)
    return zenith_hydrostatic * hydrostatic_mapping + zenith_wet * wet_mapping
end

"""
    day_of_year(epoch::TAITime) -> Int

Day of year (1–366) of `epoch`, the seasonal argument of [`tropospheric_delay`](@ref);
date a time of week as [`week_start_epoch`](@ref)`(measurement) + time_of_week`. Read off
the TAI calendar label, whose 37 s lead on UTC is negligible here.
"""
day_of_year(epoch::TAITime) = dayofyear(DateTime(epoch))

