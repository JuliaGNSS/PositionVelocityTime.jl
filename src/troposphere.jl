# ===========================================================================
#  Tropospheric correction — Saastamoinen model
#
#  The tropospheric delay is non-dispersive — the same on every GNSS frequency.
#  Unlike the ionosphere there are no broadcast coefficients: the delay is a
#  "blind" function of the user height and a standard atmosphere, mapped to the
#  line of sight by the satellite elevation. This is the Saastamoinen model as
#  implemented in RTKLIB's `tropmodel` (and used by GNSS-SDR / PocketSDR for
#  single-point positioning), with the relative humidity fixed at 70 % and a
#  simple 1/cos(z) obliquity mapping.
# ===========================================================================

const _DEFAULT_RELATIVE_HUMIDITY = 0.7  # GNSS-SDR / RTKLIB default

# Lowest elevation the obliquity mapping is evaluated at. The 1/cos(z) mapping
# diverges towards the horizon, so a satellite a few tenths of a degree up would be
# credited a delay of hundreds of metres — enough to drag a least-squares fix off.
# Below this bound the delay computed *at* the bound is reused, which keeps the
# correction finite, monotone and continuous across the horizon instead of jumping
# to zero there.
#
# This is a numerical guard, not an estimate: at the bound it yields ≈ 48 m where a
# curved-atmosphere mapping gives ≈ 37 m, while a genuinely grazing ray is ≈ 90 m.
# It is defensible only because 1/sin(E) is unusable this low anyway (see #62) and
# such observations should be masked or de-weighted rather than trusted.
#
# The value is taken from Orekit's
# `ModifiedSaastamoinenModel.DEFAULT_LOW_ELEVATION_THRESHOLD`, but note their bound
# guards a different failure: Orekit carries Saastamoinen's −B·tan²z curvature term,
# which goes quadratically negative and inverts the sign of the total delay below
# E ≈ 1.9°, so 0.05 rad sits just above *that* pathology. Ours has no such term and
# fails by diverging instead, so the value is borrowed rather than derived — it
# happens to coincide with the 3° lower validity limit that the mapping-function
# literature (Niell 1996, VMF3, IERS Conventions ch. 9) settles on.
const _LOW_ELEVATION_THRESHOLD = 0.05  # rad ≈ 2.87°, caps the mapping at 1/sin ≈ 20

"""
    tropospheric_delay(elevation, lla; humidity = 0.7) -> Float64

Slant tropospheric group delay in metres for a satellite at `elevation` (radians)
seen from the user geodetic position `lla` (a `Geodesy.LLA`), using the
Saastamoinen model driven by a standard atmosphere. The delay is non-dispersive
(frequency independent). `humidity` is the relative humidity (0…1) used for the
wet component. Elevations below 0.05 rad (≈ 2.87°), including negative ones, are
treated as 0.05 rad, so the delay saturates near the horizon instead of diverging.
Returns `0.0` when the user height is outside the model's valid range
(−100 m … 10 km).

The geometry is taken precomputed so a whole-epoch correction shares the user
geodetic conversion and the elevation with [`ionospheric_delay`](@ref).
"""
function tropospheric_delay(elevation, lla; humidity = _DEFAULT_RELATIVE_HUMIDITY)
    return saastamoinen_delay(deg2rad(lla.lat), lla.alt, elevation, humidity)
end

"""
    saastamoinen_delay(latitude, height, elevation, humidity) -> Float64

Saastamoinen slant tropospheric delay in metres. `latitude`/`elevation` are in
radians and `height` is the geodetic (ellipsoidal) height in metres. The standard
atmosphere (pressure, temperature, water-vapour partial pressure) is derived from
the height with a 15 °C sea-level temperature, and the hydrostatic and wet zenith
delays are mapped to the slant direction by `1/cos(z)` with `z = π/2 − elevation`.
Mirrors RTKLIB's `tropmodel`, except that the elevation is bounded from below by
`_LOW_ELEVATION_THRESHOLD` (as in Orekit's `ModifiedSaastamoinenModel`): the mapping
would otherwise grow without bound towards the horizon, so elevations below the
threshold — including satellites below the horizon — reuse the delay at the
threshold, roughly 20 × the zenith delay.

The flat-slab `1/cos(z)` mapping over-predicts at low elevation, by ≈ 14 % (≈ 3 m) at
5° and ≈ 4 % at 10°, since it ignores Earth curvature and ray bending; the bound caps
the divergence but not that bias. See #62.
"""
function saastamoinen_delay(latitude, height, elevation, humidity)
    (height < -100.0 || height > 1.0e4) && return 0.0
    hgt = height < 0.0 ? 0.0 : height
    # Standard atmosphere
    pressure = 1013.25 * (1.0 - 2.2557e-5 * hgt)^5.2568          # hPa
    temperature = 15.0 - 6.5e-3 * hgt + 273.16                  # K (15 °C at sea level)
    e = 6.108 * humidity * exp((17.15 * temperature - 4684.0) / (temperature - 38.45))  # hPa
    # Saastamoinen hydrostatic and wet slant delays (1/cos(z) obliquity mapping),
    # with the elevation bounded from below so the mapping cannot diverge
    z = π / 2 - max(elevation, _LOW_ELEVATION_THRESHOLD)
    trph =
        0.0022768 * pressure /
        (1.0 - 0.00266 * cos(2.0 * latitude) - 0.00028 * hgt / 1.0e3) / cos(z)
    trpw = 0.002277 * (1255.0 / temperature + 0.05) * e / cos(z)
    return trph + trpw
end
