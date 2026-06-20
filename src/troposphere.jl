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

"""
    tropospheric_delay(elevation, lla; humidity = 0.7) -> Float64

Slant tropospheric group delay in metres for a satellite at `elevation` (radians)
seen from the user geodetic position `lla` (a `Geodesy.LLA`), using the
Saastamoinen model driven by a standard atmosphere. The delay is non-dispersive
(frequency independent). `humidity` is the relative humidity (0…1) used for the
wet component. Returns `0.0` when the satellite is at or below the horizon or the
user height is outside the model's valid range (−100 m … 10 km).

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
Mirrors RTKLIB's `tropmodel`.
"""
function saastamoinen_delay(latitude, height, elevation, humidity)
    (height < -100.0 || height > 1.0e4 || elevation <= 0.0) && return 0.0
    hgt = height < 0.0 ? 0.0 : height
    # Standard atmosphere
    pressure = 1013.25 * (1.0 - 2.2557e-5 * hgt)^5.2568          # hPa
    temperature = 15.0 - 6.5e-3 * hgt + 273.16                  # K (15 °C at sea level)
    e = 6.108 * humidity * exp((17.15 * temperature - 4684.0) / (temperature - 38.45))  # hPa
    # Saastamoinen hydrostatic and wet slant delays (1/cos(z) obliquity mapping)
    z = π / 2 - elevation
    trph =
        0.0022768 * pressure /
        (1.0 - 0.00266 * cos(2.0 * latitude) - 0.00028 * hgt / 1.0e3) / cos(z)
    trpw = 0.002277 * (1255.0 / temperature + 0.05) * e / cos(z)
    return trph + trpw
end
