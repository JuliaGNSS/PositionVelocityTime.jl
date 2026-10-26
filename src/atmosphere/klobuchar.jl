# The two Klobuchar algorithms: GPS (IS-GPS-200) and BeiDou's variant
# (BDS-SIS-ICD-B1I).

"""
    klobuchar_group_delay(φ_u, λ_u, elevation, azimuth, time_of_week, α, β)

Klobuchar single-frequency ionospheric group delay for GPS L1 (IS-GPS-200N,
Fig. 20-4), returned in **seconds**. All angles are in **semicircles**:
`φ_u`/`λ_u` are the user geodetic latitude/longitude, `elevation`/`azimuth` (the ICD's
`E`/`A`) the satellite's. `time_of_week` is the GPS time of week in seconds (the ICD's
`t`). `α`/`β` are the 4-element Klobuchar coefficient tuples (SI units).
"""
function klobuchar_group_delay(φ_u, λ_u, elevation, azimuth, time_of_week, α, β)
    # Earth-centred angle between user and ionospheric pierce point (semicircles)
    ψ = 0.0137 / (elevation + 0.11) - 0.022
    # Geodetic latitude of the ionospheric pierce point (IPP), clamped per ICD
    φ_i = clamp(φ_u + ψ * cos(azimuth * π), -0.416, 0.416)
    # Geodetic longitude of the IPP
    λ_i = λ_u + ψ * sin(azimuth * π) / cos(φ_i * π)
    # Geomagnetic latitude of the IPP
    φ_m = φ_i + 0.064 * cos((λ_i - 1.617) * π)
    # Local time at the IPP (seconds), wrapped to [0, 86400)
    t = mod(4.32e4 * λ_i + time_of_week, 86400.0)
    # Obliquity / slant factor
    F = 1.0 + 16.0 * (0.53 - elevation)^3
    # Amplitude (s) and period (s) of the cosine model, with ICD floors
    AMP = max(α[1] + φ_m * (α[2] + φ_m * (α[3] + φ_m * α[4])), 0.0)
    PER = max(β[1] + φ_m * (β[2] + φ_m * (β[3] + φ_m * β[4])), 72000.0)
    x = 2π * (t - 50400.0) / PER
    return abs(x) < 1.57 ? F * (5.0e-9 + AMP * (1 - x^2 / 2 + x^4 / 24)) : F * 5.0e-9
end

"""
    beidou_klobuchar_group_delay(φ_u, λ_u, elevation, azimuth, time_of_week, α, β)

BeiDou's Klobuchar-variant ionospheric group delay along the B1I propagation
path (BDS-SIS-ICD-B1I-3.0 §5.2.4.7), returned in **seconds**. All angles are in
**radians**, the ICD's working unit: `φ_u`/`λ_u` are the user geodetic
latitude/longitude, `elevation`/`azimuth` (the ICD's `E`/`A`) the satellite's.
`time_of_week` is the BDT second of week (the ICD's `t_E`). `α`/`β` are the 4-element
coefficient tuples (SI units, per semicircle like the GPS set).

The half-cosine idea of `klobuchar_group_delay`, but the ICD replaces its fitted
approximations with exact spherical geometry and drops the geomagnetic frame:

- The Earth-central angle to the ionospheric pierce point is the exact
  `π/2 − E − arcsin(R/(R+h)·cos E)` at R = 6378 km, h = 375 km, and the
  obliquity the matching exact `1/√(1 − (R/(R+h)·cos E)²)`, instead of the GPS
  model's `0.0137/(E+0.11) − 0.022` and `1 + 16(0.53−E)³` fits.
- The latitude feeding the A₂/A₄ polynomials is the pierce point's *geographic*
  latitude from the spherical triangle — no `+0.064·cos(λ−1.617)` geomagnetic
  conversion, no ±0.416-semicircle clamp — and its longitude uses `arcsin`
  rather than a flat-Earth division.
- The period A₄ is clamped from above at 172800 s as well as from below at
  72000 s, and the cosine is evaluated exactly rather than by the two-term
  Taylor expansion inside |x| < 1.57.

The 5·10⁻⁹ s night floor and the |t − 50400| < A₄/4 day-time window are those
of the GPS model.
"""
function beidou_klobuchar_group_delay(φ_u, λ_u, elevation, azimuth, time_of_week, α, β)
    # Sine of the pierce point's zenith angle, R/(R+h)·cos E (R = 6378 km, h = 375 km).
    sinz = 6378.0 / (6378.0 + 375.0) * cos(elevation)
    # Earth-central angle between user and pierce point (radians), exact
    ψ = π / 2 - elevation - asin(sinz)
    # Geographic — not geomagnetic — latitude and longitude of the pierce point
    φ_M = asin(sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth))
    λ_M = λ_u + asin(sin(ψ) * sin(azimuth) / cos(φ_M))
    # Local time at the pierce point (seconds), wrapped to [0, 86400)
    t = mod(time_of_week + λ_M * 43200.0 / π, 86400.0)
    # Amplitude (s) and period (s) of the cosine model, with the ICD's clamps
    φ = φ_M / π
    A_2 = max(α[1] + φ * (α[2] + φ * (α[3] + φ * α[4])), 0.0)
    A_4 = clamp(β[1] + φ * (β[2] + φ * (β[3] + φ * β[4])), 72000.0, 172800.0)
    # Vertical delay I′z (s): a cosine within a quarter period of 14:00, else the floor
    x = 2π * (t - 50400.0) / A_4
    I_z = abs(t - 50400.0) < A_4 / 4 ? 5.0e-9 + A_2 * cos(x) : 5.0e-9
    # Slant delay along the B1I path, with the exact obliquity
    return I_z / sqrt(1.0 - sinz^2)
end
