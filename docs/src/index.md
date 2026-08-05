# PositionVelocityTime.jl

Calculates position, velocity, and time from GNSS satellite measurements.

## Features

- Estimation of user position, velocity, and time
- Satellite position and velocity calculation from orbital parameters (LNAV and
  CNAV/CNAV-2 quasi-Keplerian ephemerides)
- Dilution of Precision (DOP) computation
- Weighted least squares from the reported C/N₀ and the satellite elevation, so a
  marginal satellite contributes its geometry without dragging the fix, plus a formal
  accuracy in metres alongside the (purely geometric) DOP
- Support for GPS (L1 C/A, L2C, L5, L1C) and Galileo (E1B, E5a), including combined
  multi-GNSS solutions. Each measurement is one satellite-band pseudorange; the
  group-delay/ISC correction is selected by the signal the range was generated on
  (which may be a pilot such as GPS L1C-P or Galileo E1C), while the ephemeris and
  clock come from the band's data-component decoder.

## Installation

```julia
using Pkg
Pkg.add("PositionVelocityTime")
```

## Quick Start

Decoded data and code phase of a satellite must be combined in the [`SatelliteState`](@ref) struct:

```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder
using Unitful: dBHz

gpsl1 = GPSL1CA()
sat_state = SatelliteState(
    decoder = decoder,
    system = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase,  # optional, in radians
    cn0 = 42.0dBHz,                 # optional, enables weighted least squares
)
```

`code_phase` is in chips and `carrier_phase` in radians, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

`cn0` is optional too, and matches `Tracking.estimate_cn0`. Supplying it lets
[`calc_pvt`](@ref) weight each satellite by the measurement uncertainty its C/N₀ and
elevation imply (see [`PositionVelocityTime.pseudorange_variance`](@ref)) rather than
treating a marginal satellite as being as precise as a strong one. Supply
`pseudorange_variance` (e.g. `(3.0m)^2`) instead to use your own error model; with
neither, the solve is the plain unweighted one and behaves exactly as before.

Alternatively, pass a `Tracking.TrackedSat` directly — `tracked_sat` is what
`Tracking.get_sat_state` returns for a tracked satellite, and the code phase, carrier
Doppler, carrier phase, and C/N₀ are read off it:

```julia
using Tracking
sat_state = SatelliteState(decoder, gpsl1, tracked_sat)
```

Compute the PVT solution:

```julia
pvt = calc_pvt(sat_states)
lla = get_LLA(pvt)  # latitude, longitude, altitude
```

Satellites from different constellations may be passed together. Because each GNSS
references its broadcasts to its own system time, [`calc_pvt`](@ref) estimates one
receiver clock bias per GNSS time system, so a combined fix needs at least `3 + M`
satellites for `M` distinct systems. The per-system clock offsets are reported as
`pvt.inter_system_biases` relative to `pvt.reference_system`.

If too few healthy satellites are tracked to solve the constellation — or the geometry
turns out to be degenerate — [`calc_pvt`](@ref) returns the `prev_pvt` it was given (the
origin solution by default) rather than throwing, so a receiver can hand it whatever it
currently tracks each epoch and carry the last solution forward.

## Solution quality

Two quantities describe how good a fix is, and they are deliberately separate:

- `pvt.dop` ([`PositionVelocityTime.DOP`](@ref)) is the **geometry** alone,
  `(HᵀH)⁻¹` — unitless, and unchanged by measurement weighting, so it means what every
  other receiver means by GDOP/PDOP/HDOP/VDOP/TDOP.
- `pvt.accuracy` ([`PositionVelocityTime.FormalAccuracy`](@ref)) is the **formal 1σ
  accuracy in metres**, `(HᵀWH)⁻¹`, i.e. the geometry with the a-priori measurement
  uncertainties folded in. With no uncertainty reported it reduces to DOP × a nominal
  UERE.

Per satellite, `pvt.sats` reports the post-fit residual (raw metres, weighted or not)
and the σ that satellite was weighted by, whose ratio is the normalised residual used
for fault detection.
