# PositionVelocityTime.jl

Calculates position, velocity, and time from GNSS satellite measurements.

## Features

- Estimation of user position, velocity, and time
- Satellite position and velocity calculation from orbital parameters (LNAV and
  CNAV/CNAV-2 quasi-Keplerian ephemerides)
- Dilution of Precision (DOP) computation
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

gpsl1 = GPSL1()
sat_state = SatelliteState(
    decoder = decoder,
    system = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase,  # optional
)
```

Alternatively, pass a `Tracking.SatState` directly:

```julia
using Tracking
sat_state = SatelliteState(decoder, gpsl1, sat_state)
```

With at least 4 satellite states, compute the PVT solution:

```julia
pvt = calc_pvt(sat_states)
lla = get_LLA(pvt)  # latitude, longitude, altitude
```

Satellites from different constellations may be passed together. Because each GNSS
references its broadcasts to its own system time, [`calc_pvt`](@ref) estimates one
receiver clock bias per GNSS time system, so a combined fix needs at least `3 + M`
satellites for `M` distinct systems. The per-system clock offsets are reported as
`pvt.inter_system_biases` relative to `pvt.reference_system`.
