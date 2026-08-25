# PositionVelocityTime.jl

Calculates position, velocity, and time from GNSS satellite measurements.

## Features

- Estimation of user position, velocity, and time
- Satellite position and velocity calculation from orbital parameters: the
  directly-broadcast Keplerian ephemerides (GPS LNAV, Galileo I/NAV and F/NAV,
  BeiDou D1/D2), the quasi-Keplerian ones (GPS CNAV/CNAV-2, BeiDou B-CNAV1/2/3),
  and the BeiDou GEO reference frame
- Dilution of Precision (DOP) computation
- Support for GPS (L1 C/A, L2C, L5, L1C), Galileo (E1B, E5a, E5b) and BeiDou
  (B1I, B3I, B1C, B2a, B2b), including combined multi-GNSS solutions. Each
  measurement is one satellite-band pseudorange; the group-delay/ISC correction is
  selected by the signal the range was generated on (which may be a pilot such as
  GPS L1C-P, Galileo E1C or BeiDou B2a-Q), while the ephemeris and clock come from
  the band's data-component decoder.
- One receiver clock bias per GNSS time system (GPST, GST, BDT), collapsed onto GPS
  Time using the broadcast offset — Galileo's GGTO or BeiDou's BGTO — when the
  geometry is too weak to estimate it

Galileo E6-B is decoded by `GNSSDecoder` but carries no ephemeris of its own: its
C/NAV message is the Galileo High Accuracy Service correction channel. An E6-B
satellite state can be passed to `calc_pvt` safely — it is recognised and excluded
from the solve — but it contributes no pseudorange, and the HAS corrections it
carries are not yet applied here.

## Installation

```julia
using Pkg
Pkg.add("PositionVelocityTime")
```

## Quick Start

Decoded data and code phase of a satellite must be combined in the [`SatelliteState`](@ref) struct:

```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder

gpsl1 = GPSL1CA()
sat_state = SatelliteState(
    decoder = decoder,
    system = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase,  # optional, in radians
)
```

`code_phase` is in chips and `carrier_phase` in radians, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

Alternatively, pass a `Tracking.TrackedSat` directly — `tracked_sat` is what
`Tracking.get_sat_state` returns for a tracked satellite, and the code phase, carrier
Doppler, and carrier phase are read off it:

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
