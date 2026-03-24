# PositionVelocityTime.jl

Calculates position, velocity, and time from GNSS satellite measurements.

## Features

- Estimation of user position, velocity, and time
- Satellite position and velocity calculation from orbital parameters
- Dilution of Precision (DOP) computation
- Support for GPS L1 and Galileo E1B

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
