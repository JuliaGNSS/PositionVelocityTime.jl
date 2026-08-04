
# PositionVelocityTime.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/dev/)

Calculates position and time by using GNSS data

## Features
* Estimation of user position and time
* Calculates satellite position
* Precision estimation (GDOP)

## Preparing

### Install
```julia
julia> ]
pkg> add PositionVelocityTime
```

Decoded data and code phase of satellite must be combined in the provided `SatelliteState` struct. 
```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder
# decode satellite
gpsl1 = GPSL1CA()
sat_state = SatelliteState(
    decoder = decoder,
    system = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase # optional, in radians
)
```
The declaration of `carrier_phase` is optional due to its small effect on the user position.
`code_phase` is in chips and `carrier_phase` in radians, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

Alternatively, a `Tracking.TrackedSat` can be passed to `SatelliteState` instead of
`code_phase`, `carrier_doppler` and `carrier_phase` — `tracked_sat` below is what
`Tracking.get_sat_state` returns for a tracked satellite:
```julia
using Tracking
# track and decode satellite
sat_state = SatelliteState(decoder, gpsl1, tracked_sat)
```

## Usage

### User position calculation
The function 
```julia
calc_pvt(sat_states)
``` 
provides a complete position calculation. A fix needs at least 4 healthy, fully decoded
satellites (more for a multi-GNSS or multi-band set); when the epoch cannot be solved,
the previous solution is returned unchanged instead of an error, so a receiver can pass
whatever it currently tracks.