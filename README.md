
# PositionVelocityTime.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/dev/)

Calculates position and time by using GNSS data

## Features
* Estimation of user position, velocity and time
* Calculates satellite position and velocity
* Precision estimation (GDOP)
* GPS (L1 C/A, L2C, L5, L1C), Galileo (E1B, E5a, E5b) and BeiDou (B1I, B3I, B1C,
  B2a, B2b), combined in one multi-GNSS solution

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
Satellite states are handed over grouped by their ranging signal, one
`PositionVelocityTime.SignalGroup` per signal:
```julia
using PositionVelocityTime: SignalGroup

calc_pvt((
    gps = SignalGroup(GPSL1CA(), gps_sat_states),
    galileo = SignalGroup(GalileoE1B(), galileo_sat_states),
))
```
provides a complete position calculation. A single group can be passed on its own, and
with `Tracking` loaded `PositionVelocityTime.signal_groups(track_state, decoders)` builds
a whole epoch's groups from a `TrackState`. Grouping is what keeps the solve type-stable
across constellations — see the migration note in the documentation if you are coming
from 5.x, where `calc_pvt` took a flat vector of satellite states directly.

A fix needs at least 4 healthy, fully decoded
satellites (more for a multi-GNSS or multi-band set); when the epoch cannot be solved,
the previous solution is returned unchanged instead of an error, so a receiver can pass
whatever it currently tracks.

To solve epoch after epoch without allocating, overwrite one solution in place with
`calc_pvt!`, reusing a `PVTWorkspace` for the scratch buffers:
```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
calc_pvt!(pvt, workspace, groups, pvt)  # writes only `pvt`; seeded from itself
```
The estimated time `pvt.time` is a `TAITime` — whole TAI seconds since J2000 plus a
fraction. With AstroTime loaded, `TAIEpoch(pvt.time)` converts it exactly.

### Trimmed executables
The solver compiles into a standalone executable with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)'s `juliac --trim=safe` (Julia 1.12+);
`test/trim` holds such an app and a check that builds it and compares its output with a
regular Julia session.
