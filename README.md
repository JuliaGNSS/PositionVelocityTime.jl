# PositionVelocityTime.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/PositionVelocityTime.jl/dev/)

Calculates position, velocity and time from GNSS satellite measurements.

## Features
* Estimation of user position, velocity and time, with the receiver clock bias and drift
* Satellite position and velocity from the broadcast ephemerides
* Dilution of precision (GDOP, PDOP, HDOP, VDOP, TDOP)
* Ionospheric (NTCM-G, BDGIM, Klobuchar) and tropospheric (Saastamoinen with the Niell
  mapping functions) corrections from the navigation data alone
* GPS (L1 C/A, L2C, L5, L1C), Galileo (E1B, E5a, E5b) and BeiDou (B1I, B3I, B1C,
  B2a, B2b), combined in one multi-GNSS, multi-band solution

## Preparing

### Install
```julia
julia> ]
pkg> add PositionVelocityTime
```

The decoded data and the code phase of a satellite are combined in a `SignalObservation`:
```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder
# decode satellite
gpsl1 = GPSL1CA()
sat_state = SignalObservation(
    decoder = decoder,
    signal = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase # optional, in cycles
)
```
The `carrier_phase` is optional, as its effect on the user position is small.
`code_phase` is in chips and `carrier_phase` in cycles, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

With Tracking.jl loaded, a `Tracking.TrackedSat` can be passed to `SignalObservation`
instead of `code_phase`, `carrier_doppler` and `carrier_phase` — `tracked_sat` below
is what `Tracking.get_sat_state` returns for a tracked satellite:
```julia
using Tracking
# track and decode satellite
sat_state = SignalObservation(decoder, gpsl1, tracked_sat)
```

## Usage

### Position, velocity and time
Satellite states are handed over grouped by their ranging signal, one
`PositionVelocityTime.SignalGroup` per signal:
```julia
using PositionVelocityTime: SignalGroup

pvt = calc_pvt((
    gps = SignalGroup(GPSL1CA(), gps_sat_states),
    galileo = SignalGroup(GalileoE1B(), galileo_sat_states),
))
lla = get_lla(pvt)  # latitude, longitude, altitude
```
A single group can be passed on its own, and with Tracking.jl loaded
`PositionVelocityTime.signal_groups(track_state, decoders)` builds a whole epoch's groups
from a `TrackState`. Grouping is what keeps the solve type-stable across constellations.

Besides `position`, `velocity` and `time`, the solution reports the receiver `clock_bias`
(m) and `clock_drift` (m/s), the `dop`, the per-observation residuals in `observations`, and,
for a multi-GNSS or multi-band fix, the `inter_system_biases` relative to the
`reference_system` and the `inter_frequency_biases`. The inter-system biases are effective
receiver clock differences under a shared per-band receiver-delay model; unseparated
receiver delays, including reference-band differences in disconnected coverage, enter
those values. `broadcast_clock_constraints` records any broadcast offsets used to
constrain receiver-clock differences under that model.

A fix needs at least 4 healthy, fully decoded satellites (more for a multi-GNSS or
multi-band set, unless a broadcast time offset such as Galileo's GGTO ties two systems'
clocks together); when the epoch cannot be solved, the previous solution is returned
unchanged instead of an error, so a receiver can pass whatever it currently tracks.

To solve epoch after epoch without allocating, use `calc_pvt!`, which reuses the
containers of the solution passed as its first argument and a `PVTWorkspace` for the
scratch buffers. `PVTSolution` is immutable, so keep the returned solution:
```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
pvt = calc_pvt!(pvt, workspace, groups, pvt)  # reuses only `pvt`; seeded from itself
```
The estimated time `pvt.time` is GNSSSignals' `TAITime` — whole TAI seconds since J2000
plus a fraction. With AstroTime loaded, `TAIEpoch(pvt.time)` converts it exactly.

### Trimmed executables
The solver compiles into a standalone executable with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)'s `juliac --trim=safe` (Julia 1.12+);
`test/trim` holds such an app and a check that builds it and compares its output with a
regular Julia session.
