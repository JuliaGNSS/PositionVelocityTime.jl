# PositionVelocityTime.jl

Calculates position, velocity and time from GNSS satellite measurements.

## Features

- Estimation of user position, velocity and time, with the receiver clock bias and drift
- Satellite position and velocity from the broadcast ephemerides: the directly broadcast
  Keplerian ones (GPS LNAV, Galileo I/NAV and F/NAV, BeiDou D1/D2), the quasi-Keplerian
  ones (GPS CNAV/CNAV-2, BeiDou B-CNAV1/2/3), and the BeiDou GEO reference frame
- Dilution of precision (GDOP, PDOP, HDOP, VDOP, TDOP)
- Ionospheric (NTCM-G, BDGIM, Klobuchar) and tropospheric (Saastamoinen with the Niell
  mapping functions) corrections from the navigation data alone
- Support for GPS (L1 C/A, L2C, L5, L1C), Galileo (E1B, E5a, E5b) and BeiDou
  (B1I, B3I, B1C, B2a, B2b), including combined multi-GNSS solutions. Each
  measurement is one satellite-band pseudorange; the group-delay/ISC correction is
  selected by the signal the range was generated on (which may be a pilot such as
  GPS L1C-P, Galileo E1C or BeiDou B2a-Q), while the ephemeris and clock come from
  the band's data-component decoder.
- One receiver clock bias per GNSS time system (GPST, GST, BDT). Where the satellites
  cannot separate two of them, a broadcast offset between the systems — Galileo's GGTO,
  BeiDou's BGTO, or GPS's CNAV/CNAV-2 GGTO — ties them together as an extra
  observation

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
satellite_state = SatelliteState(
    decoder = decoder,
    signal = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase,  # optional, in radians
)
```

`code_phase` is in chips and `carrier_phase` in radians, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

With Tracking.jl 9.1 or later loaded, you can pass a `Tracking.TrackedSat` directly — `tracked_satellite` is what
`Tracking.get_sat_state` returns for a tracked satellite, and the code phase, carrier
Doppler, and carrier phase are read off it:

```julia
using Tracking
satellite_state = SatelliteState(decoder, gpsl1, tracked_satellite)
```

Group the satellites by their ranging signal and compute the PVT solution:

```julia
using PositionVelocityTime: SignalGroup

pvt = calc_pvt((
    gps = SignalGroup(GPSL1CA(), gps_satellite_states),
    galileo = SignalGroup(GalileoE1B(), galileo_satellite_states),
))
lla = get_lla(pvt)  # latitude, longitude, altitude
```

Each group holds the satellites tracked on one signal, as a `Dictionary` keyed by PRN or
as a plain vector. A single group needs no NamedTuple around it —
`calc_pvt(SignalGroup(GPSL1CA(), gps_satellite_states))` is a complete one-constellation solve.
With Tracking.jl 9.1 or later loaded, a whole `TrackState` and its decoders become groups in one call:

```julia
using Tracking
pvt = calc_pvt(PositionVelocityTime.signal_groups(track_state, decoders))
```

Satellites from different constellations may be combined this way. Because each GNSS
references its broadcasts to its own system time, [`calc_pvt`](@ref) estimates one
receiver clock bias per GNSS time system, so a combined fix needs at least `3 + M`
satellites for `M` distinct systems — unless a broadcast offset between two of them
(GGTO, BGTO) ties their clocks together, which brings a mixed fix down to 4 satellites.
The per-system clock offsets are reported as `pvt.inter_system_biases` relative to
`pvt.reference_system`.

If too few healthy satellites are tracked to solve the constellation — or the geometry
turns out to be degenerate — [`calc_pvt`](@ref) returns the `prev_pvt` it was given (the
origin solution by default) rather than throwing, so a receiver can hand it whatever it
currently tracks each epoch and carry the last solution forward.

### Solving without allocating

`calc_pvt` builds a new solution, containers and all, every epoch. A receiver solving a
steady stream of epochs can instead hand [`calc_pvt!`](@ref) the previous solution as
the one to reuse, and a [`PVTWorkspace`](@ref) for its scratch; once both have held an
epoch of that size, a solve allocates nothing:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for groups in epochs
    pvt = calc_pvt!(pvt, workspace, groups, pvt)   # reuse `pvt`, seeded from itself
end
```

A [`PVTSolution`](@ref) is immutable, so `calc_pvt!` returns the new fix; what it reuses
are the `satellites` and bias containers of its first argument, which the returned solution
takes over. That argument is the only thing overwritten — `prev_pvt` is only read, so
passing a different solution there keeps it intact — and it is spent by the call: keep
the returned solution, and `copy(pvt)` any solution you want to keep beyond the next
call. An epoch that cannot be solved returns `prev_pvt`, the answer `calc_pvt` would
have returned.
