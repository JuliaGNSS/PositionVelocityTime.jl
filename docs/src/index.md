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
signal observation can be passed to `calc_pvt` safely — it is recognised and excluded
from the solve — but it contributes no pseudorange, and the HAS corrections it
carries are not yet applied here.

## Installation

```julia
using Pkg
Pkg.add("PositionVelocityTime")
```

## Quick Start

A satellite's decoder and code phase are combined in a [`SignalObservation`](@ref):

```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder

gpsl1 = GPSL1CA()
observation = SignalObservation(
    decoder = decoder,
    signal = gpsl1,
    code_phase = code_phase,
    carrier_doppler = carrier_doppler,
    carrier_phase = carrier_phase,  # optional, in cycles
)
```

`code_phase` is in chips and `carrier_phase` in cycles, matching `Tracking`'s
`get_code_phase` and `get_carrier_phase`.

With Tracking.jl loaded, a `Tracking.TrackedSat` (what
`Tracking.get_sat_state` returns) supplies the code phase, carrier Doppler and phase:

```julia
using Tracking
observation = SignalObservation(;
    decoder,
    signal = gpsl1,
    code_phase = Tracking.get_code_phase(tracked_sat),
    carrier_doppler = Tracking.get_carrier_doppler(tracked_sat),
    carrier_phase = Tracking.get_carrier_phase(tracked_sat),
)
```

Collect observations with the same concrete decoder and signal types in a vector or
`Dictionaries.Dictionary`. Pass one collection directly, or combine collections in a
tuple or `NamedTuple`:

```julia
previous = PVTSolution()
pvt = calc_pvt((
    gps = gps_observations,
    galileo = galileo_observations,
), previous)
lla = get_lla(pvt)  # latitude, longitude, altitude
```

Each collection must have one concrete `SignalObservation` element type, including a
`Dictionary`'s value type. This keeps preparation statically dispatched across a mixed
epoch; a pooled mixed vector with an abstract element type is refused at entry. A single
collection needs no tuple around it: `calc_pvt(gps_observations, previous)`. Collection
names are ignored. For an empty collection, preserve its concrete element type with
`gps_observations[1:0]` or supply that concrete type explicitly when constructing it.

Satellites from different constellations may be combined this way. One receiver clock
bias is estimated per time system, reported as effective `pvt.inter_system_biases` relative to
`pvt.reference_system`, so a fix over `M` time systems needs `3 + M` satellites, unless a
broadcast offset between two systems (Galileo's GGTO, BeiDou's BGTO) ties their clocks,
which brings a mixed fix down to 4. An epoch that cannot be solved returns the `prev_pvt`
[`calc_pvt`](@ref) was given rather than throwing, so a receiver can hand over whatever it
tracks.

These effective clock differences assume one receiver delay shared by all systems on
each band. Residual receiver delays can enter them, and disconnected constellation ×
band coverage also folds differences between reference-band delays into them. A
broadcast clock constraint fixes a difference under the same receiver-delay assumption;
it does not measure the receiver's differential delays.

### Solving without allocating

A receiver solving a stream of epochs can hand [`calc_pvt!`](@ref) the solution to reuse
and a [`PVTWorkspace`](@ref); once both have held an epoch of that size, a solve allocates
nothing:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for epoch in epochs
    pvt = calc_pvt!(pvt, workspace, epoch, pvt)   # reuse `pvt`, seeded from itself
end
```

Keep the returned solution; [`calc_pvt!`](@ref) says what it reuses and overwrites.
