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

Group the satellites by their ranging signal and compute the PVT solution:

```julia
using PositionVelocityTime: SignalGroup

pvt = calc_pvt((
    gps = SignalGroup(GPSL1CA(), gps_sat_states),
    galileo = SignalGroup(GalileoE1B(), galileo_sat_states),
))
lla = get_LLA(pvt)  # latitude, longitude, altitude
```

Each group holds the satellites tracked on one signal, as a `Dictionary` keyed by PRN or
as a plain vector. A single group needs no NamedTuple around it —
`calc_pvt(SignalGroup(GPSL1CA(), gps_sat_states))` is a complete one-constellation solve.
With `Tracking` loaded, a whole `TrackState` and its decoders become groups in one call:

```julia
using Tracking
pvt = calc_pvt(PositionVelocityTime.signal_groups(track_state, decoders))
```

Satellites from different constellations may be combined this way. Because each GNSS
references its broadcasts to its own system time, [`calc_pvt`](@ref) estimates one
receiver clock bias per GNSS time system, so a combined fix needs at least `3 + M`
satellites for `M` distinct systems. The per-system clock offsets are reported as
`pvt.inter_system_biases` relative to `pvt.reference_system`.

!!! note "Migrating from 5.x"

    `calc_pvt` used to take a flat `AbstractVector{<:SatelliteState}`. Pooling several
    constellations into one vector made its element type abstract, which turned every
    per-satellite call inside the solve into a dynamic dispatch — a mixed 14-satellite
    epoch cost about five times a single-constellation one. Groups fix that at the
    source, and behind them the solver now works on one flat, parameter-free
    measurement row and compiles once for every constellation mix.

    There is deliberately no conversion function from the old input: `calc_pvt` refuses
    a pooled vector with an error saying what to build. A converter would have had to
    infer the grouping from whatever signals the vector happened to hold, which is a
    runtime property — so the conversion itself would be inference-blind, reintroducing
    in one place the dispatch the grouping removes everywhere else. Build the groups
    where the satellites are tracked, which is the shape `Tracking.jl` and
    `GNSSReceiver.jl` already carry, and hand the same groups to `calc_pvt`.

    The measurement-model surface moved with it: its per-satellite functions now take
    the [`PositionVelocityTime.SatelliteMeasurement`](@ref) rows that
    [`PositionVelocityTime.collect_measurements`](@ref) produces, rather than
    `SatelliteState`s plus parallel classification vectors. See
    [The Measurement-Model Surface](@ref).

    `pvt.time` is a [`TAITime`](@ref) rather than an AstroTime `TAIEpoch`, and AstroTime
    is no longer loaded with this package: `using AstroTime` and `TAIEpoch(pvt.time)`
    converts it exactly. `position` and `velocity` are `ECEF{Float64}`, and the keys of
    `inter_system_biases` (and `reference_system`) are
    [`PositionVelocityTime.SupportedTimeSystem`](@ref)s. These are what let the solver
    compile into a trimmed executable — see [Trimmed Executables](@ref).

If too few healthy satellites are tracked to solve the constellation — or the geometry
turns out to be degenerate — [`calc_pvt`](@ref) returns the `prev_pvt` it was given (the
origin solution by default) rather than throwing, so a receiver can hand it whatever it
currently tracks each epoch and carry the last solution forward.

### Solving without allocating

`calc_pvt` returns a new solution every epoch. A receiver solving a steady stream of
epochs can instead hand [`calc_pvt!`](@ref) the solution to overwrite and a reusable
[`PVTWorkspace`](@ref) for its scratch; once both have held an epoch of that size, a
solve allocates nothing:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for groups in epochs
    calc_pvt!(pvt, workspace, groups, pvt)   # overwrite `pvt`, seeded from itself
end
```

The overwrite is explicit — only the first argument is written, and `prev_pvt` is only
read, so passing a different solution there keeps it intact. An epoch that cannot be
solved leaves `pvt` holding a copy of `prev_pvt`, the answer `calc_pvt` would have
returned.
