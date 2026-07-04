# API Reference

## Types

```@docs
SatelliteState
PVTSolution
SatInfo
```

## PVT Computation

```@docs
calc_pvt
```

## Position and Velocity

```@docs
get_LLA
calc_satellite_position
calc_satellite_position_and_velocity
get_sat_enu
```

## Atmospheric Corrections

These corrections are applied automatically by [`calc_pvt`](@ref); they are
documented here for reference and for diagnostic use.

```@docs
PositionVelocityTime.select_ionospheric_correction
PositionVelocityTime.ionospheric_delay
PositionVelocityTime.tropospheric_delay
PositionVelocityTime.KlobucharParams
PositionVelocityTime.NTCMGParams
PositionVelocityTime._elevation_azimuth
```

## Dilution of Precision

The DOP values are read from the `dop` field of a [`PVTSolution`](@ref), e.g.
`pvt.dop.GDOP`.

```@docs
PositionVelocityTime.DOP
```

## Utilities

```@docs
get_sat_info
```

## Multi-GNSS Classification

When [`calc_pvt`](@ref) combines constellations and bands it classifies each satellite by
three keys, all provided by GNSSSignals (3.3+) and read from a satellite's ranging signal:
`GNSSSignals.get_time_system` (a `GNSSSignals.TimeSystem`, i.e. `GPST()`/`GST()`) drives
receiver-clock grouping (one bias per time system); `GNSSSignals.get_band_id` (e.g. `:L1`,
`:L5`) drives inter-frequency-bias grouping (one bias per band); and
`GNSSSignals.get_signal_id` (e.g. `:GPSL1CA`) is the per-signal identity used in the `sats`
key of [`PVTSolution`](@ref). This package forwards each of the latter two to a
[`SatelliteState`](@ref).
