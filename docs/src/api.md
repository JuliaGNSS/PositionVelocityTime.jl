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
get_num_used_sats
get_sat_info
get_frequency_offset
```

## Multi-GNSS Internals

These helpers classify satellites when [`calc_pvt`](@ref) combines constellations and
bands; documented for reference. `get_time_system` drives receiver-clock grouping (one
bias per time system) and `get_frequency_band` drives inter-frequency-bias grouping (one
bias per band), while `get_signal_id` is the per-signal identity used in the `sats` key
of [`PVTSolution`](@ref).

```@docs
PositionVelocityTime.get_time_system
PositionVelocityTime.get_frequency_band
PositionVelocityTime.get_signal_id
```
