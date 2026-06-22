# API Reference

## Types

```@docs
SatelliteState
PVTSolution
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

```@docs
get_gdop
get_pdop
get_hdop
get_vdop
get_tdop
```

## Utilities

```@docs
get_num_used_sats
get_frequency_offset
```
