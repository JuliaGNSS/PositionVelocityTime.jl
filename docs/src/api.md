# API Reference

## Types

```@docs
SatelliteState
PVTSolution
SatInfo
InterFrequencyBias
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
PositionVelocityTime.saastamoinen_zenith_delays
PositionVelocityTime.niell_mapping_functions
PositionVelocityTime.KlobucharParams
PositionVelocityTime.klobuchar_group_delay
PositionVelocityTime.BeiDouKlobucharParams
PositionVelocityTime.beidou_klobuchar_group_delay
PositionVelocityTime.NTCMGParams
PositionVelocityTime.BDGIMParams
PositionVelocityTime.klobuchar_params
PositionVelocityTime.ntcm_g_params
PositionVelocityTime.bdgim_params
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

## The Measurement-Model Surface

The pieces [`calc_pvt`](@ref) is assembled from, documented so a consumer that
runs its own estimator over the same measurement model — a navigation filter
closing tracking loops, for instance — reuses this package's model as a
stated contract. The names are deliberately not exported (they are solver
internals, not part of every user's vocabulary); bind them explicitly with
`using PositionVelocityTime: calc_corrected_time, …`. Times returned or taken
here are seconds-of-week counts on each satellite's own system scale unless a
function says otherwise.

```@docs
PositionVelocityTime.SPEED_OF_LIGHT
PositionVelocityTime.calc_corrected_time
PositionVelocityTime.calc_satellite_clock_drift
PositionVelocityTime.get_sat_position
PositionVelocityTime.get_sat_velocity
PositionVelocityTime.fold_week_crossover
PositionVelocityTime.BiasColumns
PositionVelocityTime.decide_bias_layout
PositionVelocityTime.calc_ρ_hat!
PositionVelocityTime.calc_H
PositionVelocityTime.calc_line_of_sight
PositionVelocityTime.calc_DOP
PositionVelocityTime.time_scale_offset_to_gpst
PositionVelocityTime.calc_time_scale_offsets
PositionVelocityTime.gpst_offset_available
PositionVelocityTime.calc_gpst_offset
PositionVelocityTime.get_week
PositionVelocityTime.system_start_epoch
PositionVelocityTime.day_of_year
PositionVelocityTime.predict_atmospheric_delays
PositionVelocityTime.calc_course_over_ground
```

The unexported internals the surface above links to, documented for reference:

```@docs
PositionVelocityTime.calc_H!
PositionVelocityTime.user_position
PositionVelocityTime.num_lsq_params
PositionVelocityTime.calc_gpst_range_offsets
PositionVelocityTime.positive_definite_cholesky
```

## Multi-GNSS Classification

When [`calc_pvt`](@ref) combines constellations and bands it classifies each satellite by
three keys, all provided by GNSSSignals (3.3+) and read from a satellite's ranging signal:
`GNSSSignals.get_time_system` (a `GNSSSignals.TimeSystem`, i.e. `GPST()`/`GST()`/`BDT()`) drives
receiver-clock grouping (one bias per time system); `GNSSSignals.get_band_id` (e.g. `:L1`,
`:L5`) drives inter-frequency-bias grouping (one bias per band); and
`GNSSSignals.get_signal_id` (e.g. `:GPSL1CA`) is the per-signal identity used in the `sats`
key of [`PVTSolution`](@ref). This package forwards each of the latter two to a
[`SatelliteState`](@ref).

The receiver inter-frequency biases and their reference bands (reported per
[`InterFrequencyBias`](@ref)) are laid out from the constellation × band coverage graph:

```@docs
PositionVelocityTime.band_ifb_layout
```
