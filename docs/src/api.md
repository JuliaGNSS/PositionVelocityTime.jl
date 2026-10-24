# API Reference

## Types

```@docs
SatelliteState
PVTSolution
TAITime
SatelliteInfo
InterFrequencyBias
```

With AstroTime loaded, `TAIEpoch(pvt.time)` and `TAITime(epoch)` convert between
[`TAITime`](@ref) and AstroTime's `TAIEpoch` exactly.

## Signal Groups

One epoch's measurements reach [`calc_pvt`](@ref) as *signal groups*: a `NamedTuple` of
[`PositionVelocityTime.SignalGroup`](@ref)s, one per ranging signal, mirroring
Tracking.jl's. Within a group every satellite shares one concrete type, so
[`PositionVelocityTime.collect_measurements`](@ref) dispatches statically; behind it the
solver sees one flat [`PositionVelocityTime.SatelliteMeasurement`](@ref) row per
satellite and compiles once for every constellation mix.

```@docs
PositionVelocityTime.SignalGroup
PositionVelocityTime.SignalGroups
PositionVelocityTime.signal_groups
PositionVelocityTime.SatelliteMeasurement
PositionVelocityTime.BroadcastTimeOffset
PositionVelocityTime.collect_measurements
PositionVelocityTime.collect_measurements!
PositionVelocityTime.SUPPORTED_TIME_SYSTEMS
PositionVelocityTime.SupportedTimeSystem
```

## PVT Computation

```@docs
calc_pvt
```

### Allocation-free solving

[`calc_pvt!`](@ref) is the same solve, reusing the containers of the solution it is
handed and the buffers of a [`PVTWorkspace`](@ref), so a stream of epochs allocates
nothing.

```@docs
calc_pvt!
Base.copy(::PVTSolution)
PVTWorkspace
PositionVelocityTime.empty_keeping_capacity!
```

## Position and Velocity

```@docs
get_lla
calc_satellite_position
calc_satellite_position_and_velocity
calc_elevation_and_azimuth
```

## Atmospheric Corrections

Applied by [`calc_pvt`](@ref); documented for reference and diagnostics.

```@docs
PositionVelocityTime.select_ionospheric_parameters
PositionVelocityTime.ionospheric_delay
PositionVelocityTime.tropospheric_delay
PositionVelocityTime.saastamoinen_zenith_delays
PositionVelocityTime.niell_mapping_functions
PositionVelocityTime.KlobucharParameters
PositionVelocityTime.klobuchar_group_delay
PositionVelocityTime.BeiDouKlobucharParameters
PositionVelocityTime.beidou_klobuchar_group_delay
PositionVelocityTime.NTCMGParameters
PositionVelocityTime.BDGIMParameters
PositionVelocityTime.get_klobuchar_parameters
PositionVelocityTime.get_ntcm_g_parameters
PositionVelocityTime.get_bdgim_parameters
```

## Dilution of Precision

Read from a [`PVTSolution`](@ref)'s `dop` field, e.g. `pvt.dop.GDOP`.

```@docs
PositionVelocityTime.DOP
```

## The Measurement-Model Surface

The pieces [`calc_pvt`](@ref) is assembled from, as a contract for a consumer running its
own estimator over the same measurement model, such as a navigation filter. They are
unexported; bind them with `using PositionVelocityTime: collect_measurements, …`.

The surface is specified on the flat row: a consumer first calls
[`PositionVelocityTime.collect_measurements`](@ref), and every function below that takes
per-satellite input takes those
`Vector{`[`PositionVelocityTime.SatelliteMeasurement`](@ref)`}` rows and returns one value
per row, in order. Only
[`calc_transmit_time`](@ref PositionVelocityTime.calc_transmit_time),
[`calc_satellite_clock_drift`](@ref PositionVelocityTime.calc_satellite_clock_drift),
[`get_week`](@ref PositionVelocityTime.get_week),
[`broadcast_time_offset`](@ref PositionVelocityTime.broadcast_time_offset) and the
decoder form of [`time_offset_available`](@ref PositionVelocityTime.time_offset_available)
read a decoder, for a consumer building or amending a row itself. Times are
seconds-of-week counts on each satellite's own system unless a function says otherwise.

```@docs
PositionVelocityTime.SPEED_OF_LIGHT
PositionVelocityTime.EARTH_ROTATION_RATE
PositionVelocityTime.calc_transmit_time
PositionVelocityTime.calc_satellite_clock_drift
PositionVelocityTime.fold_week_crossover
PositionVelocityTime.calc_pseudoranges
PositionVelocityTime.calc_pseudoranges!
PositionVelocityTime.BiasColumns
PositionVelocityTime.ClockConstraint
PositionVelocityTime.decide_bias_layout
PositionVelocityTime.decide_bias_layout!
PositionVelocityTime.BiasLayout
PositionVelocityTime.fit_columns
PositionVelocityTime.dop_columns
PositionVelocityTime.select_reference_system
PositionVelocityTime.predict_pseudoranges!
PositionVelocityTime.calc_design_matrix
PositionVelocityTime.calc_dop
PositionVelocityTime.calc_dop!
PositionVelocityTime.calc_pseudorange_rate
PositionVelocityTime.estimate_velocity
PositionVelocityTime.estimate_velocity!
PositionVelocityTime.leap_second_offset_to_gpst
PositionVelocityTime.calc_leap_second_corrections
PositionVelocityTime.time_offset_available
PositionVelocityTime.broadcast_time_offset
PositionVelocityTime.calc_steering_offset
PositionVelocityTime.broadcast_inter_system_bias
PositionVelocityTime.get_week
PositionVelocityTime.system_start_epoch
PositionVelocityTime.week_start_epoch
PositionVelocityTime.day_of_year
PositionVelocityTime.predict_atmospheric_delays
PositionVelocityTime.predict_atmospheric_delays!
PositionVelocityTime.IonosphericModel
PositionVelocityTime.calc_course_over_ground
```

The unexported internals the surface above links to, documented for reference:

```@docs
PositionVelocityTime.calc_design_matrix!
PositionVelocityTime.estimate_position
PositionVelocityTime.estimate_position!
PositionVelocityTime.num_parameters
PositionVelocityTime.num_observations
PositionVelocityTime.unique_time_systems
PositionVelocityTime.time_system_index
```

## Multi-GNSS Classification

Each satellite is classified by three GNSSSignals keys of its ranging signal, carried on
every [`PositionVelocityTime.SatelliteMeasurement`](@ref): `time_system`
(`get_time_system`, e.g. `GPST()`) groups the receiver clock biases, `band_id`
(`get_band_id`, e.g. `:L5`) the inter-frequency biases, and `signal_id` (`get_signal_id`,
e.g. `:GPSL1CA`) keys [`PVTSolution`](@ref)'s `satellites`. The inter-frequency biases are
laid out from the constellation × band coverage graph:

```@docs
PositionVelocityTime.inter_frequency_bias_layout
```

## Trimmed Executables

The solver compiles into a `juliac --trim=safe` executable (Julia 1.12+ with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)); `test/trim` builds one for every
navigation-data type and a mixed GPS + Galileo + BeiDou epoch and checks it against a
regular session. Two dependencies cannot be trimmed, so the package carries trim-safe
stand-ins with their interface, [`PositionVelocityTime.TAITimes`](@ref) for AstroTime's
`TAIEpoch` and [`PositionVelocityTime.LevenbergMarquardt`](@ref) for LsqFit's
`curve_fit`.

```@docs
PositionVelocityTime.TAITimes
PositionVelocityTime.LevenbergMarquardt
PositionVelocityTime.LevenbergMarquardt.curve_fit!
PositionVelocityTime.LevenbergMarquardt.LMWorkspace
PositionVelocityTime.LevenbergMarquardt.grown
PositionVelocityTime.LevenbergMarquardt.LMResult
```
