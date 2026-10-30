# API Reference

## Types

```@docs
SatelliteState
PVTSolution
SatelliteInfo
InterFrequencyBias
BroadcastInterSystemBias
DOP
SupportedTimeSystem
```

A solution's `time` is GNSSSignals' `TAITime`, a trim-safe TAI epoch; with AstroTime
loaded, `TAIEpoch(pvt.time)` and `TAITime(epoch)` convert between it and AstroTime's
`TAIEpoch` exactly.

## Signal Groups

One epoch's measurements reach [`calc_pvt`](@ref) as *signal groups*: a `NamedTuple` of
[`PositionVelocityTime.SignalGroup`](@ref)s, one per ranging signal, mirroring
Tracking.jl's. Within a group every satellite shares one concrete type, so
[`PositionVelocityTime.collect_measurements`](@ref) dispatches statically; behind it the
solver sees one [`PositionVelocityTime.SatelliteMeasurement`](@ref) per satellite, a
single concrete type, and compiles once for every constellation mix.

```@docs
PositionVelocityTime.SignalGroup
PositionVelocityTime.signal_groups
PositionVelocityTime.SatelliteMeasurement
PositionVelocityTime.collect_measurements
PositionVelocityTime.collect_measurements!
PositionVelocityTime.SUPPORTED_TIME_SYSTEMS
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
PositionVelocityTime.calc_course_over_ground
calc_satellite_position_and_velocity
calc_elevation_and_azimuth
```

## Atmospheric Corrections

Applied by [`calc_pvt`](@ref); documented for reference and diagnostics. The broadcast
coefficient sets the ionospheric models take, `KlobucharParameters`,
`BeiDouKlobucharParameters`, `NeQuickGParameters` and `BDGIMParameters`, are GNSSDecoder's,
read off a decoder by its `get_klobuchar_parameters`, `get_nequick_g_parameters` and
`get_bdgim_parameters`.

```@docs
PositionVelocityTime.select_ionospheric_parameters
PositionVelocityTime.ionospheric_delay
PositionVelocityTime.tropospheric_delay
```

## The Measurement-Model Surface

The pieces [`calc_pvt`](@ref) is assembled from, as a contract for a consumer running its
own estimator over the same measurement model, such as a navigation filter. They are
unexported; bind them with `using PositionVelocityTime: collect_measurements, …`.

The surface is specified on the satellite measurements: a consumer first calls
[`PositionVelocityTime.collect_measurements`](@ref), and every per-satellite input and
output below follows the order of the
`Vector{`[`PositionVelocityTime.SatelliteMeasurement`](@ref)`}` it returns. The functions
that read several fields of a measurement take the measurements themselves (the bias
layout, the atmospheric delays). The others take arrays drawn from them, one entry per
measurement: the satellite positions (a `3 × N` matrix or a vector of `SVector{3}`s) and
the pseudoranges. The least-squares
problem has one row per satellite and, after them, one per
[`ClockConstraint`](@ref PositionVelocityTime.ClockConstraint), so the pseudoranges and
the design matrix carry those extra rows.
[`calc_pseudorange_rate`](@ref PositionVelocityTime.calc_pseudorange_rate) takes one
satellite's values. Only
[`calc_transmit_time`](@ref PositionVelocityTime.calc_transmit_time) and
[`calc_uncorrected_transmit_time`](@ref PositionVelocityTime.calc_uncorrected_transmit_time)
read a decoder, through a [`SatelliteState`](@ref), and
[`calc_satellite_clock_drift`](@ref PositionVelocityTime.calc_satellite_clock_drift) takes
GNSSDecoder's `BroadcastEphemeris`, for a consumer building or amending a measurement
itself; the decoder's week, time of week, ephemeris, group delays and broadcast time
offsets are GNSSDecoder's `get_broadcast_week`, `get_decoded_time_of_week`,
`get_broadcast_ephemeris`, `get_group_delay` and `get_broadcast_time_offset`. Times are seconds-of-week counts on each satellite's own
system unless a function says otherwise;
[`calc_time_scale_offset`](@ref PositionVelocityTime.calc_time_scale_offset)
moves one onto another system's count and
[`calc_pseudorange`](@ref PositionVelocityTime.calc_pseudorange) differences two across a
week crossover, and GNSSSignals' `get_week_start_epoch` dates a week.

```@docs
PositionVelocityTime.SPEED_OF_LIGHT
PositionVelocityTime.calc_transmit_time
PositionVelocityTime.calc_uncorrected_transmit_time
PositionVelocityTime.calc_satellite_clock_drift
PositionVelocityTime.calc_time_scale_offset
PositionVelocityTime.calc_pseudorange
PositionVelocityTime.BiasColumns
PositionVelocityTime.num_parameters
PositionVelocityTime.ClockConstraint
PositionVelocityTime.decide_bias_layout
PositionVelocityTime.decide_bias_layout!
PositionVelocityTime.BiasLayout
PositionVelocityTime.fit_bias_columns
PositionVelocityTime.dop_bias_columns
PositionVelocityTime.select_reference_system
PositionVelocityTime.predict_pseudoranges!
PositionVelocityTime.calc_design_matrix
PositionVelocityTime.calc_design_matrix!
PositionVelocityTime.calc_fit_covariance!
PositionVelocityTime.calc_dop
PositionVelocityTime.calc_dop!
PositionVelocityTime.calc_pseudorange_rate
PositionVelocityTime.broadcast_time_offset
PositionVelocityTime.has_broadcast_time_offset
PositionVelocityTime.broadcast_inter_system_bias
PositionVelocityTime.predict_atmospheric_delays
PositionVelocityTime.predict_atmospheric_delays!
PositionVelocityTime.IonosphericModel
```

The unexported internals the surface above links to, documented for reference:

```@docs
PositionVelocityTime.estimate_position
PositionVelocityTime.estimate_position!
PositionVelocityTime.estimate_velocity
PositionVelocityTime.estimate_velocity!
PositionVelocityTime.num_observations
```

## Multi-GNSS Classification

Each satellite is classified by three GNSSSignals keys of its ranging signal, carried on
every [`PositionVelocityTime.SatelliteMeasurement`](@ref): `time_system`
(`get_time_system`, e.g. `GPST()`) groups the receiver clock biases, `band_id`
(`get_band_id`, e.g. `:L5`) the inter-frequency biases, and `signal_id` (`get_signal_id`,
e.g. `:GPSL1CA`) keys [`PVTSolution`](@ref)'s `satellites`. The inter-frequency biases are
laid out from the constellation × band coverage graph:

```@docs
PositionVelocityTime.ifb_layout
```

## Trimmed Executables

The solver compiles into a `juliac --trim=safe` executable (Julia 1.12+ with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)); `test/trim` builds one for every
navigation-data type and a mixed GPS + Galileo + BeiDou epoch and checks it against a
regular session. LsqFit cannot be trimmed, so the package carries a trim-safe stand-in
with its interface, [`PositionVelocityTime.LevenbergMarquardt`](@ref) for its
`curve_fit`; AstroTime cannot either, which is why a solution is dated in GNSSSignals'
trim-safe `TAITime` rather than AstroTime's `TAIEpoch`.

```@docs
PositionVelocityTime.LevenbergMarquardt
PositionVelocityTime.LevenbergMarquardt.curve_fit!
PositionVelocityTime.LevenbergMarquardt.LMWorkspace
PositionVelocityTime.LevenbergMarquardt.grown
PositionVelocityTime.LevenbergMarquardt.LMResult
```
