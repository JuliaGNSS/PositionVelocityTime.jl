# API Reference

## Types

```@docs
SatelliteState
PVTSolution
SatelliteInfo
InterFrequencyBias
BroadcastInterSystemBias
```

A solution's `time` is GNSSSignals' `TAITime`, a trim-safe TAI epoch; with AstroTime
loaded, `TAIEpoch(pvt.time)` and `TAITime(epoch)` convert between it and AstroTime's
`TAIEpoch` exactly.

## Signal Groups

One epoch's measurements reach [`calc_pvt`](@ref) as *signal groups*: a `NamedTuple` of
[`PositionVelocityTime.SignalGroup`](@ref)s, one per ranging signal. Deliberately the
same names `Tracking.jl` uses, and deliberately unexported, so
`using Tracking, PositionVelocityTime` does not clash and both types print fully
qualified.

Grouping is what makes the solve type-stable: within a group every satellite shares one
concrete [`SatelliteState`](@ref) type, so the per-satellite work of
[`PositionVelocityTime.collect_measurements`](@ref) is statically dispatched. Behind that
collection pass the solver sees one flat, parameter-free
[`PositionVelocityTime.SatelliteMeasurement`](@ref) row per satellite and therefore
compiles **once** for every constellation mix.

Groups are built where the satellites are tracked, not derived from a pooled vector at
solve time — [`PositionVelocityTime.signal_groups`](@ref) builds a whole epoch's groups
from a `Tracking.TrackState`. There is deliberately no conversion from a flat
`Vector{SatelliteState}`; `calc_pvt` refuses one with an error saying what to build.

```@docs
PositionVelocityTime.SignalGroup
PositionVelocityTime.SignalGroups
PositionVelocityTime.signal_groups
PositionVelocityTime.SatelliteMeasurement
PositionVelocityTime.collect_measurements
PositionVelocityTime.collect_measurements!
PositionVelocityTime.collect_measurement_rows!
PositionVelocityTime.SUPPORTED_TIME_SYSTEMS
PositionVelocityTime.SupportedTimeSystem
```

## PVT Computation

```@docs
calc_pvt
```

### Allocation-free solving

[`calc_pvt!`](@ref) is the same solve, returning a new (immutable) [`PVTSolution`](@ref)
that reuses the containers of the solution it is handed, with its scratch in a reusable
[`PVTWorkspace`](@ref). Once both have held an epoch of a given size, a solve allocates
nothing. The overwrite is explicit: only the containers of the `solution` argument are
reused, and passing the same solution as `prev_pvt` — the receiver loop,
`pvt = calc_pvt!(pvt, workspace, groups, pvt)` — is allowed.

```@docs
calc_pvt!
Base.copy(::PVTSolution)
PVTWorkspace
PositionVelocityTime.empty_keeping_capacity!
```

## Position and Velocity

```@docs
get_lla
calc_satellite_position_and_velocity
calc_elevation_and_azimuth
```

## Atmospheric Corrections

These corrections are applied automatically by [`calc_pvt`](@ref); they are
documented here for reference and for diagnostic use. The broadcast coefficient sets the
ionospheric models take, `KlobucharParameters`, `BeiDouKlobucharParameters`,
`NeQuickGParameters` and `BDGIMParameters`, are GNSSDecoder's, read off a decoder by its
`get_klobuchar_parameters`, `get_nequick_g_parameters` and `get_bdgim_parameters`.

```@docs
PositionVelocityTime.select_ionospheric_parameters
PositionVelocityTime.ionospheric_delay
PositionVelocityTime.tropospheric_delay
PositionVelocityTime.saastamoinen_zenith_delays
PositionVelocityTime.niell_mapping_functions
PositionVelocityTime.klobuchar_group_delay
PositionVelocityTime.beidou_klobuchar_group_delay
```

## Dilution of Precision

The DOP values are read from the `dop` field of a [`PVTSolution`](@ref), e.g.
`pvt.dop.GDOP`.

```@docs
PositionVelocityTime.DOP
```

## The Measurement-Model Surface

The pieces [`calc_pvt`](@ref) is assembled from, documented so a consumer that
runs its own estimator over the same measurement model — a navigation filter
closing tracking loops, for instance — reuses this package's model as a
stated contract. The names are deliberately not exported (they are solver
internals, not part of every user's vocabulary); bind them explicitly with
`using PositionVelocityTime: collect_measurements, …`.

**The surface is specified on the measurements.** A consumer's first call is
[`PositionVelocityTime.collect_measurements`](@ref), which turns its signal groups into
a `Vector{`[`PositionVelocityTime.SatelliteMeasurement`](@ref)`}`; every function below
that takes per-satellite input takes those measurements, in that order, and returns one
value per measurement in the same order. A measurement is a single concrete type with no
type parameters, so a
consumer's own estimator compiles once over it as well, however many constellations it
mixes.

Only [`calc_transmit_time`](@ref PositionVelocityTime.calc_transmit_time) reads a
decoder, through a [`SatelliteState`](@ref), and
[`calc_satellite_clock_drift`](@ref PositionVelocityTime.calc_satellite_clock_drift) takes
GNSSDecoder's `BroadcastEphemeris`, for a consumer building or amending a measurement
itself; the decoder's week, ephemeris, group delays and broadcast time offsets are
GNSSDecoder's `get_week`, `get_broadcast_ephemeris`, `get_group_delay` and
`get_broadcast_time_offset`.

Times returned or taken here are seconds-of-week counts on each satellite's own system
scale unless a function says otherwise. The week arithmetic is GNSSSignals':
`get_week_start_epoch` dates a week, `get_leap_second_offset` relates two systems'
counts, and `fold_week_crossover` folds a difference of two counts across a week
crossover.

```@docs
PositionVelocityTime.SPEED_OF_LIGHT
PositionVelocityTime.EARTH_ROTATION_RATE
PositionVelocityTime.calc_transmit_time
PositionVelocityTime.calc_satellite_clock_drift
PositionVelocityTime.BiasColumns
PositionVelocityTime.ClockConstraint
PositionVelocityTime.decide_bias_layout
PositionVelocityTime.decide_bias_layout!
PositionVelocityTime.BiasLayout
PositionVelocityTime.fit_bias_columns
PositionVelocityTime.dop_bias_columns
PositionVelocityTime.select_reference_system
PositionVelocityTime.predict_pseudoranges!
PositionVelocityTime.calc_design_matrix
PositionVelocityTime.calc_dop
PositionVelocityTime.calc_dop!
PositionVelocityTime.calc_pseudorange_rate
PositionVelocityTime.estimate_velocity
PositionVelocityTime.estimate_velocity!
PositionVelocityTime.has_broadcast_time_offset
PositionVelocityTime.broadcast_time_offset
PositionVelocityTime.broadcast_inter_system_bias
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

When [`calc_pvt`](@ref) combines constellations and bands it classifies each satellite by
three keys, all provided by GNSSSignals (3.3+) and read from a satellite's ranging signal:
`GNSSSignals.get_time_system` (a `GNSSSignals.TimeSystem`, i.e. `GPST()`/`GST()`/`BDT()`) drives
receiver-clock grouping (one bias per time system); `GNSSSignals.get_band_id` (e.g. `:L1`,
`:L5`) drives inter-frequency-bias grouping (one bias per band); and
`GNSSSignals.get_signal_id` (e.g. `:GPSL1CA`) is the per-signal identity used in the `satellites`
key of [`PVTSolution`](@ref), and the key a
[`PositionVelocityTime.SignalGroup`](@ref) is one signal's worth of. All three are read
off the ranging signal once, by the collection pass, and carried on every
[`PositionVelocityTime.SatelliteMeasurement`](@ref) as `time_system`, `band_id` and
`signal_id`.

The receiver inter-frequency biases and their reference bands (reported per
[`InterFrequencyBias`](@ref)) are laid out from the constellation × band coverage graph:

```@docs
PositionVelocityTime.ifb_layout
```

## Trimmed Executables

The solver compiles into a `juliac --trim=safe` executable (Julia 1.12+ with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)); `test/trim` holds an app that does
so for every navigation-data type and a mixed GPS + Galileo + BeiDou epoch, and checks
its output against a regular session. Signal groups are what make that possible for any
constellation mix — each group is concretely typed, and the solver behind them sees one
concrete measurement type.

LsqFit, which the package would otherwise use, cannot be trimmed, so it carries a
trim-safe stand-in for the part it uses, a self-contained submodule with the original's
interface — to be dropped once the original can be trimmed:
[`PositionVelocityTime.LevenbergMarquardt`](@ref) replaces LsqFit's `curve_fit`, which
calls the model through abstractly typed fields. AstroTime's `TAIEpoch` cannot be
trimmed either, which is why a solution is dated in GNSSSignals' trim-safe `TAITime`.

```@docs
PositionVelocityTime.LevenbergMarquardt
PositionVelocityTime.LevenbergMarquardt.curve_fit!
PositionVelocityTime.LevenbergMarquardt.LMWorkspace
PositionVelocityTime.LevenbergMarquardt.grown
PositionVelocityTime.LevenbergMarquardt.LMResult
```
