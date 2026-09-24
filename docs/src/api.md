# API Reference

## Types

```@docs
SatelliteState
PVTSolution
TAITime
SatInfo
InterFrequencyBias
```

With AstroTime loaded, a package extension converts between [`TAITime`](@ref) and
AstroTime's `TAIEpoch` both ways — `TAIEpoch(pvt.time)` and `TAITime(epoch)` — exactly,
as the two share their layout.

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
PositionVelocityTime.BroadcastTimeOffset
PositionVelocityTime.collect_measurements
PositionVelocityTime.collect_measurements!
PositionVelocityTime.CANDIDATE_HUB_SYSTEMS
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
PVTWorkspace
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
`using PositionVelocityTime: collect_measurements, …`.

**The surface is specified on the flat measurement row.** A consumer's first call is
[`PositionVelocityTime.collect_measurements`](@ref), which turns its signal groups into
a `Vector{`[`PositionVelocityTime.SatelliteMeasurement`](@ref)`}`; every function below
that takes per-satellite input takes those rows, in that order, and returns one value
per row in the same order. A row is a single concrete type with no type parameters, so a
consumer's own estimator compiles once over it as well, however many constellations it
mixes.

Only [`calc_corrected_time`](@ref PositionVelocityTime.calc_corrected_time),
[`calc_satellite_clock_drift`](@ref PositionVelocityTime.calc_satellite_clock_drift),
[`get_week`](@ref PositionVelocityTime.get_week) and
[`broadcast_time_offset`](@ref PositionVelocityTime.broadcast_time_offset) still read a
decoder — they are what the collection pass is built from, and are listed so a consumer
can build or amend a row itself.

Times returned or taken here are seconds-of-week counts on each satellite's own system
scale unless a function says otherwise.

```@docs
PositionVelocityTime.SPEED_OF_LIGHT
PositionVelocityTime.calc_corrected_time
PositionVelocityTime.calc_satellite_clock_drift
PositionVelocityTime.get_sat_position
PositionVelocityTime.get_sat_velocity
PositionVelocityTime.fold_week_crossover
PositionVelocityTime.BiasColumns
PositionVelocityTime.decide_bias_layout
PositionVelocityTime.decide_bias_layout!
PositionVelocityTime.BiasLayoutWorkspace
PositionVelocityTime.bias_layout
PositionVelocityTime.BiasLayout
PositionVelocityTime.calc_ρ_hat!
PositionVelocityTime.calc_H
PositionVelocityTime.calc_line_of_sight
PositionVelocityTime.calc_DOP
PositionVelocityTime.calc_DOP!
PositionVelocityTime.calc_user_velocity_and_clock_drift
PositionVelocityTime.calc_user_velocity_and_clock_drift!
PositionVelocityTime.time_scale_offset_to_gpst
PositionVelocityTime.calc_time_scale_offsets
PositionVelocityTime.time_offset_available
PositionVelocityTime.broadcast_time_offset
PositionVelocityTime.calc_steering_offset
PositionVelocityTime.get_week
PositionVelocityTime.system_start_epoch
PositionVelocityTime.day_of_year
PositionVelocityTime.predict_atmospheric_delays
PositionVelocityTime.predict_atmospheric_delays!
PositionVelocityTime.IonosphericModel
PositionVelocityTime.calc_course_over_ground
```

The unexported internals the surface above links to, documented for reference:

```@docs
PositionVelocityTime.calc_H!
PositionVelocityTime.user_position
PositionVelocityTime.user_position!
PositionVelocityTime.num_lsq_params
PositionVelocityTime.calc_hub_range_offsets
PositionVelocityTime.positive_definite_cholesky
PositionVelocityTime.unique_time_systems
PositionVelocityTime.time_system_index
```

## Multi-GNSS Classification

When [`calc_pvt`](@ref) combines constellations and bands it classifies each satellite by
three keys, all provided by GNSSSignals (3.3+) and read from a satellite's ranging signal:
`GNSSSignals.get_time_system` (a `GNSSSignals.TimeSystem`, i.e. `GPST()`/`GST()`/`BDT()`) drives
receiver-clock grouping (one bias per time system); `GNSSSignals.get_band_id` (e.g. `:L1`,
`:L5`) drives inter-frequency-bias grouping (one bias per band); and
`GNSSSignals.get_signal_id` (e.g. `:GPSL1CA`) is the per-signal identity used in the `sats`
key of [`PVTSolution`](@ref), and the key a
[`PositionVelocityTime.SignalGroup`](@ref) is one signal's worth of. All three are read
off the ranging signal once, by the collection pass, and carried on every
[`PositionVelocityTime.SatelliteMeasurement`](@ref) as `time_system`, `band_id` and
`signal_id`.

The receiver inter-frequency biases and their reference bands (reported per
[`InterFrequencyBias`](@ref)) are laid out from the constellation × band coverage graph:

```@docs
PositionVelocityTime.band_ifb_layout
PositionVelocityTime.band_ifb_layout!
```

## Trimmed Executables

The solver compiles into a `juliac --trim=safe` executable (Julia 1.12+ with
[JuliaC](https://github.com/JuliaLang/JuliaC.jl)); `test/trim` holds an app that does
so for every navigation-data type and a mixed GPS + Galileo + BeiDou epoch, and checks
its output against a regular session. Signal groups are what make that possible for any
constellation mix — each group is concretely typed, and the solver behind them sees one
concrete row type.

Two dependencies the package would otherwise have cannot be trimmed, so it carries
trim-safe stand-ins for the part of each it uses, as self-contained submodules with the
original's interface — to be dropped once the originals can be trimmed:

- [`PositionVelocityTime.TAITimes`](@ref) replaces AstroTime's `TAIEpoch` as the type
  of `PVTSolution.time` (AstroTime's own `__init__` cannot be trimmed), and
- [`PositionVelocityTime.LevenbergMarquardt`](@ref) replaces LsqFit's `curve_fit`,
  which calls the model through abstractly typed fields.

```@docs
PositionVelocityTime.TAITimes
PositionVelocityTime.LevenbergMarquardt
PositionVelocityTime.LevenbergMarquardt.curve_fit
PositionVelocityTime.LevenbergMarquardt.curve_fit!
PositionVelocityTime.LevenbergMarquardt.LMWorkspace
PositionVelocityTime.LevenbergMarquardt.grown
PositionVelocityTime.LevenbergMarquardt.LMResult
```
