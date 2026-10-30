# The GNSS time systems: broadcast offsets between them, and the whole seconds between
# their counts.

"""
    broadcast_time_offset(measurement::SatelliteMeasurement, target) -> Union{Nothing,Float64}

The broadcast offset (s) from `measurement`'s own GNSS time system to `target` at its
transmit time, as GNSSDecoder's `get_broadcast_time_offset` evaluated it: the steering
residual between two atomic scales, tens of nanoseconds, without the whole seconds between
their counts, which [`calc_time_scale_offset`](@ref) applies to the transmit times.
`nothing` where none is broadcast or `target` is outside
[`SUPPORTED_TIME_SYSTEMS`](@ref). Toward GPS Time that is Galileo's GGTO and BeiDou's BGTO;
toward Galileo System Time, BeiDou's BGTO variant and the GGTO GPS broadcasts on
CNAV/CNAV-2. An offset can tie two receiver clocks (see [`decide_bias_layout`](@ref)).
"""
broadcast_time_offset(measurement::SatelliteMeasurement, ::GPST) =
    measurement.time_offset_to_gpst
broadcast_time_offset(measurement::SatelliteMeasurement, ::GST) =
    measurement.time_offset_to_gst
broadcast_time_offset(measurement::SatelliteMeasurement, ::BDT) =
    measurement.time_offset_to_bdt
broadcast_time_offset(::SatelliteMeasurement, ::GNSSSignals.TimeSystem) = nothing

"""
    has_broadcast_time_offset(measurement::SatelliteMeasurement, target) -> Bool

Whether `measurement` carries a broadcast offset toward `target`
([`broadcast_time_offset`](@ref)). GNSSDecoder screens every way one can be missing: none
broadcast toward `target` (GPS L1 C/A and Galileo E6-B broadcast none), not decoded yet,
decoded toward another system, the ICD's "not available" sentinel, or on Galileo a
reference week not yet decoded.
"""
has_broadcast_time_offset(
    measurement::SatelliteMeasurement,
    target::GNSSSignals.TimeSystem,
) = !isnothing(broadcast_time_offset(measurement, target))

"""
    calc_time_scale_offset(measurement, reference_system) -> Float64

Seconds to add to `measurement`'s transmit time to express it in `reference_system`'s
count (GNSSSignals' `get_time_scale_offset`), so that transmit times of different systems
may be differenced ([`calc_pseudorange`](@ref)); `0.0` for systems that count alike. Only
`measurement.time_system` is read. Applied only to the times handed to the differencing;
`SatelliteInfo.transmit_time` and the broadcast offsets keep the unconverted value.

`calc_pvt` corrects toward the solve's reference system, not GPS Time, because its count
also dates the reported epoch: anchoring on GPST would leave a BeiDou-referenced solve 14 s
off in `PVTSolution.time`. A consumer with an epoch of its own may anchor on any system.

!!! note "A measurement correction, not a time-scale shift"

    The satellite is propagated from its transmit time against the message's own `t_0e`,
    on the broadcasting system's scale, so shifting the transmit time would move the 14 s
    into the ephemeris — about 55 km of along-track error at BeiDou MEO velocities. The
    transmit time stays on its own scale and the correction lands on the pseudorange.
"""
calc_time_scale_offset(measurement, reference_system) =
    ustrip(s, GNSSSignals.get_time_scale_offset(measurement.time_system, reference_system))

"""
    calc_pseudorange(reference_time, transmit_time, correction = 0.0) -> Float64

The pseudorange (m) of a satellite that sent `transmit_time` (s, a time of week on its own
count) against `reference_time` (s, a time of week on the reference count): how far the
transmit time, moved onto the reference count by `correction`
([`calc_time_scale_offset`](@ref)), trails the reference time, folded across a week
crossover (GNSSSignals' `fold_week_crossover`) and scaled by [`SPEED_OF_LIGHT`](@ref).

A crossover can fall between the two times: within the spread of transmit times, or, in a
mixed solve, between BeiDou times on the GPS count (`SOW + 14`, past 604800 for the first
14 s of a GPS week) and GPS times already wrapped. `reference_time` need not be a receive
time: a common one cancels into the receiver clock biases.
"""
calc_pseudorange(reference_time, transmit_time, correction = 0.0) =
    fold_week_crossover(reference_time - (correction + transmit_time)) * SPEED_OF_LIGHT
