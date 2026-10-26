# The GNSS time systems: broadcast offsets between them, and the leap seconds between
# their counts.

"""
    broadcast_time_offset(measurement::SatelliteMeasurement, target) -> Union{Nothing,Float64}

The broadcast offset (s) from `measurement`'s own GNSS time system to `target` at its
transmit time, as GNSSDecoder's `get_broadcast_time_offset` evaluated it: the steering
residual between two atomic scales, tens of nanoseconds, without the whole seconds between
their counts, which `calc_leap_second_corrections!` applies to the transmit times.
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
    calc_leap_second_corrections!(corrections, measurements, reference_system)
        -> corrections

Seconds to add to each satellite's transmit time to express it in `reference_system`'s
count (GNSSSignals' `get_leap_second_offset`), into `corrections` resized to the satellite
count, so that the transmit times may be differenced; `0.0` for systems that count alike.
Applied only to the times handed to the differencing; `SatelliteInfo.transmit_time` and the
broadcast offsets keep the unconverted value.

Relative to the reference system, not GPS Time, because its count also dates the reported
epoch: anchoring on GPST would leave a BeiDou-referenced solve 14 s off in
`PVTSolution.time`.

!!! note "A measurement correction, not a time-scale shift"

    The satellite is propagated from its transmit time against the message's own `t_0e`,
    on the broadcasting system's scale, so shifting the transmit time would move the 14 s
    into the ephemeris — about 55 km of along-track error at BeiDou MEO velocities. The
    transmit time stays on its own scale and the correction lands on the pseudorange.
"""
function calc_leap_second_corrections!(corrections, measurements, reference_system)
    corrections = resize!(corrections, length(measurements))
    for (j, measurement) in enumerate(measurements)
        corrections[j] =
            ustrip(s, get_leap_second_offset(measurement.time_system, reference_system))
    end
    corrections
end
