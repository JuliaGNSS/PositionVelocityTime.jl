# The GNSS time systems: broadcast offsets between them, the offsets between their
# counts, and the TAI epochs and weeks that date them.

"""
    broadcast_time_offset(measurement::SatelliteMeasurement, target) -> Union{Nothing,Float64}

The broadcast offset (s) from `measurement`'s own GNSS time system to `target` at its
transmit time, as GNSSDecoder's `get_broadcast_time_offset` evaluated it: the steering
residual between two atomic scales, tens of nanoseconds, without the whole seconds between
their counts, which [`calc_time_scale_offsets`](@ref) applies to the transmit times.
`nothing` where none is broadcast or `target` is outside
[`CANDIDATE_HUB_SYSTEMS`](@ref). Toward GPS Time that is Galileo's GGTO and BeiDou's BGTO;
toward Galileo System Time, BeiDou's BGTO variant and the GGTO GPS broadcasts on
CNAV/CNAV-2. Such an offset lets that constellation's measurements be expressed on the
target system's clock, which makes a fix possible when the geometry is too weak to
estimate an independent clock bias for it — see [`decide_bias_layout`](@ref), which
collapses onto a hub system this way.
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
    time_scale_offset_to_gpst(time_system) -> Float64

Signed offset of a GNSS time system's *count* against GPS Time's for the same
instant — `get_tai_offset(GPST) − get_tai_offset(time_system)`, so **negative**
where the system's count reads lower (BDT: `19 − 33 = −14.0`), and `0.0` for a
system that counts alike. `calc_time_scale_offsets` adds the *difference of two
of these* to a transmit time, which is what puts the +14 s onto a BeiDou time
in a GPS-referenced solve.

This is structural, not a bias: it follows from the time scales' definitions, not
from either system's steering. Both GPST and GST are `TAI − 19 s`, so they count
identically and this is `0.0` for every GPS and Galileo satellite. BDT is
`TAI − 33 s`, so a BeiDou second-of-week reads 14 s lower than the GPS
time-of-week for the same instant — BDT week 0 second 0 *is* GPS week 1356 time
of week 14, the 14 leap seconds that had accrued between the two epochs.

Derived from `GNSSSignals.get_tai_offset` rather than tabulated, so a constellation
added later is covered without touching this.

!!! note "Why this is a measurement correction and not a time-scale shift"

    It would be tempting to fold the 14 s into `get_time_of_week` and be done. That
    would be wrong: the satellite position is propagated from the *same* transmit
    time against the message's own `t_0e`, which is on the broadcasting system's
    scale, so shifting the reported time would move the 14 s into the ephemeris —
    about 55 km of along-track error at BeiDou MEO velocities. The transmit time
    must stay on its own scale and the correction must land on the pseudorange.
"""
time_scale_offset_to_gpst(time_system::GNSSSignals.TimeSystem) =
    ustrip(s, get_tai_offset(GPST()) - get_tai_offset(time_system))

"""
    calc_time_scale_offsets(measurements, reference_system) -> Vector{Float64}

Seconds to add to each satellite's transmit time to express it in
`reference_system`'s count, so that `calc_pseudoranges` may difference them.
`0.0` for every satellite of a system that counts alike — which is every GPS and
Galileo satellite, and all of them in a single-constellation epoch.

Relative to the *reference* system rather than to GPS Time, because the reference
system's count is also what dates the reported epoch: `reference_time` comes out
of the same differencing, and is combined with the reference system's week and
`system_start_epoch`. Anchoring on GPST instead leaves a BeiDou-referenced solve
reporting a GPS-count time of week dated from the BDT epoch — a 14 s error in
`PVTSolution.time` traded for the one this removes from the pseudoranges.

This is applied to the times *handed to the differencing*, and nowhere else. The
satellite position must still be propagated from the untouched transmit time
against the message's own `t_0e`, which is on the broadcasting system's scale, so
shifting the time itself would move the offset into the ephemeris — about 55 km
of along-track error at BeiDou MEO velocities. `SatelliteInfo.transmit_time` and the
broadcast time offsets likewise keep the unconverted value.
"""
calc_time_scale_offsets(measurements, reference_system) = calc_time_scale_offsets!(
    Vector{Float64}(undef, length(measurements)), measurements, reference_system)

# `calc_time_scale_offsets` into `offsets`, resized to the satellite count.
function calc_time_scale_offsets!(offsets, measurements, reference_system)
    offsets = resize!(offsets, length(measurements))
    reference = time_scale_offset_to_gpst(reference_system)
    # Each measurement already carries its own anchor (`count_offset_to_gpst`), precomputed by
    # `collect_measurements`, so this reads a number rather than branching on the
    # `time_system` field per satellite.
    for (j, measurement) in enumerate(measurements)
        offsets[j] = reference - measurement.count_offset_to_gpst
    end
    offsets
end

"""
    system_start_epoch(signal) -> TAITime

Absolute TAI epoch of a ranging signal's GNSS time-scale origin (week 0, time of
week 0), from GNSSSignals' `get_tai_system_start_time` — the epoch already
labelled on the atomic scale (GPS `1980-01-06T00:00:19` TAI, Galileo
`1999-08-22T00:00:19` TAI, BeiDou `2006-01-01T00:00:33` TAI), so no leap-second
table is consulted. This used to go through `get_system_start_time` (a UTC
label) and AstroTime's leap-aware conversion; GNSSSignals 4.1 states the TAI
labels itself, which is also what makes the value safe to derive at
precompile time.
"""
system_start_epoch(signal) = TAITime(get_tai_system_start_time(signal))
