# The GNSS time systems: broadcast offsets between them, the leap seconds between their
# counts, and the TAI epochs and weeks that date them.

"""
    time_offset_available(measurement::SatelliteMeasurement, target) -> Bool
    time_offset_available(decoder::GNSSDecoderState, target) -> Bool

Whether a usable broadcast offset from this satellite's own GNSS time system to `target`
is available, which can tie the two clocks (see [`decide_bias_layout`](@ref)). Toward GPS
Time that is Galileo's GGTO and BeiDou's BGTO; toward Galileo System Time, BeiDou's BGTO
variant and the GGTO GPS broadcasts on CNAV/CNAV-2.

The measurement form reads the flag [`collect_measurements`](@ref) flattened into
`measurement.broadcast_time_offsets`; it is `false` for any `target` outside
[`SUPPORTED_TIME_SYSTEMS`](@ref).

The decoder form delegates to GNSSDecoder's `get_time_offset`, which screens every way an
offset can be missing: none broadcast toward `target` (GPS L1 C/A and Galileo E6-B
broadcast none), not decoded yet, decoded toward another system, the ICD's "not
available" sentinel (a `GNSS_ID` of 0; Galileo's all-ones GGTO), or on Galileo a
reference week not yet decoded. Always `false` when `target` is the decoder's own system.
"""
function time_offset_available(
    measurement::SatelliteMeasurement,
    target::GNSSSignals.TimeSystem,
)
    index = time_offset_index(target)
    isnothing(index) ? false : measurement.broadcast_time_offsets[index].available
end

time_offset_available(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem,
) = !isnothing(get_time_offset(decoder, target))

"""
    calc_steering_offset(offset::BroadcastTimeOffset, t) -> Float64

The broadcast *steering* offset (s) between a constellation's time scale and the target
scale of `offset`, at the own-scale time of week `t`; subtract it to convert a transmit
time to the target system's time. Only defined where `offset.available`.

Tens of nanoseconds: the residual between two atomic scales. GNSSDecoder's
`GNSSTimeOffset.A_0` also folds in the defined whole-second offset, which
`offset.leap_second_offset_to_target` takes back out, as
[`calc_leap_second_corrections`](@ref) already applies it to the transmit times. Both use
[`leap_second_offset_to_gpst`](@ref), so the halves cancel exactly. Subtracting 14 s
costs about 1.8e-15 s of rounding (0.5 µm of range), hence the tests' explicit tolerance.
"""
function calc_steering_offset(offset::BroadcastTimeOffset, t)
    # A `NaN` `t_0` marks a message without a reference epoch (BeiDou D1/D2), whose
    # polynomial argument is the time of week itself. Not folded: `week_term` adds the
    # whole weeks since the reference week, `SECONDS_PER_WEEK * (own week − WN_0)`.
    Δτ = isnan(offset.t_0) ? t : t - offset.t_0 + offset.week_term
    offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2 -
    offset.leap_second_offset_to_target
end

# The broadcast offsets toward every supported time system, flattened. `map` over the
# tuple specialises per element, so each `get_time_offset` call dispatches statically.
function broadcast_time_offsets(decoder, week, own_leap_second_offset)
    map(SUPPORTED_TIME_SYSTEMS) do target
        offset = get_time_offset(decoder, target)
        isnothing(offset) && return NO_TIME_OFFSET
        BroadcastTimeOffset(
            true,
            offset.A_0,
            offset.A_1,
            offset.A_2,
            isnothing(offset.t_0) ? NaN : Float64(offset.t_0),
            isnothing(offset.WN_0) ? 0.0 : SECONDS_PER_WEEK * (week - offset.WN_0),
            own_leap_second_offset - leap_second_offset_to_gpst(target),
        )
    end
end

"""
    broadcast_time_offset(decoder, target; approximate_year = year(now(UTC)))
        -> BroadcastTimeOffset

The [`BroadcastTimeOffset`](@ref) from `decoder`'s own GNSS time system toward `target`,
as [`collect_measurements`](@ref) flattens it onto every row, for use outside a fix.
Unavailable for a `target` outside [`SUPPORTED_TIME_SYSTEMS`](@ref) or a message without
such an offset.
"""
function broadcast_time_offset(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem;
    approximate_year::Integer = year(now(UTC)),
)
    index = time_offset_index(target)
    isnothing(index) && return NO_TIME_OFFSET
    broadcast_time_offsets(
        decoder,
        get_week(decoder; approximate_year),
        leap_second_offset_to_gpst(get_time_system(decoder)),
    )[index]
end

"""
    leap_second_offset_to_gpst(time_system) -> Float64

The signed leap-second offset of a GNSS time system's *count* against GPS Time's for the
same instant: `get_tai_offset(GPST) − get_tai_offset(time_system)`, negative where the
system's count reads lower (BDT: `19 − 33 = −14.0`). GPST and GST are both `TAI − 19 s`,
so it is `0.0` for every GPS and Galileo satellite; BDT week 0 second 0 *is* GPS week 1356
time of week 14. Structural, not a bias, and derived from `GNSSSignals.get_tai_offset`, so
a later constellation is covered.

!!! note "A measurement correction, not a time-scale shift"

    The satellite is propagated from its transmit time against the message's own `t_0e`,
    on the broadcasting system's scale, so shifting the transmit time would move the 14 s
    into the ephemeris — about 55 km of along-track error at BeiDou MEO velocities. The
    transmit time stays on its own scale and the correction lands on the pseudorange
    ([`calc_leap_second_corrections`](@ref)).
"""
leap_second_offset_to_gpst(time_system::GNSSSignals.TimeSystem) =
    ustrip(s, get_tai_offset(GPST()) - get_tai_offset(time_system))

"""
    calc_leap_second_corrections(measurements, reference_system) -> Vector{Float64}

Seconds to add to each satellite's transmit time to express it in `reference_system`'s
count, so that [`calc_pseudoranges`](@ref) may difference them; `0.0` for systems that
count alike. Applied only to the times handed to the differencing (see
[`leap_second_offset_to_gpst`](@ref)); `SatelliteInfo.transmit_time` and
[`calc_steering_offset`](@ref) keep the unconverted value.

Relative to the reference system, not GPS Time, because its count also dates the reported
epoch: anchoring on GPST would leave a BeiDou-referenced solve 14 s off in
`PVTSolution.time`.
"""
calc_leap_second_corrections(measurements, reference_system) =
    calc_leap_second_corrections!(Vector{Float64}(undef, length(measurements)),
        measurements, reference_system)

# `calc_leap_second_corrections` into `corrections`, resized to the satellite count.
function calc_leap_second_corrections!(corrections, measurements, reference_system)
    corrections = resize!(corrections, length(measurements))
    reference = leap_second_offset_to_gpst(reference_system)
    # The row's precomputed offset, rather than a branch on `time_system` per satellite.
    for (j, measurement) in enumerate(measurements)
        corrections[j] = reference - measurement.leap_second_offset_to_gpst
    end
    corrections
end

"""
    system_start_epoch(signal) -> TAITime

The TAI epoch of a ranging signal's GNSS time-scale origin (week 0, time of week 0), from
GNSSSignals' `get_tai_system_start_time` (GPS `1980-01-06T00:00:19` TAI, Galileo
`1999-08-22T00:00:19` TAI, BeiDou `2006-01-01T00:00:33` TAI). No leap-second table is
consulted, so it is safe to derive at precompile time.
"""
system_start_epoch(signal) = TAITime(get_tai_system_start_time(signal))

"""
    week_start_epoch(signal_or_time_system, week) -> TAITime
    week_start_epoch(measurement::SatelliteMeasurement) -> TAITime

The TAI epoch at which `week` of a GNSS time system begins: [`system_start_epoch`](@ref)
plus `week` whole weeks, exactly. A time of week `t` dates as `week_start_epoch(…) + t`
with full `Float64` resolution in the fraction. The measurement form reads the row's own
start epoch and week.
"""
week_start_epoch(signal_or_time_system, week) =
    system_start_epoch(signal_or_time_system) + week * SECONDS_PER_WEEK
week_start_epoch(measurement::SatelliteMeasurement) =
    measurement.system_start_epoch + measurement.week * SECONDS_PER_WEEK

"""
    get_week(decoder::GNSSDecoderState{<:GPSL1CAData}; approximate_year)

The absolute GPS week number of a GPS L1 C/A decoder. LNAV broadcasts only a 10-bit week
number (IS-GPS-200, §20.3.3.3), so the 1024-week rollover is resolved against
`approximate_year`: a cycle is ~19.6 years, so any year within ±9 years of the
observation selects the right one (boundaries 1999-08-22, 2019-04-07, 2038-11-21, …).
"""
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData};
    approximate_year::Integer = year(now(UTC)),
)
    # The cycle that puts `WN` closest to mid-`approximate_year` (GPS week 0: 1980-01-06).
    days_at_anchor = Date(approximate_year, 6, 30) - Date(1980, 1, 6)
    weeks_at_anchor = Dates.value(days_at_anchor) ÷ 7
    n_cycles = round(Int, (weeks_at_anchor - decoder.data.WN) / 1024)
    return n_cycles * 1024 + decoder.data.WN
end

# Every other message broadcasts the absolute week number (Galileo I/NAV and F/NAV 12-bit;
# GPS CNAV, CNAV-2, BeiDou D1/D2 and B-CNAV 13-bit), so `approximate_year` is unused.
# `AbstractGalileoEphemerisData` rather than `AbstractGalileoData`, because
# `GalileoE6BData` carries no week number (C/NAV stamps a time of hour).
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{
        <:Union{
            GNSSDecoder.AbstractGalileoEphemerisData,
            GNSSDecoder.AbstractBeiDouData,
            AbstractGPSCNAVData,
        },
    };
    approximate_year::Integer = year(now(UTC)),
)
    decoder.data.WN
end
