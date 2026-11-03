# The GNSS time systems: broadcast offsets between them, the offsets between their
# counts, and the TAI epochs and weeks that date them.

"""
    time_offset_available(measurement::SatelliteMeasurement, target) -> Bool
    time_offset_available(decoder::GNSSDecoderState, target) -> Bool

Whether a usable broadcast offset from this satellite's own GNSS time system to
`target` is available. Such an offset lets that constellation's measurements be
expressed on the target system's clock, which makes a fix possible when the geometry
is too weak to estimate an independent clock bias for it — see
[`decide_bias_layout`](@ref), which collapses onto a hub system this way. Toward GPS
Time that is Galileo's GGTO and BeiDou's BGTO; toward Galileo System Time it is
BeiDou's BGTO variant and the GGTO GPS itself broadcasts on CNAV/CNAV-2.

The measurement form reads the flag [`collect_measurements`](@ref) already flattened
into `measurement.time_offsets` (equivalently `measurement.time_offsets[i].available`),
and is `false` for any `target` outside [`CANDIDATE_HUB_SYSTEMS`](@ref), which no
message broadcasts an offset toward.

The decoder form is the underlying predicate, and delegates to GNSSDecoder's
`get_time_offset`, which screens every way a signal can fail to have one: it
broadcasts none toward `target` (GPS L1 C/A and Galileo E6-B broadcast none at all),
it does but has not decoded one yet, the decoded one names a different target system,
the ICD's "not available" sentinel is set (a `GNSS_ID` of 0; Galileo's all-ones GGTO),
or — on Galileo — the reference week has not arrived yet, since word type 10 can be
decoded before the week number is. Always `false` when `target` is the decoder's own
system: the offset from a scale to itself is not a broadcast quantity.
"""
function time_offset_available(
    measurement::SatelliteMeasurement,
    target::GNSSSignals.TimeSystem,
)
    index = hub_system_index(target)
    isnothing(index) ? false : measurement.time_offsets[index].available
end

time_offset_available(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem,
) = !isnothing(get_time_offset(decoder, target))

"""
    calc_steering_offset(offset::BroadcastTimeOffset, t) -> Float64

The broadcast *steering* offset between a constellation's time scale and the target
scale `offset` was built against, in seconds, evaluated at the own-scale time of week
`t`. Subtract it to convert a transmit time to the target system's time. Only defined
where `offset.available` (see [`time_offset_available`](@ref)).

Tens of nanoseconds: it is the residual between two atomic scales, not the whole
difference between their counts. GNSSDecoder's `GNSSTimeOffset.A_0` folds in the
*defined* whole-second offset as well — so that `t_target = t_own − Δt` holds for the
seconds — and `offset.count_anchor` takes that part back out, because
[`calc_time_scale_offsets`](@ref) has already applied it to the transmit times before
they were differenced. Both sides are the same
[`time_scale_offset_to_gpst`](@ref) expression, so the two halves compose exactly
rather than approximately.

!!! note "The subtraction costs about seven digits of the residual"

    Recovering a ~1e-8 s residual by subtracting 14 s from a `Float64` leaves roughly
    1.8e-15 s of rounding — one ULP at 14. That is 0.5 µm of range, so it is
    irrelevant to a fix, but it is why the tests here compare the steering term with
    an explicit tolerance rather than the default `≈`.
"""
function calc_steering_offset(offset::BroadcastTimeOffset, t)
    # A `NaN` `t_0` marks a message that carries no reference epoch at all (BeiDou
    # D1/D2); there the polynomial argument is the time of week itself. Where it is
    # present, `week_term` already carries `SECONDS_PER_WEEK * (own week − WN_0)`.
    Δτ = isnan(offset.t_0) ? t : t - offset.t_0 + offset.week_term
    offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2 - offset.count_anchor
end

# The broadcast offsets toward every candidate hub, flattened. `map` over the
# heterogeneous `CANDIDATE_HUB_SYSTEMS` tuple specialises per element, so each
# `get_time_offset` call is statically dispatched on both the decoder and the target.
function broadcast_time_offsets(decoder, week, count_offset_to_gpst)
    map(CANDIDATE_HUB_SYSTEMS) do target
        offset = get_time_offset(decoder, target)
        isnothing(offset) && return NO_TIME_OFFSET
        BroadcastTimeOffset(
            true,
            offset.A_0,
            offset.A_1,
            offset.A_2,
            isnothing(offset.t_0) ? NaN : Float64(offset.t_0),
            isnothing(offset.WN_0) ? 0.0 : SECONDS_PER_WEEK * (week - offset.WN_0),
            count_offset_to_gpst - time_scale_offset_to_gpst(target),
        )
    end
end

"""
    broadcast_time_offset(decoder, target; approximate_year = year(now(UTC)))
        -> BroadcastTimeOffset

The single [`BroadcastTimeOffset`](@ref) from `decoder`'s own GNSS time system toward
`target`, as [`collect_measurements`](@ref) flattens it onto every
[`SatelliteMeasurement`](@ref) — the standalone form, for a consumer evaluating a
broadcast time offset outside a fix. Unavailable (`available == false`) for a `target`
outside [`CANDIDATE_HUB_SYSTEMS`](@ref), and for any message that carries no such
offset (see [`time_offset_available`](@ref)).
"""
function broadcast_time_offset(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem;
    approximate_year::Integer = year(now(UTC)),
)
    index = hub_system_index(target)
    isnothing(index) && return NO_TIME_OFFSET
    broadcast_time_offsets(
        decoder,
        get_week(decoder; approximate_year),
        time_scale_offset_to_gpst(get_time_system(decoder)),
    )[index]
end

"""
    time_scale_offset_to_gpst(time_system) -> Float64

Signed offset of a GNSS time system's *count* against GPS Time's for the same
instant — `get_tai_offset(GPST) − get_tai_offset(time_system)`, so **negative**
where the system's count reads lower (BDT: `19 − 33 = −14.0`), and `0.0` for a
system that counts alike. `calc_time_scale_offsets` adds the *difference of two
of these* to a transmit time, which is what puts the +14 s onto a BeiDou time
in a GPS-primary solve.

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
    calc_time_scale_offsets(measurements, primary_system) -> Vector{Float64}

Seconds to add to each satellite's transmit time to express it in
`primary_system`'s count, so that `calc_pseudoranges` may difference them.
`0.0` for every satellite of a system that counts alike — which is every GPS and
Galileo satellite, and all of them in a single-constellation epoch.

Relative to the *primary* system rather than to GPS Time, because the primary
system's count is also what dates the reported epoch: `reference_time` comes out
of the same differencing, and is combined with the primary system's week and
`system_start_epoch`. Anchoring on GPST instead leaves a BeiDou-primary solve
reporting a GPS-count time of week dated from the BDT epoch — a 14 s error in
`PVTSolution.time` traded for the one this removes from the pseudoranges.

This is applied to the times *handed to the differencing*, and nowhere else. The
satellite position must still be propagated from the untouched transmit time
against the message's own `t_0e`, which is on the broadcasting system's scale, so
shifting the time itself would move the offset into the ephemeris — about 55 km
of along-track error at BeiDou MEO velocities. `SatInfo.time` and
[`calc_steering_offset`](@ref) likewise keep the unconverted value.
"""
calc_time_scale_offsets(measurements, primary_system) = calc_time_scale_offsets!(
    Vector{Float64}(undef, length(measurements)), measurements, primary_system)

# `calc_time_scale_offsets` into `offsets`, resized to the satellite count.
function calc_time_scale_offsets!(offsets, measurements, primary_system)
    offsets = resize!(offsets, length(measurements))
    primary = time_scale_offset_to_gpst(primary_system)
    # Each row already carries its own anchor (`count_offset_to_gpst`), precomputed by
    # `collect_measurements`, so this reads a number rather than branching on the
    # `time_system` field per satellite.
    for (j, measurement) in enumerate(measurements)
        offsets[j] = primary - measurement.count_offset_to_gpst
    end
    offsets
end

"""
    system_start_epoch(system) -> TAITime

Absolute TAI epoch of a ranging signal's GNSS time-scale origin (week 0, time of
week 0), from GNSSSignals' `get_tai_system_start_time` — the epoch already
labelled on the atomic scale (GPS `1980-01-06T00:00:19` TAI, Galileo
`1999-08-22T00:00:19` TAI, BeiDou `2006-01-01T00:00:33` TAI), so no leap-second
table is consulted. This used to go through `get_system_start_time` (a UTC
label) and AstroTime's leap-aware conversion; GNSSSignals 4.1 states the TAI
labels itself, which is also what makes the value safe to derive at
precompile time.
"""
system_start_epoch(system) = TAITime(get_tai_system_start_time(system))

"""
    get_week(decoder::GNSSDecoderState{<:GPSL1CAData}; approximate_year)

Return the absolute GPS week number for a GPSL1 decoder, resolving the
1024-week rollover ambiguity using `approximate_year` as a calendar
anchor.

The legacy GPS L1 C/A LNAV message broadcasts only a 10-bit week number
(0–1023) modulo 1024, so the receiver cannot determine which 1024-week
cycle the recording is in from the data alone (IS-GPS-200, §20.3.3.3).
Each cycle is ~19.6 years, so any anchor within ±9 years of the true
observation date selects the correct cycle.

GPS week 0 is 1980-01-06; cycle boundaries fall on 1999-08-22,
2019-04-07, 2038-11-21, 2058-07-08, …

For Galileo, the broadcast WN is 12 bits and does not need this
treatment in any practical operational scenario.
"""
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData};
    approximate_year::Integer = year(now(UTC)),
)
    # GPS week 0 begins 1980-01-06. Compute the integer week count from
    # there to mid-`approximate_year`, then choose the cycle base such
    # that `cycle_base + WN` is closest to that anchor.
    days_at_anchor = Date(approximate_year, 6, 30) - Date(1980, 1, 6)
    weeks_at_anchor = Dates.value(days_at_anchor) ÷ 7
    n_cycles = round(Int, (weeks_at_anchor - decoder.data.WN) / 1024)
    return n_cycles * 1024 + decoder.data.WN
end

# Every other navigation message broadcasts the absolute week number — Galileo I/NAV
# and F/NAV 12-bit, GPS CNAV (L5, L2C) and CNAV-2 (L1C) 13-bit, BeiDou D1/D2 and
# B-CNAV 13-bit — so there is no rollover to resolve, `approximate_year` is unused,
# and one method serves all of them. GPS L1 C/A above is the sole exception, and the
# reason this function takes `approximate_year` at all.
#
# The Galileo half of the union is `AbstractGalileoEphemerisData`, not the wider
# `AbstractGalileoData`: `GalileoE6BData` is Galileo data with no week number at all
# (C/NAV stamps a time of hour instead), so the wider bound would put a method on it
# that raises a `FieldError` on the field it goes looking for. That is precisely the
# distinction GNSSDecoder introduced the narrower supertype to draw.
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
