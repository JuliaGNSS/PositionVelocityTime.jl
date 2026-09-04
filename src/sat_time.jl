"""
    fold_week_crossover(t) -> Float64

Fold a difference of two seconds-of-week counts modulo the week, into
`(-302400, 302400]`: near a week rollover the two operands sit on opposite
sides of the wrap and their raw difference is off by ±604800 s. Every time
difference in this package — and in a consumer differencing transmit times —
must pass through this.
"""
function fold_week_crossover(t)
    half_week = SECONDS_PER_WEEK / 2
    t + (t > half_week ? -2 * half_week : (t < -half_week ? 2 * half_week : 0.0))
end

function calc_uncorrected_time(state::SatelliteState)
    system = state.system
    # Time of week at the navigation-frame reference epoch, from GNSSDecoder's own
    # accessor. Most signals broadcast seconds outright, but GPS L1C-D counts
    # two-hour intervals plus 18 s TOI steps and BeiDou B1C counts hours plus 18 s
    # SOH steps, with no seconds-of-week field at all — reconstructions that belong
    # next to the ICDs that define them rather than here.
    t_tow = get_time_of_week(state.decoder)
    # Bit-count term uses the decoder's nav-message symbol rate — the rate the bits
    # were counted at — reported by `get_data_frequency` on the decoder
    # state (GNSSDecoder 3.6). It is the decoder's *data* signal, not the tracked
    # `system`: with joint same-band tracking the pseudorange may be generated on a
    # dataless pilot (e.g. GPS L1C-P, Galileo E1C) whose own data frequency is 0 Hz,
    # while the bits come from the data component the decoder ran on. The code- and
    # carrier-phase terms use the tracked ranging `system`.
    t_bits =
        state.decoder.num_bits_after_valid_syncro_sequence /
        get_data_frequency(state.decoder) * Hz
    # The code phase is shared across the same-band signals of a jointly tracked
    # satellite and can run up to a whole data symbol of chips. That whole-symbol part
    # is already counted by `t_bits` (which uses the decoder's data-symbol rate), so
    # reduce the code phase to the residual within one data symbol to avoid
    # double-counting it. A no-op for a single-signal sat, whose code phase is already
    # below one data symbol of chips.
    chips_per_symbol =
        ustrip(Hz, get_code_frequency(system)) / ustrip(Hz, get_data_frequency(state.decoder))
    t_code_phase = mod(state.code_phase, chips_per_symbol) / get_code_frequency(system) * Hz
    # `carrier_phase` is in radians (that is what `Tracking.get_carrier_phase` reports and
    # what the Tracking extension passes on), so convert to cycles before dividing by the
    # centre frequency — unlike the code phase above, which is already in chips.
    t_carrier_phase = state.carrier_phase / 2π / get_center_frequency(system) * Hz

    t_tow + t_bits + t_code_phase + t_carrier_phase
end

function calc_relativistic_correction(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - data.t_0e)
    # √A from the effective elements: the broadcast `sqrt_A` directly for LNAV/Galileo,
    # `√(A_REF + ΔA)` for CNAV/CNAV-2 (which carry no `sqrt_A` field).
    el = orbital_elements(data, decoder.constants.μ, time_from_ephemeris_reference_epoch)
    E = calc_eccentric_anomaly(decoder, t)
    decoder.constants.F * data.e * el.sqrt_A * sin(E)
end

function correct_clock(decoder::GNSSDecoder.GNSSDecoderState, system, t)
    Δtr = calc_relativistic_correction(decoder, t)
    data = decoder.data
    # Folded modulo the week like every other time difference here: near a week
    # rollover `t` and `t_0c` sit on opposite sides of the wrap, and the raw
    # difference of ±604800 s puts ~a_f1·604800 ≈ microseconds (kilometres of
    # range) into a polynomial whose real argument is seconds.
    Δt_from_reference = fold_week_crossover(t - data.t_0c)
    Δt =
        data.a_f0 +
        data.a_f1 * Δt_from_reference +
        data.a_f2 * Δt_from_reference^2 +
        Δtr
    t - correct_by_group_delay(decoder, system, Δt)
end

"""
    calc_satellite_clock_drift(decoder::GNSSDecoderState, t) -> Float64

Rate of the satellite clock correction (s/s): the time derivative of the clock
polynomial `correct_clock` applies, so it must be evaluated about the same clock
reference epoch `t_0c`. Using `t` instead would offset the `a_f2` term by a
spurious `2·a_f2·t_0c` — zero whenever `a_f2 = 0` (the usual GPS broadcast) but
up to ~4e-9 s/s (≈1.3 m/s of range rate) at the edge of the broadcast range late
in the week, and per-satellite, so it would leak into the estimated velocity
rather than the common clock drift.

The rate of the relativistic periodic term `Δt_rel = F·e·√A·sin(E)` that
`correct_clock` includes is not modelled here; it is at most ~1 mm/s.
"""
function calc_satellite_clock_drift(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    data.a_f1 + 2 * data.a_f2 * fold_week_crossover(t - data.t_0c)
end

# Group-delay / inter-signal correction, selected by the *ranging* signal `system`
# while the values come from the `decoder`'s navigation message. The two can differ:
# with joint same-band tracking a band's pseudorange may be generated on a pilot or
# secondary signal (e.g. GPS L1C-P) while the ephemeris/clock are decoded from the
# data component (L1C-D). For GPS this selects the right ISC per signal
# (IS-GPS-705J §20.3.3.3.1.2 / IS-GPS-800J §3.5.4.1); for Galileo the broadcast group
# delay is per band, so the decoder's message alone determines it.

# GPS L1 C/A via LNAV: only T_GD (LNAV carries no ISCs).
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData},
    ::GNSSSignals.GPSL1CA,
    t,
) = t - decoder.data.T_GD

# The CNAV / CNAV-2 group-delay terms are ~metre-level inter-signal corrections a
# satellite may not have broadcast yet when its ephemeris and clock are already
# complete, so any of them can be `nothing` on a usable decoder — apply each when
# present, treat a missing one as zero, never throw.
#
# T_GD is a single per-SV group-delay differential (L1 P(Y)−L2 P(Y), IS-GPS-200
# §20.3.3.3.3.2): the same value on every GPS signal that carries it — only the ISC
# differs per signal. A missing T_GD could therefore be sourced from another GPS
# decoder for the same PRN, ideally an L1 one where it is broadcast far more often
# (L1 C/A every 30 s subframe, L1C every 18 s subframe 2, versus up to 288 s / 144 s
# on L2C / L5). Decoders are independent per (signal, PRN) here, so for now a missing
# term is simply taken as zero — a ~metre bias, never a blocked fix.
group_delay_term(x) = something(x, 0.0)

# GPS CNAV rides on one shared data container for L5 and L2C; the ISC is picked by
# the ranging signal — L5 I5 for an L5 range, L2C for an L2C-M range
# (IS-GPS-705J §20.3.3.3.1.2 / IS-GPS-200N §30.3.3.3.1.1). -T_GD + ISC.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::GNSSSignals.GPSL5I,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L5I5)
# The dataless pilots range off the same CNAV message: the L5 pilot has its own
# ISC_L5Q5, while L2 CL shares the single broadcast ISC_L2C with CM (IS-GPS-200N
# §30.3.3.3.1.1 defines one correction for the L2C signal).
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::GNSSSignals.GPSL5Q,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L5Q5)
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::Union{GNSSSignals.GPSL2CM,GNSSSignals.GPSL2CL},
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L2C)

# GPS CNAV-2 (L1C) carries the full L1 ISC set, so one decoder serves a range
# generated on the L1C data (L1C-D), the L1C pilot (L1C-P) or C/A — pick the ISC of
# the signal the range was actually generated on.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1C_D,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L1CD)
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1C_P,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L1CP)
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1CA,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L1CA)

# Galileo: which broadcast group delay applies is a property of the *band* the range
# was generated on, not of the data/pilot split — E1B and E1C share one correction,
# E5a-I and E5a-Q another, E5b-I and E5b-Q a third.
#
# Two facts drive the rules below (Galileo OS SIS ICD, Issue 2.2, §5.1.5).
#
# First, *which* BGD: the broadcast clock polynomial is referred to an ionosphere-free
# dual-frequency combination, and which one depends on the message. I/NAV (E1-B and
# E5b-I) is referred to E1/E5b, so it pairs with BGD(E1,E5b); F/NAV (E5a-I) is referred
# to E1/E5a and pairs with BGD(E1,E5a). I/NAV broadcasts *both* BGDs, but only its own
# is ever wanted — see the note below the scaling rule.
#
# Second, *how much* of it: a single-frequency user on the combination's first band
# (E1) applies the BGD as broadcast, but one on its second band (E5a or E5b) applies it
# scaled by `(f_E1/f_band)²`. That factor is 1.79 on E5a and 1.70 on E5b — dropping it
# leaves ~0.8 of a BGD, a few nanoseconds, so metre-level and per-satellite rather than
# a common clock offset the solve would absorb.
#
# I/NAV broadcasts *both* BGDs, but only ever needs its own: an E5a range carries an
# E5a (F/NAV) decoder, never an I/NAV one.
#
# `galileo_group_delay_scaling` is that factor, derived from the carrier frequencies
# rather than tabulated, so it is exactly 1 on E1 without a special case.
#
# Note the absence of `group_delay_term` below. Both Galileo messages gate their BGDs
# in `is_decoding_completed_for_positioning` — I/NAV on both of them, F/NAV on its
# only one — so a decoder that reached this point cannot carry a `nothing` here, and
# treating one as zero would hide a decoder contract violation as a few-nanosecond
# bias. `group_delay_term` is reserved for the terms a message may legitimately not
# have broadcast yet, which for Galileo is none of them.
galileo_group_delay_scaling(system) =
    (get_center_frequency(GalileoE1B) / get_center_frequency(system))^2

# One method per (message, ranging signal) that can actually occur, and no fallback —
# the same shape as the GPS and BeiDou methods above, so the whole table can be
# checked by introspection (see `test/galileo_e5b_e6b.jl`).
#
# A range is generated on the band its own data component was decoded from, so the
# ranging band is fixed by the decoder: I/NAV is decoded on E1-B or E5b-I and can be
# asked for E1 or E5b signals, F/NAV only on E5a-I and so only for E5a signals. Both
# I/NAV bands take the E1/E5b BGD its clock is referred to; F/NAV's takes E1/E5a.
#
# The absent combinations matter as much as the present ones. There is no fallback, so
# an E6 range — E6 carries no BGD at all, so there is no right answer to give — raises
# instead of returning a plausible-looking number scaled from the E5b BGD, and the
# cross-band I/NAV + E5a and F/NAV + E1 pairings raise for the same reason. I/NAV does
# broadcast BGD(E1,E5a), but nothing reads it: an E5a range carries an F/NAV decoder.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoINAVData},
    system::Union{
        GNSSSignals.GalileoE1B,
        GNSSSignals.GalileoE1C,
        GNSSSignals.GalileoE1B_BOC11,
        GNSSSignals.GalileoE1C_BOC11,
        GNSSSignals.GalileoE5bI,
        GNSSSignals.GalileoE5bQ,
    },
    t,
) =
    t -
    galileo_group_delay_scaling(system) * decoder.data.BGD_E1_E5b

correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE5aData},
    system::Union{
        GNSSSignals.GalileoE5aI,
        GNSSSignals.GalileoE5aQ,
        GNSSSignals.GalileoE5aQP,
    },
    t,
) =
    t -
    galileo_group_delay_scaling(system) * decoder.data.BGD_E1_E5a

# BeiDou references its broadcast clock polynomial to the B3I signal, on both the
# legacy and the BDS-3 messages, and then publishes a group delay per signal to
# get from there to the signal a range was actually generated on
# (BDS-SIS-ICD-B1I-3.0 §5.2.4.10, BDS-SIS-ICD-B1C-1.0 §7.6 and its B2a/B2b
# counterparts). B3I itself therefore needs no correction at all — the one
# ranging signal in this package whose clock is already referred to it. The data
# components carry an additional inter-signal correction on top of their band's
# pilot group delay (`ISC_B1Cd`, `ISC_B2ad`), exactly like the GPS L1C and L5
# ISCs above.
#
# `group_delay_term` appears on B2a's terms and on no others, which is not an
# oversight but the decoder's gating read back. `is_decoding_completed_for_positioning`
# requires a signal's own single-band correction wherever it rides in the block the
# gate already waits for, and excludes it where requiring it would mean waiting on a
# message the ICD does not schedule:
#
#   - D1/D2 gates `T_GD1`, the correction for the one band it is decoded on that has
#     a ranging signal here (not `T_GD2` — see the B1I method below); B1C gates
#     `T_GD_B1Cp` and `ISC_B1Cd` (bits
#     546-569 of the same CRC-protected subframe 2 as the ephemeris); B2b gates
#     `T_GD_B2bI` (same MT30 as the clock). All read unwrapped — a `nothing` there
#     means a broken decoder, and taking it as zero would hide that as a metre of
#     range.
#   - B2a alone is excluded: `T_GD_B2ap`, `ISC_B2ad` and `T_GD_B1Cp` ride in MT30
#     only, while the clock and IODC are in all of MT30-34, and BDS-SIS-ICD-B2a-1.0
#     §6.2 schedules only the MT10/11 pair. Gating on MT30 would gate on an interval
#     nothing bounds, so those really can be absent and are taken as zero.

# Legacy D1/D2 on B1I: T_GD1 is the B1I group delay differential. The message also
# carries T_GD2, for BDS-2's B2I, and nothing reads it: B2I has no signal definition
# in GNSSSignals and will not get one (JuliaGNSS/GNSSSignals.jl#156 — 13 of the 15
# remaining BDS-2 satellites were decommissioned in April 2026, leaving two IGSO
# transmitters and no successor, since B2a and B2b replaced it on BDS-3). So this is
# a field with no consumer rather than one awaiting a ranging signal.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
    ::GNSSSignals.BeiDouB1I,
    t,
) = t - decoder.data.T_GD1
# ... and on B3I: none, the clock is already B3I-referenced.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
    ::GNSSSignals.BeiDouB3I,
    t,
) = t

# B-CNAV1 (B1C) carries the B1C pilot group delay plus the data-component ISC. It
# also broadcasts the B2a pilot's group delay, but no method reads it: a B2a range is
# generated on the B2a component and so carries a B2a decoder, never a B1C one.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB1CData},
    ::GNSSSignals.BeiDouB1C_D,
    t,
) = t - decoder.data.T_GD_B1Cp - decoder.data.ISC_B1Cd
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB1CData},
    ::GNSSSignals.BeiDouB1C_P,
    t,
) = t - decoder.data.T_GD_B1Cp
# B-CNAV2 (B2a) mirrors it: the B2a pilot group delay plus the B2a data ISC. It too
# carries the other band's pilot delay (`T_GD_B1Cp`) that nothing here reads.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2aData},
    ::GNSSSignals.BeiDouB2aI,
    t,
) =
    t - group_delay_term(decoder.data.T_GD_B2ap) -
    group_delay_term(decoder.data.ISC_B2ad)
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2aData},
    ::GNSSSignals.BeiDouB2aQ,
    t,
) = t - group_delay_term(decoder.data.T_GD_B2ap)
# B-CNAV3 (B2b) carries only its own: B2b_I has no pilot and no sibling ISC.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2bData},
    ::GNSSSignals.BeiDouB2bI,
    t,
) = t - decoder.data.T_GD_B2bI

# There is deliberately no catch-all method, and the methods above — GPS, Galileo and
# BeiDou alike — cover exactly the ranging signals of the band each message is decoded
# on, data component and pilot. Anything else is a cross-band pairing, which cannot
# arise: a range is generated on the band its own data component was decoded from.
# Returning `t` for such a call, as an earlier version did, is the one answer that
# cannot be told apart from a correct zero correction; a `MethodError` says what
# actually happened. (D1/D2 spans two bands because one message is broadcast on both
# B1I and B3I, so both are same-band pairings with their own decoder — as are GPS
# CNAV's L5 and L2C signals, for the same reason.)

"""
    calc_corrected_time(state::SatelliteState) -> Float64

The satellite's transmit time as a seconds-of-week count on its own system's
scale: the decoder's time of week plus the bit, code-phase and carrier-phase
fractions (`calc_uncorrected_time`), corrected by the SV clock polynomial, the
relativistic term and the ranging signal's group delay (`correct_clock`). This
is the `t` every ephemeris and clock evaluation takes; only pseudorange
*differencing* converts it to another scale (see
[`time_scale_offset_to_gpst`](@ref)), never the value itself.
"""
function calc_corrected_time(state::SatelliteState)
    approximated_time = calc_uncorrected_time(state)
    correct_clock(state.decoder, state.system, approximated_time)
end

"""
    time_offset_available(decoder, target::GNSSSignals.TimeSystem) -> Bool

Whether `decoder` carries a usable broadcast offset from its own GNSS time system
to `target`. Such an offset lets that constellation's measurements be expressed on
the target system's clock, which makes a fix possible when the geometry is too
weak to estimate an independent clock bias for it — see [`decide_bias_layout`](@ref),
which collapses onto a hub system this way. Toward GPS Time that is Galileo's GGTO
and BeiDou's BGTO; toward Galileo System Time it is BeiDou's BGTO variant and the
GGTO GPS itself broadcasts on CNAV/CNAV-2.

Delegates to GNSSDecoder's `get_time_offset`, which screens every way a signal can
fail to have one: it broadcasts none toward `target` (GPS L1 C/A and Galileo E6-B
broadcast none at all), it does but has not decoded one yet, the decoded one names
a different target system, the ICD's "not available" sentinel is set (a `GNSS_ID`
of 0; Galileo's all-ones GGTO), or — on Galileo — the reference week has not
arrived yet, since word type 10 can be decoded before the week number is. Always
`false` when `target` is the decoder's own system: the offset from a scale to
itself is not a broadcast quantity.
"""
time_offset_available(decoder::GNSSDecoder.GNSSDecoderState, target::GNSSSignals.TimeSystem) =
    !isnothing(get_time_offset(decoder, target))

"""
    calc_steering_offset(decoder, target::GNSSSignals.TimeSystem, t) -> Float64

The broadcast *steering* offset between this constellation's time scale and
`target`'s, in seconds, at its own time of week `t`. Subtract it to convert a
transmit time to the target system's time. Only defined where
[`time_offset_available`](@ref) is `true` for the same `target`.

Tens of nanoseconds: it is the residual between two atomic scales, not the
whole difference between their counts. GNSSDecoder's `GNSSTimeOffset.A_0` folds
in the *defined* whole-second offset as well — so that `t_target = t_own − Δt`
holds for the seconds — and this takes that part back out (as the difference of
two [`time_scale_offset_to_gpst`](@ref) anchors, which is target-independent),
because [`calc_time_scale_offsets`](@ref) has already applied it to the transmit
times before they were differenced. Both use the same `get_tai_offset`
expression, so the two halves compose exactly rather than approximately.

!!! note "The subtraction costs about seven digits of the residual"

    Recovering a ~1e-8 s residual by subtracting 14 s from a `Float64` leaves
    roughly 1.8e-15 s of rounding — one ULP at 14. That is 0.5 µm of range, so it
    is irrelevant to a fix, but it is why the tests here compare the steering term
    with an explicit tolerance rather than the default `≈`.
"""
function calc_steering_offset(
    decoder::GNSSDecoder.GNSSDecoderState,
    target::GNSSSignals.TimeSystem,
    t,
)
    offset = get_time_offset(decoder, target)
    # `t_0`/`WN_0` are absent on BeiDou D1/D2, whose two-term offset carries no
    # reference epoch at all; there Δτ is the time of week itself. Where they are
    # present, `WN_0` has already been lifted into the decoder's own week numbering
    # by GNSSDecoder, so this subtraction is a real week count on every signal —
    # Galileo broadcasts the field in 6 bits against a 12-bit week number.
    Δτ =
        isnothing(offset.t_0) ? t :
        t - offset.t_0 + SECONDS_PER_WEEK * (get_week(decoder) - offset.WN_0)
    total = offset.A_0 + offset.A_1 * Δτ + offset.A_2 * Δτ^2
    total - (
        time_scale_offset_to_gpst(get_time_system(decoder)) -
        time_scale_offset_to_gpst(target)
    )
end
