"""
    fold_week_crossover(t) -> Float64

Fold a difference of two seconds-of-week counts into `(-302400, 302400]`: across a week
crossover (IS-GPS-200's term) the raw difference is off by ±604800 s. Every difference of
two counts within a week of each other, here and in a consumer, passes through this; one
across whole weeks must not ([`calc_steering_offset`](@ref) adds the weeks back itself).
"""
function fold_week_crossover(t)
    half_week = SECONDS_PER_WEEK / 2
    t + (t > half_week ? -2 * half_week : (t < -half_week ? 2 * half_week : 0.0))
end

"""
    calc_pseudoranges(times) -> (pseudoranges, reference_time)

The pseudoranges (m) of transmit `times` received at one instant: how far each trails the
latest of them, `reference_time`, in metres and folded across a week crossover
([`fold_week_crossover`](@ref)). The common reference cancels into the receiver clock
biases, which then date `reference_time`. The times must be on one count, as
[`calc_leap_second_corrections`](@ref) puts them.
"""
function calc_pseudoranges(times)
    pseudoranges = collect(Float64, times)
    reference_time = calc_pseudoranges!(pseudoranges)
    pseudoranges, reference_time
end

"""
    calc_pseudoranges!(times) -> reference_time

[`calc_pseudoranges`](@ref) in place: `times` is overwritten with the pseudoranges, and only
the reference time is returned. Allocates nothing.
"""
function calc_pseudoranges!(times)
    t_ref = maximum(times)
    # A crossover can fall between the times: within the transmit-time spread, or, in a
    # mixed solve, between BeiDou times on the GPS count (`SOW + 14`, past 604800 for
    # the first 14 s of a GPS week) and GPS times already wrapped.
    for j in eachindex(times)
        times[j] = fold_week_crossover(t_ref - times[j]) * SPEED_OF_LIGHT
    end
    return t_ref
end

"""
    calc_uncorrected_transmit_time(state::SatelliteState) -> Float64

The satellite's *raw* transmit time, a seconds-of-week count on its own system's scale:
the decoder's time of week plus the bit, code-phase and carrier-phase fractions, before
the clock polynomial, relativistic and group-delay corrections that
[`calc_transmit_time`](@ref) applies. It is the time a RINEX pseudorange is measured
against, as a RINEX reader applies the broadcast clock correction itself.
"""
function calc_uncorrected_transmit_time(state::SatelliteState)
    signal = state.signal
    # Time of week at the navigation-frame reference epoch, from GNSSDecoder (which
    # reconstructs it where a message has no seconds-of-week field, as L1C-D and B1C).
    t_tow = get_time_of_week(state.decoder)
    # The bits are counted at the decoder's data rate, not the ranging `signal`'s: the
    # range may be generated on a dataless pilot (L1C-P, E1C) whose data rate is 0 Hz.
    t_bits =
        state.decoder.num_bits_after_valid_syncro_sequence /
        get_data_frequency(state.decoder) * Hz
    # With joint same-band tracking the code phase can run up to a whole data symbol,
    # which `t_bits` already counts; keep only the residual within one symbol.
    chips_per_symbol =
        ustrip(Hz, get_code_frequency(signal)) / ustrip(Hz, get_data_frequency(state.decoder))
    t_code_phase = mod(state.code_phase, chips_per_symbol) / get_code_frequency(signal) * Hz
    # `carrier_phase` is in radians, as `Tracking.get_carrier_phase` reports it.
    t_carrier_phase = state.carrier_phase / 2π / get_center_frequency(signal) * Hz

    t_tow + t_bits + t_code_phase + t_carrier_phase
end

function calc_relativistic_correction(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - decoded(data.t_0e))
    # √A of the effective elements (`√(A_REF + ΔA)` for CNAV/CNAV-2, which carry none).
    el = orbital_elements(data, decoder.constants.μ, time_from_ephemeris_reference_epoch)
    E = calc_eccentric_anomaly(decoder, t)
    decoder.constants.F * decoded(data.e) * el.sqrt_A * sin(E)
end

function correct_clock(decoder::GNSSDecoder.GNSSDecoderState, signal, t)
    Δtr = calc_relativistic_correction(decoder, t)
    data = decoder.data
    # Unfolded, a week crossover puts ~a_f1·604800 s (kilometres) into the polynomial.
    Δt_from_reference = fold_week_crossover(t - decoded(data.t_0c))
    Δt =
        decoded(data.a_f0) +
        decoded(data.a_f1) * Δt_from_reference +
        decoded(data.a_f2) * Δt_from_reference^2 +
        Δtr
    t - correct_by_group_delay(decoder, signal, Δt)
end

"""
    calc_satellite_clock_drift(decoder::GNSSDecoderState, t) -> Float64

Rate of the satellite clock correction (s/s): the derivative of the clock polynomial
about its reference epoch `t_0c`. (About `t` instead, the `a_f2` term would gain up to
~4e-9 s/s, ≈ 1.3 m/s per satellite, which the velocity would absorb.) The rate of the
relativistic term, at most ~1 mm/s, is not modelled.
"""
function calc_satellite_clock_drift(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    decoded(data.a_f1) + 2 * decoded(data.a_f2) * fold_week_crossover(t - decoded(data.t_0c))
end

# Group-delay / inter-signal correction: selected by the ranging `signal` (which may be a
# pilot such as L1C-P), the values from the `decoder`'s message (decoded on the band's data
# component). GPS has an ISC per signal (IS-GPS-705J §20.3.3.3.1.2, IS-GPS-800J §3.5.4.1);
# Galileo's group delay is per band.

# GPS L1 C/A via LNAV: only T_GD (LNAV carries no ISCs).
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData},
    ::GNSSSignals.GPSL1CA,
    t,
) = t - decoder.data.T_GD

# A CNAV/CNAV-2 group-delay term may not be broadcast yet when ephemeris and clock are
# complete, so a missing one is taken as zero: a ~metre bias, never a blocked fix. (T_GD,
# L1 P(Y)−L2 P(Y) per IS-GPS-200 §20.3.3.3.3.2, is the same on every GPS signal and could
# come from another decoder of the PRN, but decoders are independent here.)
group_delay_term(x) = something(x, 0.0)

# GPS CNAV (L5 and L2C share its container): −T_GD + the ranging signal's ISC
# (IS-GPS-705J §20.3.3.3.1.2, IS-GPS-200N §30.3.3.3.1.1).
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::GNSSSignals.GPSL5I,
    t,
) = t - group_delay_term(decoder.data.T_GD) + group_delay_term(decoder.data.ISC_L5I5)
# The L5 pilot has its own ISC_L5Q5; L2 CL shares ISC_L2C with CM (IS-GPS-200N
# §30.3.3.3.1.1 defines one for L2C).
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

# GPS CNAV-2 (L1C) carries the full L1 ISC set: L1C-D, L1C-P or C/A.
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

# Galileo: the broadcast group delay is per band (E1B and E1C share one), and two facts
# set it (Galileo OS SIS ICD, Issue 2.2, §5.1.5):
#
# - Which BGD: the clock is referred to an ionosphere-free pair, E1/E5b for I/NAV (E1-B,
#   E5b-I) and E1/E5a for F/NAV (E5a-I). I/NAV broadcasts both BGDs, but an E5a range
#   carries an F/NAV decoder, so it only ever needs its own.
# - How much: on the pair's second band the BGD scales by `(f_E1/f_band)²`, 1.79 on E5a
#   and 1.70 on E5b; dropping it leaves a per-satellite metre-level error. Derived from
#   the carriers, the factor is exactly 1 on E1.
#
# No `group_delay_term`: both messages gate their BGDs in
# `is_decoding_completed_for_positioning`, so a `nothing` here is a decoder bug that a
# zero would hide.
galileo_group_delay_scaling(signal) =
    (get_center_frequency(GalileoE1B) / get_center_frequency(signal))^2

# One method per (message, ranging signal) that can occur, and no fallback, as for GPS
# and BeiDou (checked by introspection in `test/signals/galileo_e5b_e6b.jl`): I/NAV serves E1
# and E5b ranges, F/NAV E5a ones. An E6 range (E6 has no BGD) or a cross-band pairing
# raises rather than returning a plausible number.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoINAVData},
    signal::Union{
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
    galileo_group_delay_scaling(signal) * decoder.data.BGD_E1_E5b

correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE5aData},
    signal::Union{
        GNSSSignals.GalileoE5aI,
        GNSSSignals.GalileoE5aQ,
        GNSSSignals.GalileoE5aQP,
    },
    t,
) =
    t -
    galileo_group_delay_scaling(signal) * decoder.data.BGD_E1_E5a

# BeiDou refers its clock to B3I on every message and publishes a group delay per signal
# from there (BDS-SIS-ICD-B1I-3.0 §5.2.4.10, BDS-SIS-ICD-B1C-1.0 §7.6 and the B2a/B2b
# counterparts), so B3I needs none. Data components add an ISC to their band's pilot
# delay (`ISC_B1Cd`, `ISC_B2ad`), as GPS's do.
#
# Only B2a's terms use `group_delay_term`, mirroring the decoder's gating:
#
#   - D1/D2 gates `T_GD1` (not `T_GD2`; see the B1I method), B1C gates `T_GD_B1Cp` and
#     `ISC_B1Cd` (bits 546-569 of subframe 2, with the ephemeris), B2b gates `T_GD_B2bI`
#     (same MT30 as the clock). They are read unwrapped: a `nothing` is a decoder bug.
#   - B2a's terms ride in MT30 only, while the clock is in all of MT30-34 and
#     BDS-SIS-ICD-B2a-1.0 §6.2 schedules only MT10/11, so gating on them would wait on an
#     unbounded interval; they can be absent and are taken as zero.

# Legacy D1/D2 on B1I: T_GD1. Its T_GD2, for BDS-2's B2I, has no consumer: B2I has no
# signal definition in GNSSSignals and will not get one (JuliaGNSS/GNSSSignals.jl#156).
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

# B-CNAV1 (B1C): the B1C pilot group delay plus the data ISC. (Its B2a pilot delay is
# unused: a B2a range carries a B2a decoder.)
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
# B-CNAV2 (B2a) mirrors it.
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

# No catch-all: the methods cover exactly the ranging signals of the band each message is
# decoded on, and a cross-band pairing cannot arise. A `MethodError` says so where
# returning `t` would pass for a correct zero.

"""
    calc_transmit_time(state::SatelliteState) -> Float64

The satellite's transmit time as a seconds-of-week count on its own system's scale: the
decoder's time of week plus the bit, code-phase and carrier-phase fractions, corrected by
the clock polynomial, the relativistic term and the ranging signal's group delay. Every
ephemeris and clock evaluation takes it unconverted (see
[`leap_second_offset_to_gpst`](@ref)).
"""
function calc_transmit_time(state::SatelliteState)
    approximated_time = calc_uncorrected_transmit_time(state)
    correct_clock(state.decoder, state.signal, approximated_time)
end
