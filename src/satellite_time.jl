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
    # The bits are counted at the decoder's bit rate, not the ranging `signal`'s:
    # the range may be generated on a dataless pilot (L1C-P, E1C) whose rate is 0 Hz.
    # GNSSDecoder's `num_symbols_since_sync` counts these tracked bits.
    t_bits =
        state.decoder.num_symbols_since_sync / get_bit_rate(state.decoder) * Hz
    # With joint same-band tracking the code phase can run up to a whole bit,
    # which `t_bits` already counts; keep only the residual within one bit.
    chips_per_bit =
        ustrip(Hz, get_chip_rate(signal)) / ustrip(Hz, get_bit_rate(state.decoder))
    t_code_phase = mod(state.code_phase, chips_per_bit) / get_chip_rate(signal) * Hz
    # `carrier_phase` is in radians, as `Tracking.get_carrier_phase` reports it.
    t_carrier_phase = state.carrier_phase / 2π / get_center_frequency(signal) * Hz

    t_tow + t_bits + t_code_phase + t_carrier_phase
end

# The relativistic clock term `F · e · √A · sin E` at `t`.
function calc_relativistic_correction(ephemeris::BroadcastEphemeris, t)
    E = calc_eccentric_anomaly(ephemeris, t)
    ephemeris.F * ephemeris.e * sqrt(ephemeris.A) * sin(E)
end

# The group delay of `signal`'s range, taken as zero until GNSSDecoder has every term of
# it: GPS CNAV/CNAV-2 and BeiDou B2a broadcast theirs apart from the clock, on their own
# cadence, and a single-frequency fix without one is off by about a metre rather than
# waiting on an interval the ICDs do not bound.
group_delay(decoder, signal) = something(get_group_delay(decoder, signal), 0.0)

# The time `t` corrected by the satellite clock: the polynomial, the relativistic term
# and the ranging signal's `group_delay`.
function correct_clock(ephemeris::BroadcastEphemeris, group_delay, t)
    # Unfolded, a week crossover puts ~a_f1·604800 s (kilometres) into the polynomial.
    Δt_from_reference = fold_week_crossover(t - ephemeris.t_0c)
    Δt =
        ephemeris.a_f0 + ephemeris.a_f1 * Δt_from_reference +
        ephemeris.a_f2 * Δt_from_reference^2 + calc_relativistic_correction(ephemeris, t)
    t - (Δt - group_delay)
end

"""
    calc_satellite_clock_drift(ephemeris::BroadcastEphemeris, t) -> Float64

Rate of the satellite clock correction (s/s) of GNSSDecoder's `BroadcastEphemeris`: the
derivative of the clock polynomial about its reference epoch `t_0c`. (About `t` instead,
the `a_f2` term would gain up to ~4e-9 s/s, ≈ 1.3 m/s per satellite, which the velocity
would absorb.) The rate of the relativistic term, at most ~1 mm/s, is not modelled.
"""
calc_satellite_clock_drift(ephemeris::BroadcastEphemeris, t) =
    ephemeris.a_f1 + 2 * ephemeris.a_f2 * fold_week_crossover(t - ephemeris.t_0c)

"""
    calc_transmit_time(state::SatelliteState) -> Float64

The satellite's transmit time as a seconds-of-week count on its own system's scale: the
decoder's time of week plus the bit, code-phase and carrier-phase fractions, corrected
by the clock polynomial, the relativistic term and the ranging signal's group delay
(GNSSDecoder's `get_broadcast_ephemeris` and `get_group_delay`; a group delay not yet
decoded in full counts as zero, a metre or so on one signal). Every ephemeris and clock
evaluation takes it unconverted (see `calc_leap_second_corrections!`). The decoder must be
complete for positioning.
"""
calc_transmit_time(state::SatelliteState) = correct_clock(
    something(get_broadcast_ephemeris(state.decoder)),
    group_delay(state.decoder, state.signal),
    calc_uncorrected_transmit_time(state),
)
