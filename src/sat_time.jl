function correct_week_crossovers(t)
    half_week = 302400  #Half of Week in Seconds
    t + (t > half_week ? -2 * half_week : (t < -half_week ? 2 * half_week : 0.0))
end

# Time of week (seconds) at the navigation-frame reference epoch. GPS LNAV, GPS
# CNAV (L5) and Galileo broadcast it directly as TOW; GPS CNAV-2 (L1C) instead
# broadcasts the two-hour interval count ITOW and the 18 s time-of-interval TOI
# (IS-GPS-800J §3.5.3), so reconstruct TOW = ITOW·7200 + TOI·18.
get_tow(decoder::GNSSDecoder.GNSSDecoderState) = decoder.data.TOW
get_tow(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData}) =
    decoder.data.ITOW * 7200 + decoder.data.toi * 18

function calc_uncorrected_time(state::SatelliteState)
    system = state.system
    t_tow = get_tow(state.decoder)
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
    t_code_phase = state.code_phase / get_code_frequency(system) * Hz
    t_carrier_phase = state.carrier_phase / get_center_frequency(system) * Hz

    t_tow + t_bits + t_code_phase + t_carrier_phase
end

function calc_relativistic_correction(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - data.t_0e)
    # √A from the effective elements: the broadcast `sqrt_A` directly for LNAV/Galileo,
    # `√(A_REF + ΔA)` for CNAV/CNAV-2 (which carry no `sqrt_A` field).
    el = orbital_elements(data, decoder.constants.μ, time_from_ephemeris_reference_epoch)
    E = calc_eccentric_anomaly(decoder, t)
    decoder.constants.F * data.e * el.sqrt_A * sin(E)
end

function correct_clock(decoder::GNSSDecoder.GNSSDecoderState, system, t)
    Δtr = calc_relativistic_correction(decoder, t)
    Δt =
        decoder.data.a_f0 +
        decoder.data.a_f1 * (t - decoder.data.t_0c) +
        decoder.data.a_f2 * (t - decoder.data.t_0c)^2 +
        Δtr
    t - correct_by_group_delay(decoder, system, Δt)
end

function calc_satellite_clock_drift(decoder::GNSSDecoder.GNSSDecoderState, t)
    decoder.data.a_f1 +
    decoder.data.a_f2 * t * 2
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

# Galileo: the broadcast group delay is per band (E1 vs E5a), so it depends only on
# the decoder's message, not on whether the range came from the data or the pilot
# component — E1B/E1C share BGD(E1,E5b); E5a-I/E5a-Q share BGD(E1,E5a).
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
    ::AbstractGNSSSignal,
    t,
) = t - decoder.data.broadcast_group_delay_e1_e5b
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE5aData},
    ::AbstractGNSSSignal,
    t,
) = t - decoder.data.broadcast_group_delay_e1_e5a

function calc_corrected_time(state::SatelliteState)
    approximated_time = calc_uncorrected_time(state)
    correct_clock(state.decoder, state.system, approximated_time)
end

"""
    ggto_available(decoder) -> Bool

Return `true` if `decoder` carries a complete Galileo–GPS Time Offset (GGTO)
record (Galileo word type 10: `A_0G`, `A_1G`, `t_0G`, `WN_0G`). The GGTO lets
the receiver express Galileo System Time in GPS time, which makes it possible to
combine GPS and Galileo satellites when the geometry is too weak to estimate an
independent Galileo clock bias. Always `false` for non-Galileo systems.
"""
ggto_available(::GNSSDecoder.GNSSDecoderState) = false
function ggto_available(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGalileoData})
    data = decoder.data
    !isnothing(data.A_0G) &&
        !isnothing(data.A_1G) &&
        !isnothing(data.t_0G) &&
        !isnothing(data.WN_0G)
end

"""
    calc_ggto_offset(decoder, t) -> Float64

Galileo–GPS Time Offset `Δt_systems = GST − GPST` in seconds at Galileo time of
week `t`, per the Galileo OS SIS ICD (word type 10):

    Δt_systems = A_0G + A_1G · (t − t_0G + 604800 · ((WN − WN_0G) mod 64))

`WN_0G` is the 6-bit GGTO reference week, so the week difference is taken modulo
64. To convert a Galileo system time to GPS time, subtract this offset.
"""
function calc_ggto_offset(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGalileoData}, t)
    data = decoder.data
    Δweek = mod(data.WN - data.WN_0G, 64)
    data.A_0G + data.A_1G * (t - data.t_0G + 604800 * Δweek)
end
