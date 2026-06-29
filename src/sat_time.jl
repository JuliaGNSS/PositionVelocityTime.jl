function correct_week_crossovers(t)
    half_week = 302400  #Half of Week in Seconds
    t + (t > half_week ? -2 * half_week : (t < -half_week ? 2 * half_week : 0.0))
end

# Time of week (seconds) at the navigation-frame reference epoch. GPS LNAV, GPS
# CNAV (L5) and Galileo broadcast it directly as TOW; GPS CNAV-2 (L1C) instead
# broadcasts the two-hour interval count ITOW and the 18 s time-of-interval TOI
# (IS-GPS-800J §3.5.3), so reconstruct TOW = ITOW·7200 + TOI·18.
# NOTE: the L1C-D reconstruction is per-ICD but not yet validated against a real
# CNAV-2 capture (no L1C fixture available); the frame/bit alignment may need
# refinement against real data.
get_tow(decoder::GNSSDecoder.GNSSDecoderState) = decoder.data.TOW
get_tow(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData}) =
    decoder.data.ITOW * 7200 + decoder.data.toi * 18

# Symbol rate of the navigation message the decoder demodulated — the rate that
# converts `num_bits_after_valid_syncro_sequence` to seconds. Taken from the
# *decoder's* data signal, not the tracked `system`: with joint same-band tracking a
# band's pseudorange may be generated on a dataless pilot (e.g. GPS L1C-P,
# Galileo E1C) whose own `get_data_frequency` is 0 Hz, while the navigation bits come
# from the data component the decoder ran on. For a data signal this equals
# `get_data_frequency(system)`, so existing single-signal results are unchanged.
nav_data_frequency(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData}) =
    GNSSSignals.get_data_frequency(GNSSSignals.GPSL1CA())
# GPS CNAV shares one data container across L5-I (100 Hz symbols) and L2C-M (50 Hz),
# so the symbol rate is selected by the decoder's constants type, not the data type.
nav_data_frequency(
    ::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData,<:GNSSDecoder.GPSL5IConstants},
) = GNSSSignals.get_data_frequency(GNSSSignals.GPSL5I())
nav_data_frequency(
    ::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData,<:GNSSDecoder.GPSL2CMConstants},
) = GNSSSignals.get_data_frequency(GNSSSignals.GPSL2CM())
nav_data_frequency(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData}) =
    GNSSSignals.get_data_frequency(GNSSSignals.GPSL1C_D())
nav_data_frequency(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData}) =
    GNSSSignals.get_data_frequency(GNSSSignals.GalileoE1B())
nav_data_frequency(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE5aData}) =
    GNSSSignals.get_data_frequency(GNSSSignals.GalileoE5aI())

function calc_uncorrected_time(state::SatelliteState)
    system = state.system
    t_tow = get_tow(state.decoder)
    # Bit-count term uses the decoder's nav-message symbol rate (the rate the bits
    # were counted at); the code- and carrier-phase terms use the tracked ranging
    # `system`, which may be a pilot distinct from the decoder's data signal.
    t_bits =
        state.decoder.num_bits_after_valid_syncro_sequence /
        nav_data_frequency(state.decoder) * Hz
    t_code_phase = state.code_phase / GNSSSignals.get_code_frequency(system) * Hz
    t_carrier_phase = state.carrier_phase / GNSSSignals.get_center_frequency(system) * Hz

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

# GPS CNAV rides on one shared data container for L5 and L2C; the ISC is picked by
# the ranging signal — L5 I5 for an L5 range, L2C for an L2C-M range
# (IS-GPS-705J §20.3.3.3.1.2 / IS-GPS-200N §30.3.3.3.1.1). -T_GD + ISC.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::GNSSSignals.GPSL5I,
    t,
) = t - decoder.data.T_GD + decoder.data.ISC_L5I5
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSCNAVData},
    ::GNSSSignals.GPSL2CM,
    t,
) = t - decoder.data.T_GD + decoder.data.ISC_L2C

# GPS CNAV-2 (L1C) carries the full L1 ISC set, so one decoder serves a range
# generated on the L1C data (L1C-D), the L1C pilot (L1C-P) or C/A — pick the ISC of
# the signal the range was actually generated on.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1C_D,
    t,
) = t - decoder.data.T_GD + decoder.data.ISC_L1CD
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1C_P,
    t,
) = t - decoder.data.T_GD + decoder.data.ISC_L1CP
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1C_DData},
    ::GNSSSignals.GPSL1CA,
    t,
) = t - decoder.data.T_GD + decoder.data.ISC_L1CA

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
function ggto_available(decoder::GNSSDecoder.GNSSDecoderState{<:GalileoNavData})
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
function calc_ggto_offset(decoder::GNSSDecoder.GNSSDecoderState{<:GalileoNavData}, t)
    data = decoder.data
    Δweek = mod(data.WN - data.WN_0G, 64)
    data.A_0G + data.A_1G * (t - data.t_0G + 604800 * Δweek)
end
