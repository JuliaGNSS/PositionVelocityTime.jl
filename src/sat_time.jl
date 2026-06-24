function correct_week_crossovers(t)
    half_week = 302400  #Half of Week in Seconds
    t + (t > half_week ? -2 * half_week : (t < -half_week ? 2 * half_week : 0.0))
end

function calc_uncorrected_time(state::SatelliteState)
    system = state.system
    t_tow = state.decoder.data.TOW
    t_bits =
        state.decoder.num_bits_after_valid_syncro_sequence /
        GNSSSignals.get_data_frequency(system) * Hz
    t_code_phase = state.code_phase / GNSSSignals.get_code_frequency(system) * Hz
    t_carrier_phase = state.carrier_phase / GNSSSignals.get_center_frequency(system) * Hz

    t_tow + t_bits + t_code_phase + t_carrier_phase
end

function calc_relativistic_correction(decoder::GNSSDecoder.GNSSDecoderState, t)
    E = calc_eccentric_anomaly(decoder, t)
    decoder.constants.F * decoder.data.e * decoder.data.sqrt_A * sin(E)
end

function correct_clock(decoder::GNSSDecoder.GNSSDecoderState, t)
    Δtr = calc_relativistic_correction(decoder, t)
    Δt =
        decoder.data.a_f0 +
        decoder.data.a_f1 * (t - decoder.data.t_0c) +
        decoder.data.a_f2 * (t - decoder.data.t_0c)^2 +
        Δtr
    t - correct_by_group_delay(decoder, Δt)
end

function calc_satellite_clock_drift(decoder::GNSSDecoder.GNSSDecoderState, t)
    decoder.data.a_f1 +
    decoder.data.a_f2 * t * 2
end

function correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData},
    t,
)
    t - decoder.data.T_GD
end

function correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
    t,
)
    t - decoder.data.broadcast_group_delay_e1_e5b
end

function calc_corrected_time(state::SatelliteState)
    approximated_time = calc_uncorrected_time(state)
    correct_clock(state.decoder, approximated_time)
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
function ggto_available(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
)
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
function calc_ggto_offset(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
    t,
)
    data = decoder.data
    Δweek = mod(data.WN - data.WN_0G, 64)
    data.A_0G + data.A_1G * (t - data.t_0G + 604800 * Δweek)
end
