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

function calc_satellite_clock_drift(state::SatelliteState)
    approximated_time = calc_uncorrected_time(state)
    calc_satellite_clock_drift(state.decoder, approximated_time)
end

function correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data},
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