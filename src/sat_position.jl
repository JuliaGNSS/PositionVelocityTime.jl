"""
    Calculates ECI Position of SV

    $SIGNATURES
    ´decoder_state´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the SV in ECI coorinates.
    The implementation follows IS-GPS-200K
"""
function sat_position_ECI(decoder_state::GNSSDecoderState, code_phase, carrier_phase = 0)
    is_sat_healthy_and_decodable(decoder_state) || throw(BadData("SV not decoded properly"))
    
    t = calc_corrected_time(decoder_state, code_phase, carrier_phase)
    tk = t - decoder_state.data.t_oe
    tk = check_crossover(tk)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(decoder_state, code_phase, carrier_phase)
    dt_sv = code_phase_offset(decoder_state, t_sv)
    Ek = calc_eccentric_anomaly(decoder_state, t, dt_sv)
    
    
    
    
    e  = decoder_state.data.e
    vk =  2 * atan(sqrt((1+e) / (1-e)) * tan(Ek/2))
    # Argument of latitude
    Phik = vk + decoder_state.data.ω
    
    # Second harmonic perturbations
    duk = decoder_state.data.C_us * sin(2*Phik) + decoder_state.data.C_uc * cos(2*Phik);
    drk = decoder_state.data.C_rs * sin(2*Phik) + decoder_state.data.C_rc * cos(2*Phik);
    dik = decoder_state.data.C_is * sin(2*Phik) + decoder_state.data.C_ic * cos(2*Phik);
    
    
    # Corrected orbit parameters
    uk = Phik + duk;
    rk = decoder_state.data.sqrt_A^2 * (1 - e * cos(Ek)) + drk;
    ik = decoder_state.data.i_0 + dik + decoder_state.data.IDOT * tk;
    
    # Positions in orbital plane
    xks = rk * cos(uk);
    yks = rk * sin(uk);
    
    # Corrected longitude of ascending node
    Omegak = decoder_state.data.Ω_0 + decoder_state.data.Ω_dot * tk;
    
    # ECI coordinates
    swk = sin(Omegak);
    cwk = cos(Omegak);
    cik = cos(ik);
    sik = sin(ik);
            
    xk = xks * cwk - yks * cik * swk;
    yk = xks * swk + yks * cik * cwk;
    zk = yks * sik;
    
    return [xk, yk, zk]
end


"""
    Calculates ECEF Position by converting ECI position

    $SIGNATURES
    ´decoder_state´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure

    This function calculates the position of the SV in ECI coorinates and
    converts this position into ECEF coordinates.
    The implementation follows IS-GPS-200K
"""
function sat_position_ECI_2_ECEF(decoder_state::GNSSDecoderState, code_phase, carrier_phase = 0)
    pos_ECI = sat_position_ECI( decoder_state, code_phase, carrier_phase)
            
    tk = calc_corrected_time(decoder_state, code_phase, carrier_phase)
    tk = tk - decoder_state.data.t_oe
    tk = check_crossover( tk)
    theta = decoder_state.constants.Ω_dot_e * (tk + decoder_state.data.t_oe)
    cost  = cos(theta)
    sint  = sin(theta)
    R = [cost sint; -sint cost] 
    pos_ECI[1:2] = R * pos_ECI[1:2]
    return pos_ECI
end


"""
    Calculates ECEF Position of SV 

    $SIGNATURES
    ´decoder_state´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the SV in ECEF coorinates.
    The implementation follows IS-GPS-200K Table 20-IV.
"""
function sat_position_ECEF(decoder_state::GNSSDecoderState, code_phase, carrier_phase)
    
    is_sat_healthy_and_decodable(decoder_state) || throw(BadData("SV not decoded properly"))
    F = decoder_state.constants.F
    e = decoder_state.data.e
    μ = decoder_state.constants.μ    
    Ω_dot_e = decoder_state.constants.Ω_dot_e

    A = decoder_state.data.sqrt_A^2
    n_0 = sqrt(μ / (A^3))
    n = n_0 + decoder_state.data.Δn
    
    
    t = calc_corrected_time(decoder_state, code_phase, carrier_phase)
    tk = t - decoder_state.data.t_oe
    tk = check_crossover(tk)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(decoder_state, code_phase, carrier_phase)
    dt_sv = code_phase_offset(decoder_state, t_sv)
    E = calc_eccentric_anomaly(decoder_state, t, dt_sv)
    
    vk = 2 * atan(sqrt((1+e) / (1-e)) * tan(E/2))
    
    Φk = vk + decoder_state.data.ω
    
    duk = decoder_state.data.C_us * sin(2*Φk) + decoder_state.data.C_uc * cos(2*Φk);
    drk = decoder_state.data.C_rs * sin(2*Φk) + decoder_state.data.C_rc * cos(2*Φk);
    dik = decoder_state.data.C_is * sin(2*Φk) + decoder_state.data.C_ic * cos(2*Φk);
    
    uk = Φk + duk
    rk = A * (1 - e*cos(E)) + drk
    ik = decoder_state.data.i_0 + dik + decoder_state.data.IDOT * tk
    
    xks = rk * cos(uk)
    yks = rk * sin(uk)
    
    Ωk = decoder_state.data.Ω_0 + (decoder_state.data.Ω_dot - Ω_dot_e) * tk - Ω_dot_e*decoder_state.data.t_oe
    
    xk = xks * cos(Ωk) - yks * cos(ik)*sin(Ωk)
    yk = xks * sin(Ωk) + yks * cos(ik)*cos(Ωk)
    zk = yks * sin(ik)
    
    position = [xk, yk, zk]
    return position
end

"""
    Checks if satellite data is correctly tranmitted
    $SIGNATURES
    ´decoder_state´: Decoder

    Checks if satellite data contains the needed information and 
    if errors during decoding occured
"""
function is_sat_healthy_and_decodable(decoder_state::GNSSDecoderState)
    decoder_state.data.svhealth == "000000" ? decoder_state.subframes_decoded[1:3] == [1,1,1] : false
end