"""
    Calculates ECI Position of SV

    $SIGNATURES
    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the SV in ECI coorinates.
    The implementation follows IS-GPS-200K
"""
function sat_position_ECI(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    can_get_sat_position(dc) || throw(BadData("SV not decoded properly"))
    
    t = calc_corrected_time(dc, code_phase, carrier_phase)
    tk = t - dc.data.t_oe
    tk = check_crossover(tk)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(dc, code_phase, carrier_phase)
    dt_sv = calc_time_corr(dc, t_sv)
    Ek = calc_eccentric_anomaly( dc, t, dt_sv )
    
    
    
    
    e  = dc.data.e
    vk =  2 * atan(sqrt( (1+e) / (1-e) ) * tan(Ek / 2) )
    # Argument of latitude
    Phik = vk + dc.data.ω
    
    # Second harmonic perturbations
    duk = dc.data.C_us * sin( 2*Phik ) + dc.data.C_uc * cos( 2*Phik );
    drk = dc.data.C_rs * sin( 2*Phik ) + dc.data.C_rc * cos( 2*Phik );
    dik = dc.data.C_is * sin( 2*Phik ) + dc.data.C_ic * cos( 2*Phik );
    
    
    # Corrected orbit parameters
    uk = Phik + duk;
    rk = dc.data.sqrt_A^2 * (1 - e * cos(Ek)) + drk;
    ik = dc.data.i_0 + dik + dc.data.IDOT * tk;
    
    # Positions in orbital plane
    xks = rk * cos(uk);
    yks = rk * sin(uk);
    
    # Corrected longitude of ascending node
    Omegak = dc.data.Ω_0 + dc.data.Ω_dot * tk;
    
    # ECI coordinates
    swk = sin( Omegak );
    cwk = cos( Omegak );
    cik = cos( ik );
    sik = sin( ik );
            
    xk = xks * cwk - yks * cik * swk;
    yk = xks * swk + yks * cik * cwk;
    zk = yks * sik;
    
    return Vector([xk, yk, zk])
end


"""
    Calculates ECEF Position by converting ECI position

    $SIGNATURES
    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure

    This function calculates the position of the SV in ECI coorinates and
    converts this position into ECEF coordinates.
    The implementation follows IS-GPS-200K
"""
function sat_position_ECI_2_ECEF( dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    pos_ECI = satellite_position_ECI( dc, code_phase, carrier_phase)
            
    tk = calc_corrected_time(dc, code_phase, carrier_phase)
    tk = tk - dc.data.t_oe
    tk = check_crossover( tk)
    theta = dc.constants.Ω_dot_e * ( tk + dc.data.t_oe )
    cost  = cos(theta)
    sint  = sin(theta)
    R = [cost sint; -sint cost] 
    pos_ECI[1:2] = R * pos_ECI[1:2]
    return pos_ECI
end


"""
    Calculates ECEF Position of SV 

    $SIGNATURES
    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the SV in ECEF coorinates.
    The implementation follows IS-GPS-200K Table 20-IV.
"""
function sat_position_ECEF(dc::GNSSDecoderState, code_phase, carrier_phase)
    
    can_get_sat_position(dc) || throw(BadData("SV not decoded properly"))
    F = dc.constants.F
    e = dc.data.e
    μ = dc.constants.μ    
    Ω_dot_e = dc.constants.Ω_dot_e

    A = dc.data.sqrt_A ^ 2
    n_0 = sqrt(μ / (A^3))
    n = n_0 + dc.data.Δn
    
    
    t = calc_corrected_time(dc, code_phase, carrier_phase)
    tk = t - dc.data.t_oe
    tk = check_crossover(tk)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(dc, code_phase, carrier_phase)
    dt_sv = calc_time_corr(dc, t_sv)
    E = calc_eccentric_anomaly( dc, t, dt_sv )
    
    vk = 2 * atan(sqrt( (1+e) / (1-e) ) * tan(E / 2) )
    
    Φk = vk + dc.data.ω
    
    duk = dc.data.C_us * sin( 2*Φk ) + dc.data.C_uc * cos( 2*Φk );
    drk = dc.data.C_rs * sin( 2*Φk ) + dc.data.C_rc * cos( 2*Φk );
    dik = dc.data.C_is * sin( 2*Φk ) + dc.data.C_ic * cos( 2*Φk );
    
    uk = Φk + duk
    rk = A * (1 - e*cos(E)) + drk
    ik = dc.data.i_0 + dik + dc.data.IDOT * tk
    
    xks = rk * cos(uk)
    yks = rk * sin(uk)
    
    Ωk = dc.data.Ω_0 + (dc.data.Ω_dot - Ω_dot_e) * tk - Ω_dot_e*dc.data.t_oe
    
    xk = xks * cos(Ωk) - yks * cos(ik)*sin(Ωk)
    yk = xks * sin(Ωk) + yks * cos(ik)*cos(Ωk)
    zk = yks * sin(ik)
    
    position = Vector([xk, yk, zk])
    return position, tk
end

"""
    Checks if satellite data is correctly tranmitted
    $SIGNATURES
    ´dc´: Decoder

    Checks if satellite data contains the needed information and 
    if errors during decoding occured
"""
function can_get_sat_position(dc::GNSSDecoderState)
    status = true
    
    dc.data.svhealth == "000000" ? status : status = false
    dc.subframes_decoded[1:3] == [1,1,1] ? status : status = false
    return status
end