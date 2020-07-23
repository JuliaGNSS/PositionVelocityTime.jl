function satellite_position_ECI(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    
    
    t = calc_corrected_time(dc, code_phase, carrier_phase)
    tk = t - dc.data.t_oe
    tk = check_crossover(dc, tk)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(dc, code_phase, carrier_phase)
    dt_sv = code_phase_offset(dc, t_sv)
    Ek = calc_eccentric_anomaly( dc, t, dt_sv )
    
    
    
    
    e  = dc.data.e
    vk = atan( sqrt( 1 - e^2 ) * sin(Ek), cos(Ek) - e )
    # Argument of latitude
    Phik = vk + dc.data.ω
    
    # Second harmonic perturbations
    duk = dc.data.C_us .* sin( 2*Phik ) + dc.data.C_uc .* cos( 2*Phik );
    drk = dc.data.C_rs .* sin( 2*Phik ) + dc.data.C_rc .* cos( 2*Phik );
    dik = dc.data.C_is .* sin( 2*Phik ) + dc.data.C_ic .* cos( 2*Phik );
    
    
    # Corrected orbit parameters
    uk = Phik + duk;
    rk = dc.data.sqrt_A .^2 .* (1 - e .* cos(Ek)) + drk;
    ik = dc.data.i_0 + dik + dc.data.IDOT * tk;
    
    # Positions in orbital plane
    xks = rk .* cos(uk);
    yks = rk .* sin(uk);
    
    # Corrected longitude of ascending node
    Omegak = dc.data.Ω_0 + dc.data.Ω_dot * tk;
    
    # ECI coordinates
    swk = sin( Omegak );
    cwk = cos( Omegak );
    cik = cos( ik );
    sik = sin( ik );
            
    xk = xks .* cwk - yks .* cik .* swk;
    yk = xks .* swk + yks .* cik .* cwk;
    zk = yks .* sik;
            
            
    return [xk, yk, zk]
end


function satellite_position_ECI_2_ECEF( dc::GNSSDecoderState, code_phase, carrier_phase = 0 )
    # This function implements the calculation of the satellite
    # position in ECEF coordinates.
    # This implementation follows IS-GPS-200K
    pos_ECI = satellite_position_ECI( dc, code_phase, carrier_phase)
            
    tk = calc_corrected_time(dc, code_phase, carrier_phase)
    tk = tk - dc.data.t_oe
    tk = check_crossover(dc, tk)
        
    theta = dc.constants.Ω_dot_e * ( tk + dc.data.t_oe )
    cost  = cos(theta)
    sint  = sin(theta)
    R = [cost sint; -sint cost] 
    pos_ECI[1:2] = R * pos_ECI[1:2]
    return pos_ECI
end

function sat_position_ECEF(dc::GNSSDecoderState, code_phase, carrier_phase = 0)

    tk = calc_corrected_time(dc, code_phase)
    Ek = calc_eccentric_anomaly(dc, tk)
    e = dc.data.e
    A = dc.data.sqrt_A^2
    vk = atan( sqrt( 1 - e^2 ) * sin(Ek), cos(Ek) - e )
    phi = vk + dc.data.ω
    
    u = phi +
        dc.data.C_uc * cos(2*phi) + 
        dc.data.C_us * sin(2*phi)
    
    
    r = A * (1 - dc.data.e*cos(Ek)) + 
        dc.data.C_rc * cos(2*phi) + 
        dc.data.C_rs * sin(2*phi)
    
    i = dc.data.i_0 + dc.data.IDOT * tk + 
        dc.data.C_ic * cos(2*phi) + 
        dc.data.C_is * sin(2*phi)
    
    Omega = dc.data.Ω_0 + (dc.data.Ω_dot - dc.constants.Ω_dot_e) * tk - 
            dc.data.Ω_dot * dc.data.t_oe
    
    satPositions_x = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega)
    satPositions_y = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega)
    satPositions_z = sin(u)*r * sin(i)
    
    
    
    return [satPositions_x, satPositions_y, satPositions_z]
end
