"""
    Calculates ECI Position of SV

    $SIGNATURES
    ´decoder_state´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the SV in ECI coorinates.
    The implementation follows IS-GPS-200K
"""
function sat_position_ECI(sat_state::SatelliteState)
    is_sat_healthy_and_decodable(sat_state.decoder_state) || throw(BadData("SV not decoded properly"))
    
    t = calc_corrected_time(sat_state)
    tₖ = t - sat_state.decoder_state.data.t_oe
    tₖ = check_crossover(tₖ)
    
    # Excentric and true anomaly
    tₛᵥ = calc_uncorrected_time(sat_state)
    dtₛᵥ = code_phase_offset(sat_state.decoder_state, tₛᵥ)
    Eₖ = calc_eccentric_anomaly(sat_state.decoder_state, t, dtₛᵥ)
    
    
    e  = sat_state.decoder_state.data.e
    vₖ =  2 * atan(sqrt((1 + e) / (1 - e)) * tan(Eₖ / 2))
    # Argument of latitude
    Φₖ = vₖ + sat_state.decoder_state.data.ω
    
    # Second harmonic perturbations
    δuₖ = sat_state.decoder_state.data.C_us * sin(2 * Φₖ) + sat_state.decoder_state.data.C_uc * cos(2 * Φₖ);
    δrₖ = sat_state.decoder_state.data.C_rs * sin(2 * Φₖ) + sat_state.decoder_state.data.C_rc * cos(2 * Φₖ);
    δiₖ = sat_state.decoder_state.data.C_is * sin(2 * Φₖ) + sat_state.decoder_state.data.C_ic * cos(2 * Φₖ);
    
    
    # Corrected orbit parameters
    uₖ = Φₖ + δuₖ;
    rₖ = sat_state.decoder_state.data.sqrt_A^2 * (1 - e * cos(Eₖ)) + δrₖ;
    iₖ = sat_state.decoder_state.data.i_0 + δiₖ + sat_state.decoder_state.data.IDOT * tₖ;
    
    # Positions in orbital plane
    xₖs = rₖ * cos(uₖ);
    yₖs = rₖ * sin(uₖ);
    
    # Corrected longitude of ascending node
    Ωₖ = sat_state.decoder_state.data.Ω_0 + sat_state.decoder_state.data.Ω_dot * tₖ;
    
    # ECI coordinates
    sinΩₖ = sin(Ωₖ);
    cosΩₖ = cos(Ωₖ);
    cosiₖ = cos(iₖ);
    siniₖ = sin(iₖ);
            
    xₖ = xₖs * cosΩₖ - yₖs * cosiₖ * sinΩₖ;
    yₖ = xₖs * sinΩₖ + yₖs * cosiₖ * cosΩₖ;
    zₖ = yₖs * siniₖ;
    
    return [xₖ, yₖ, zₖ]
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
function sat_position_ECI_2_ECEF(sat_state::SatelliteState)
    pos_ECI = sat_position_ECI(sat_state)
            
    tₖ = calc_corrected_time(sat_state)
    tₖ = tₖ - sat_state.decoder_state.data.t_oe
    tₖ = check_crossover(tₖ)
    θ = sat_state.decoder_state.constants.Ω_dot_e * (tₖ + sat_state.decoder_state.data.t_oe)
    cosθ  = cos(θ)
    sinθ  = sin(θ)
    R = [cosθ sinθ; -sinθ cosθ] 
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
function sat_position_ECEF(sat_state)
    
    is_sat_healthy_and_decodable(sat_state.decoder_state) || throw(BadData("SV not decoded properly"))
    F = sat_state.decoder_state.constants.F
    e = sat_state.decoder_state.data.e
    μ = sat_state.decoder_state.constants.μ    
    Ω_dot_e = sat_state.decoder_state.constants.Ω_dot_e

    A = sat_state.decoder_state.data.sqrt_A^2
    n₀ = sqrt(μ / (A^3))
    n = n₀ + sat_state.decoder_state.data.Δn
    
    
    t = calc_corrected_time(sat_state)
    tₖ = t - sat_state.decoder_state.data.t_oe
    tₖ = check_crossover(tₖ)
    
    # Excentric and true anomaly
    t_sv = calc_uncorrected_time(sat_state)
    dt_sv = code_phase_offset(sat_state.decoder_state, t_sv)
    E = calc_eccentric_anomaly(sat_state.decoder_state, t, dt_sv)
    
    vₖ = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    
    Φₖ = vₖ + sat_state.decoder_state.data.ω
    
    δuₖ = sat_state.decoder_state.data.C_us * sin(2 * Φₖ) + sat_state.decoder_state.data.C_uc * cos(2 * Φₖ)
    δrₖ = sat_state.decoder_state.data.C_rs * sin(2 * Φₖ) + sat_state.decoder_state.data.C_rc * cos(2 * Φₖ)
    δiₖ = sat_state.decoder_state.data.C_is * sin(2 * Φₖ) + sat_state.decoder_state.data.C_ic * cos(2 * Φₖ)
    
    uₖ = Φₖ + δuₖ 
    rₖ = A * (1 - e * cos(E)) + δrₖ
    iₖ = sat_state.decoder_state.data.i_0 + δiₖ + sat_state.decoder_state.data.IDOT * tₖ

    xₖs = rₖ * cos(uₖ)
    yₖs = rₖ * sin(uₖ)
    
    Ωₖ = sat_state.decoder_state.data.Ω_0 + 
        (sat_state.decoder_state.data.Ω_dot - Ω_dot_e) * tₖ -
        sat_state.decoder_state.data.t_oe * Ω_dot_e
    
    xₖ = xₖs * cos(Ωₖ) - yₖs * cos(iₖ) * sin(Ωₖ)
    yₖ = xₖs * sin(Ωₖ) + yₖs * cos(iₖ) * cos(Ωₖ)
    zₖ = yₖs * sin(iₖ)
    
    position = [xₖ, yₖ, zₖ]
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