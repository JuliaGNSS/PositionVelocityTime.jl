"""
    Calculates the corrected satellite time
    $SIGNATURES

    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    Provides the implementation of the clock computations
    from raw to corrected
    
    Following instructions in IS-GPS-200K: 
    20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    Equation (1)
"""
function calc_corrected_time(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    
    t_sv = calc_uncorrected_time(dc, code_phase, carrier_phase)
    Δt_sv = code_phase_offset(dc, t_sv)
    t = t_sv - Δt_sv
end


"""
    Calculates the code Phase offset for time correction
    $SIGNATURES

    ´dc´: Decoder containing ephemeris data of satellite
    ´t´: Current time(i.e. raw time without correction)

    Provides the implementation of the code phase clock
    correction term
    
    Following instructions in IS-GPS-200K: 
    20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    Equation (2)
"""
function code_phase_offset(dc::GNSSDecoderState, t)
    Δtr = relativistic_correction(dc, t)
    a_f0 = dc.data.a_f0
    a_f1 = dc.data.a_f1
    a_f2 = dc.data.a_f2
    t_oc = dc.data.t_oc
    Δt_sv = a_f0 + a_f1 * (t - t_oc) + a_f2 * (t - t_oc)^2 + Δtr
    Δt_sv = L1_correction(dc, Δt_sv)
end

"""
    Provides the implementation of the relativistic clock
    correction term    
    
    $SIGNATURES

    ´dc´: Decoder containing ephemeris data of satellite
    ´t´: Current time(i.e. raw time without correction)

    Provides the implementation of the relativistic clock
    correction term 
    
    Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.1 User Algorithm for SV Clock Correction
"""
function relativistic_correction(dc::GNSSDecoderState, t, dt = 0)
    
    F = dc.constants.F
    e = dc.data.e
    sqrtA = dc.data.sqrt_A
    Ek = calc_eccentric_anomaly(dc, t, dt)
    Δtr = F * e * sqrtA * sin(Ek)
end

"""
    Calculates raw time data by observed times
    $SIGNATURES

    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    Calculates time by number of bits since TOW, TOW and code and carrier phases

"""
function calc_uncorrected_time(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    t_tow = dc.data.TOW * 6
    t_bits = dc.num_bits_buffered / GNSSSignals.get_data_frequency(GPSL1) * Hz
    t_code_phase = code_phase / GNSSSignals.get_code_frequency(GPSL1) * Hz
    t_carrier_phase = carrier_phase / GNSSSignals.get_center_frequency(GPSL1) * Hz
    
    t = t_tow + t_bits + t_code_phase + t_carrier_phase
end

"""
    Provides the implementation of the Group Delay correction
    $SIGNATURES
    ´dc´: Decoder
    ´Δt_sv´: code phase offset term
    
    Following instructions in IS-GPS-200K: 
    20.3.3.3.3.2 L1 or L2 Correction.
"""
function L1_correction(dc::GNSSDecoderState, Δt_sv)
   
    Δt_sv_L1 = Δt_sv - dc.data.T_GD
    return Δt_sv_L1
end

"""
Calculates the eccentric anomaly
$SIGNATURES
´dc´: Decoder storing ephemeris data
´t_sv´: time of satellite vehicle
´dt_sv´: time correction for SV

This function implements the calculation of the Eccentric
anomaly, which is required for both, the calculation of the
satellite position and the time correction.

Following instructions in IS-GPS-200K: 
Table 20-IV
"""
function calc_eccentric_anomaly( dc, t_sv, dt_sv = 0 )  
    
    t = t_sv - dt_sv
    
    tk = check_crossover(t - dc.data.t_oe)
    μ = dc.constants.μ
    A = dc.data.sqrt_A^2

    n_0 = sqrt(μ / (A^3))
    n = n_0 + dc.data.Δn
    M = dc.data.M_0 + n * tk
            
            
    e = dc.data.e
    E = M
    for k = 1:30
        Et = E
        E = M + e * sin(E)
        if abs( E - Et ) <= 1e-12 
            break;
        end
    end          
    return E
end

"""
Correction for end of week crossovers
$SIGNATURES

´t´: Time to be controlled

Following instructions in IS-GPS-200K: 
20.3.3.3.3.1 User Algorithm for SV Clock Correction.
"""
function check_crossover(t)

    
    half_week = 302400.0  #Half of Week in Seconds
    if (t > half_week)
        t = t - 2.0 * half_week
    elseif (t < -half_week)
        t = t + 2.0 * half_week
    end
    return t
end

