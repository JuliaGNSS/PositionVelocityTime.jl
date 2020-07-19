function calc_corrected_time(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    # Provides the implementation of the relativistic clock
    # correction term
    #
    # Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    # Equation (1)
    t_sv = calc_uncorrected_time(dc, code_phase, carrier_phase)
    Δt_sv = code_phase_offset(dc, t_sv)
    t = t_sv - Δt_sv
end



function code_phase_offset(dc::GNSSDecoderState, t)
    # Computes the SV code phase offset
    #
    # Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    # Equation (2)
    Δtr = relativistic_correction(dc, t)
    a_f0 = dc.data.a_f0
    a_f1 = dc.data.a_f1
    a_f2 = dc.data.a_f2
    t_oc = dc.data.t_oc
    Δt_sv = a_f0 + a_f1 * (check_crossover(dc, t - t_oc)) + a_f2 * (check_crossover(dc, t - t_oc))^2 + Δtr
    Δt_sv = L1_correction(dc, Δt_sv)
end


function relativistic_correction(dc::GNSSDecoderState, t, dt = 0)
    # Provides the implementation of the relativistic clock
    # correction term
    #
    # Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    F = dc.constants.F
    e = dc.data.e
    sqrtA = dc.data.sqrt_A
    Ek = calc_eccentric_anomaly(dc, t, dt)
    Δtr = F * e * sqrtA * sin(Ek)
end


function calc_uncorrected_time(dc::GNSSDecoderState, code_phase, carrier_phase = 0)
    t_tow = dc.data.TOW * 6
    t_bits = dc.num_bits_buffered / GNSSSignals.get_data_frequency(GPSL1) * Hz
    t_code_phase = code_phase / GNSSSignals.get_code_frequency(GPSL1) * Hz
    t_carrier_phase = carrier_phase / GNSSSignals.get_center_frequency(GPSL1) * Hz
    
    t = t_tow + t_bits + t_code_phase + t_carrier_phase
end


function L1_correction(dc::GNSSDecoderState, Δt_sv)
    # Provides the implementation of the Group Delay correction
    #
    # Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.2 L1 or L2 Correction.
    Δt_sv_L1 = Δt_sv - dc.data.T_GD
    return Δt_sv_L1
end

function calc_eccentric_anomaly( dc, t_sv, dt_sv = 0 )  
    #  This function implements the calculation of the Eccentric
    #  anomaly, which is required for both, the calculation of the
    #  satellite position and the time correction.
    #  The implementation follows IS-GPS-200K
        
    
    t = t_sv - dt_sv
    tk = t - dc.data.t_oe
    tk = check_crossover(dc, tk)
    
    
    n_0 = calc_n0(dc)
    n = n_0 + dc.data.Δn
    Mk = dc.data.M_0 + n * tk
            
            
    ec = dc.data.e
    Et = Mk
    Ek = Mk + ec .* sin(Et)
    for k = 1:50
        Et = Ek
        Ek = Mk + ec .* sin(Et)
        if abs( Et - Ek ) <= 1e-12 
            break;
        end
    end          
    return Ek
end




function calc_n0(dc::GNSSDecoderState)
    # calculates n_0 for eccentric anomaly
    #
    # Following instructions in IS-GPS-200K: 
    # Table 20-IV.  Elements of Coordinate Systems 
    μ = dc.constants.μ
    n0 = sqrt(μ / dc.data.sqrt_A ^ 6)
end

function check_crossover(dc::GNSSDecoderState, t)
    # Correction for end of week crossovers
    #
    # Following instructions in IS-GPS-200K: 
    # 20.3.3.3.3.1 User Algorithm for SV Clock Correction.
    
    #t_oc = dc.data.t_oc
    half_week = 302400.0  #Half of Week in Seconds
    if (t > half_week)
        t = t - 2.0 * half_week
    elseif (t < -half_week)
        t = t + 2.0 * half_week
    end
    return t
end

