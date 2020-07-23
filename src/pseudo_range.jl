function pseudo_range(dcs::Array{GNSSDecoderState,1},code_phases::Array{Float64, 1}, carrier_phases = -1)
    
    N_sats = length(dcs)
    c = dcs[1].constants.c
    if carrier_phases == -1
        (length(dcs) == length(code_phases)) || throw(IncompatibleData("Length of Arrays of decoder and Code Phases must be equal"))
        carrier_phases = zeros(N_sats)
    else
        (length(dcs) == length(code_phases) == length(carrier_phases)) || throw(IncompatibleData("Length of Arrays of Decoder, Code Phases and Carrier Phases must be equal"))
    end

    times = map(i -> calc_corrected_time(dcs[i], code_phases[i], carrier_phases), 1:N_sats)
    t_ref = max(times)

    τ = map(i -> t_ref - times[i], 1:N_sats)
    pseudoranges = τ .* c
    return pseudoranges
end

function pseudo_range(dc::GNSSDecoderState, dc_ref::GNSSDecoderState, code_phase::Float64, code_phase_ref::Float64, carrier_phase::Float64 = 0.0, carrier_phase_ref::Float64 = 0.0)
    c = dc.constants.c
    t = calc_corrected_time(dc, code_phase, carrier_phase)
    t_ref = calc_corrected_time(dc_ref, code_phase_ref, carrier_phase_ref)
    τ = t_ref - t 
    pseudorange = τ * c 
    return pseudorange
end