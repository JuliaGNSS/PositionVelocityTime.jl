
"""
Computes psuedo ranges 

$SIGNATURES
´dc´: Decoder containing ephemeris data of satellite
´code_phase´: Code phase of signal at time of measure
´carrier_phase´: Center frequency of carrier signal

Computes relative pseudo ranges of satellite vehicles.
The algorithm is based on the common reception method. 
"""
function pseudo_range(dcs::Vector{GNSSDecoderState} ,code_phases::Vector{Float64}, carrier_phases = -1)
    for i in 1:length(dcs)
        can_get_sat_position(dcs[i]) || throw(BadData("SV not decoded properly, #dc: " * string(i)))
    end
    N_sats = length(dcs)
    c = dcs[1].constants.c
    if carrier_phases == -1
        (length(dcs) == length(code_phases)) || throw(IncompatibleData("Length of Arrays of decoder and Code Phases must be equal"))
        carrier_phases = zeros(N_sats)
    else
        (length(dcs) == length(code_phases) == length(carrier_phases)) || throw(IncompatibleData("Length of Arrays of Decoder, Code Phases and Carrier Phases must be equal"))
    end

    times = map(i -> calc_corrected_time(dcs[i], code_phases[i], carrier_phases[i]), 1:N_sats)
    t_ref = maximum(times)

    τ = map(i -> t_ref - times[i], 1:N_sats)
    pseudoranges = τ .* c
    return Vector(pseudoranges)
end
