
"""
Computes psuedo ranges 

$SIGNATURES
´decoder_state´: Decoder containing ephemeris data of satellite
´code_phase´: Code phase of signal at time of measure
´carrier_phase´: Center frequency of carrier signal

Computes relative pseudo ranges of satellite vehicles.
The algorithm is based on the common reception method. 
"""
function pseudo_ranges(
    decoder_states::AbstractVector{GNSSDecoderState}, 
    code_phases::AbstractVector, 
    carrier_phases::AbstractVector = zeros(length(decoder_states))
    )
    
    for i in 1:length(decoder_states)
        is_sat_healthy_and_decodable(decoder_states[i]) || throw(BadData("SV not decoded properly, #decoder: " * string(i)))
    end
    (length(decoder_states) == length(code_phases) == length(carrier_phases)) || throw(IncompatibleData("Length of Arrays of Decoder, Code Phases and Carrier Phases must be equal"))


    N_sats = length(decoder_states)
    c = decoder_states[1].constants.c
    
    times = map(i -> calc_corrected_time(decoder_states[i], code_phases[i], carrier_phases[i]), 1:N_sats)
    t_ref = maximum(times)

    τ = map(i -> t_ref - times[i], 1:N_sats)
    pseudoranges = τ .* c
    return pseudoranges
end
