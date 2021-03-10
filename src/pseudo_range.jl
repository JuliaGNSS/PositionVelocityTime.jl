
"""
Computes psuedo ranges 

$SIGNATURES
<<<<<<< HEAD
´dc´: Decoder containing ephemeris data of satellite
´code_phase´: Code phase of signal at time of measure
´carrier_phase´: Center frequency of carrier signal
=======
´sat_state´: satellite state, combining decoded data, code- and carrierphase 
>>>>>>> dev

Computes relative pseudo ranges of satellite vehicles.
The algorithm is based on the common reception method. 
"""
function pseudo_ranges(
    sat_states::AbstractVector{SatelliteState{Float64}}
    )
    
    for i in 1:length(sat_states)
        is_sat_healthy_and_decodable(sat_states[i].decoder_state) || throw(BadData("SV not decoded properly, #decoder: " * string(i)))
    end


    N_sats = length(sat_states)
    c = sat_states[1].decoder_state.constants.c
    
    times = map(i -> calc_corrected_time(sat_states[i]), 1:N_sats)

    t_ref = maximum(times)

    τ = map(i -> t_ref - times[i], 1:N_sats)
    pseudoranges = τ .* c
    return pseudoranges
end
