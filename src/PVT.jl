module PVT

    using DocStringExtensions, Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra
    using Unitful: s, Hz
    
    export  calc_PVT,
            is_sat_healthy_and_decodable,
            PVTSolution,
            sat_position_ECI_2_ECEF, 
            sat_position_ECEF,
            user_position
            
    """
    Calculates ECEF position of user

    $SIGNATURES
    ´decoder_state´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the user in ECEF coordinates
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(
        decoder_states::AbstractVector{GNSSDecoderState}, 
        code_phases::AbstractVector, 
        carrier_phases::AbstractVector = zeros(length(decoder_states))
        )
        
       
        (length(decoder_states) == length(code_phases) == length(carrier_phases)) || throw(IncompatibleData("Length of Arrays of Decoder, Code Phases and Carrier Phases must be equal"))
        
        
        usable_sv = Vector{GNSSDecoderState}(undef, 0)
        usable_code_phases = Vector{Float64}(undef, 0)
        usable_carr_phases = Vector{Float64}(undef, 0)
        for i in 1:length(decoder_states)
            if is_sat_healthy_and_decodable(decoder_states[i])
                push!(usable_sv, decoder_states[i])
                push!(usable_code_phases, code_phases[i])
                push!(usable_carr_phases, carrier_phases[i])
            end
        end
        length(usable_sv) >= 4 || throw(SmallData("Not enough usable SV Data"))

        sv_positions = map( (decoder_states, cps, caps) -> sat_position_ECEF(decoder_states, cps, caps), usable_sv, usable_code_phases, usable_carr_phases)
        pseudoranges = pseudo_ranges(usable_sv, usable_code_phases, usable_carr_phases)

        userpos = user_position(sv_positions, pseudoranges)

    end
    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
