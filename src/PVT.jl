module PVT

    using DocStringExtensions, LinearAlgebra, GNSSDecoder, GNSSSignals
    using Unitful: s, Hz
    
    export  calc_PVT,
            sat_position_ECI_2_ECEF, 
            sat_position_ECEF,
            user_position,
            can_get_sat_position
    """
    Calculates ECEF position of user

    $SIGNATURES
    ´dc´: Decoder containing ephemeris data of satellite
    ´code_phase´: Code phase of signal at time of measure
    ´carrier_phase´: Center frequency of carrier signal

    This function calculates the position of the user in ECEF coordinates
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(dcs::Vector{GNSSDecoderState}, code_phases::Vector{Float64}, carrier_phases = zeros(length(dcs)))
        
       
        (length(dcs) == length(code_phases) == length(carrier_phases)) || throw(IncompatibleData("Length of Arrays of Decoder, Code Phases and Carrier Phases must be equal"))
        
        
        usable_sv = Vector{GNSSDecoderState}(undef, 0)
        usable_code_phases = Vector{Float64}(undef, 0)
        usable_carr_phases = Vector{Float64}(undef, 0)
        for i in 1:length(dcs)
            if is_sat_healthy_and_decodable(dcs[i])
                push!(usable_sv, dcs[i])
                push!(usable_code_phases, code_phases[i])
                push!(usable_carr_phases, carrier_phases[i])
            end
        end
        length(usable_sv) >= 4 || throw(SmallData("Not enough usable SV Data"))

        sv_positions = map( (dcs, cps, caps) -> sat_position_ECEF(dcs, cps, caps), usable_sv, usable_code_phases, usable_carr_phases)
        pseudoranges = pseudo_ranges(usable_sv, usable_code_phases, usable_carr_phases)

        userpos = user_position(sv_positions, pseudoranges)

    end
    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
