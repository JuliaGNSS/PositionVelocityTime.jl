module PositionVelocityTime

    using DocStringExtensions, Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, Parameters
    using Unitful: s, Hz
    
    export  calc_PVT,
            is_sat_healthy_and_decodable,
            PVTSolution,
            sat_position_ECI_2_ECEF, 
            sat_position_ECEF,
            user_position,
            SatelliteState
            
    
    """
    Struct of decoder, code- and carrierphase of satellite vehicle
    """
    @with_kw struct SatelliteState{CP <: Real}
        decoder_state::GNSSDecoderState
        code_phase::CP
        carrier_phase::CP = 0
    end


    """
    Calculates ECEF position of user

    $SIGNATURES
    ´sat_state´: satellite state, combining decoded data, code- and carrierphase 

    This function calculates the position of the user in ECEF coordinates
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(
        satellite_states::AbstractVector{SatelliteState{Float64}}
        )
        
        
        usable_satellite_states = filter(x -> is_sat_healthy_and_decodable(x.decoder_state), satellite_states)

        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data"))
        sv_positions = map(x -> sat_position_ECEF(x), usable_satellite_states)
        pseudoranges = pseudo_ranges(usable_satellite_states)

        userpos = user_position(sv_positions, pseudoranges)

    end
    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
