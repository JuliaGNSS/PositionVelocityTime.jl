module PositionVelocityTime

    using DocStringExtensions, Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, Parameters
    using Unitful: s, Hz
    
    export  calc_PVT,
            calc_PVT_old,
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
    Calculates elevation and azimuth of satellite

    $SIGNATURES
    user_pos_ecef: user position in ecef coordinates
    sat_pos_ecef: satellite position in ecef coordinates

    This function calculates the elevation and azimuth of the satellite in degree
    The implementation follows ESA_GNSS-Book_TM-23_Vol_I B.3
    """

    function get_elevation_azimuth(user_pos_ecef::ECEF, sat_pos_ecef)
        
        p = sat_pos_ecef-user_pos_ecef
        dist = norm(p)
        p̂ = p/dist
        
        user_lla = LLAfromECEF(wgs84)(ECEF(user_pos_ecef[1], user_pos_ecef[2], user_pos_ecef[3]))
        lat = user_lla.lat
        lon = user_lla.lon

        ê = [-sind(lat), cosd(lat), 0]
        n̂ = [-cosd(lon)*sind(lat), -sind(lat)*sind(lon), cosd(lat)]
        û = [cosd(lon)*cosd(lat), sind(lon)*cosd(lat), sind(lat)]

        el = asind(dot(p̂,û))
        az = atand(dot(p̂,ê),dot(p̂,n̂))

        return el, az

    end

    """
    Calculates ECEF position of user

    $SIGNATURES
    ´sat_state´: satellite state, combining decoded data, code- and carrierphase 

    This function calculates the position of the user in ECEF coordinates
    The implementation follows IS-GPS-200K Table 20-IV.
    """

    function calc_PVT_old(
        satellite_states::AbstractVector{SatelliteState{Float64}}
        )
        c = satellite_states[1].decoder_state.constants.c

        usable_satellite_states = filter(x -> is_sat_healthy_and_decodable(x.decoder_state), satellite_states)
        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data"))
        sv_positions = map(x -> sat_position_ECEF(x), usable_satellite_states)
        pseudoranges = pseudo_ranges(usable_satellite_states)

        userpos = user_position(sv_positions, pseudoranges)
    end

    """
    Calculates ECEF position of user

    $SIGNATURES
    ´sat_state´: satellite state, combining decoded data, code- and carrierphase 
    
    ´min_elevation´: the minimum elevation for a satellite, that can be used for a positioning
                    min_elevation = -100, deactivates the  elevation filter
    

    This function calculates the position of the user in ECEF coordinates
        incl. earth rotation correction

  
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(
        satellite_states::AbstractVector{SatelliteState{Float64}},
        min_elevation::Int64 = -100 #elevation filter is normally turned off
        )
        c = satellite_states[1].decoder_state.constants.c

        usable_satellite_states = filter(x -> is_sat_healthy_and_decodable(x.decoder_state), satellite_states)
        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data"))
        sv_positions = map(x -> sat_position_ECEF(x), usable_satellite_states)
        pseudoranges = pseudo_ranges(usable_satellite_states)

        userpos = user_position(sv_positions, pseudoranges)


        # Remove Satellite with low elevation angle (new)
        if min_elevation != -100  # min elevaiton filter
            out = map( x-> get_elevation_azimuth(userpos.pos, x), sv_positions)
            el, az = map(x->getindex.(out,x), 1:2)

            println(el)
            println(az)
            println("vec length sv_positions: ", length(sv_positions))
            map(x-> print(x.decoder_state.PRN,", "), usable_satellite_states)
            println("\n")

            buf_state = usable_satellite_states
            buf_sat_pos = sv_positions
            for i = 1:length(usable_satellite_states)
                if el[i] < min_elevation
                    #println("Filter PRN:", buf_state[i].decoder_state.PRN)
                    usable_satellite_states = filter(x -> buf_state[i].decoder_state.PRN != x.decoder_state.PRN, usable_satellite_states)
                    sv_positions = filter(x -> x != buf_sat_pos[i], sv_positions)
                end
            end
        end

        #new part: compensation of earth rotation

        dt =  map(x -> norm(x-userpos.pos)/c, sv_positions)
        sv_positions = map(i -> sat_position_ECEF2ECI(usable_satellite_states[i],-dt[i]), 1:length(usable_satellite_states))
        
        pseudoranges = pseudo_ranges(usable_satellite_states)
        userpos = user_position(sv_positions, pseudoranges)

    end

    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
