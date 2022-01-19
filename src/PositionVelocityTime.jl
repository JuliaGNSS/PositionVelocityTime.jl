module PositionVelocityTime

    using DocStringExtensions, Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, Parameters
    using Unitful: s, Hz
    
    export  calc_PVT,
            is_sat_healthy_and_decodable,
            PVTSolution,
            sat_position_rotation, 
            sat_position_ECI_2_ECEF, 
            sat_position_ECEF,
            user_position,
            null_PVT,
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
    PVTSolution, with used sats and el and az.
    """
    @with_kw mutable struct PVTSolution
        pos::ECEF
        receiver_time_correction::Float64
        GDOP::Union{Nothing, Float64}
        PDOP::Union{Nothing, Float64} 
        VDOP::Union{Nothing, Float64} 
        HDOP::Union{Nothing, Float64} 
        TDOP::Union{Nothing, Float64}
        #following is added for satllite observation
        num_used_sats::Int64
        used_sats::Union{Nothing,AbstractVector{Int64}} 
        elevation_sats::Union{Nothing,AbstractVector{Float64}}  
        azimuth_sats::Union{Nothing,AbstractVector{Float64}}  
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
    
    ´min_elevation´: the minimum elevation for a satellite, that can be used for a positioning
                    min_elevation = -100, deactivates the  elevation filter
    

    This function calculates the position of the user in ECEF coordinates
        incl. earth rotation correction

  
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(
        prev_pvt::PVTSolution,
        satellite_states::AbstractVector{SatelliteState{Float64}},
        min_elevation::Int64 = -100 #elevation filter is normally turned off
        )

        c = satellite_states[1].decoder_state.constants.c
        usable_satellite_states = filter(x -> is_sat_healthy_and_decodable(x.decoder_state), satellite_states)
        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data"))
        sv_positions = map(x -> sat_position_ECEF(x), usable_satellite_states)
        
        if prev_pvt.num_used_sats == 0
            pseudoranges = pseudo_ranges(usable_satellite_states)
            pvt_sol = user_position([0,0,0], sv_positions, pseudoranges)
            prev_pos = pvt_sol.pos
        else
            prev_pos = prev_pvt.pos
        end

        # Remove Satellite with low elevation angle
        if min_elevation >= 0 && min_elevation < 90 # min elevaiton filter
            out = map( x-> get_elevation_azimuth(prev_pos, x), sv_positions)
            el, az = map(x->getindex.(out,x), 1:2)

            buf_state = usable_satellite_states
            buf_sat_pos = sv_positions
            el_save = Float64[]
            az_save = Float64[]
            PRN_save = Int64[]
            for i = 1:length(usable_satellite_states)
                if el[i] < min_elevation
                    usable_satellite_states = filter(x -> buf_state[i].decoder_state.PRN != x.decoder_state.PRN, usable_satellite_states)
                    sv_positions = filter(x -> x != buf_sat_pos[i], sv_positions)
                else
                    append!(PRN_save, buf_state[i].decoder_state.PRN)
                    append!(el_save, el[i])
                    append!(az_save, az[i])
                end
            end
        end

        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data, with elevation filter."))

        #compensation of earth rotation within TOF
        dt =  map(x -> norm(x-prev_pos)/c, sv_positions)
        sv_positions = map(i -> sat_position_rotation(usable_satellite_states[i], dt[i]), 1:length(usable_satellite_states))

        pseudoranges = pseudo_ranges(usable_satellite_states)
        pvt = user_position(prev_pos, sv_positions, pseudoranges)
        pvt.used_sats = PRN_save
        pvt.elevation_sats = el_save
        pvt.azimuth_sats = az_save
        return pvt
    end

    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
