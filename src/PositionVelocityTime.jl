module PositionVelocityTime
    using CoordinateTransformations, DocStringExtensions, Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, Parameters
    using Unitful: s, Hz
    
    export  calc_PVT,
            is_sat_healthy_and_decodable,
            PVTSolution,
            sat_position_rotation, 
            sat_position_ECI_2_ECEF, 
            sat_position_ECEF,
            user_position,
            SatelliteState,
            get_num_used_sats,
            get_gdop,
            get_pdop,
            get_hdop,
            get_vdop,
            get_tdop
            
    
    """
    Struct of decoder, code- and carrierphase of satellite vehicle
    """
    @with_kw struct SatelliteState{CP <: Real}
        decoder_state::GNSSDecoderState
        code_phase::CP
        carrier_phase::CP = 0
    end

    """
    DOP, stores the Dilution of Precision values
    """
    @with_kw mutable struct DOP
        GDOP::Union{Nothing, Float64} = nothing
        PDOP::Union{Nothing, Float64} = nothing
        VDOP::Union{Nothing, Float64} = nothing
        HDOP::Union{Nothing, Float64} = nothing
        TDOP::Union{Nothing, Float64} = nothing
    end


    """
    PVTSolution, with used sats and el and az.
    """
    @with_kw mutable struct PVTSolution
        pos::ECEF = ECEF([0,0,0])
        receiver_time_correction::Float64 = 0
        DOP_val::DOP = DOP()
        #following is added for satllite observation
        used_sats::AbstractVector{Int64} = [] 
        elevation_sats::AbstractVector{Float64} = []
        azimuth_sats::AbstractVector{Float64} = [] 
    end

    function get_num_used_sats(pvt_solution::PVTSolution)
        length(pvt_solution.used_sats)
    end

    #Get methods for single DOP values
    function get_gdop(pvt_sol::PVTSolution)
        return pvt_sol.DOP_val.GDOP
    end

    function get_pdop(pvt_sol::PVTSolution)
        return pvt_sol.DOP_val.PDOP
    end

    function get_vdop(pvt_sol::PVTSolution)
        return pvt_sol.DOP_val.VDOP
    end

    function get_hdop(pvt_sol::PVTSolution)
        return pvt_sol.DOP_val.HDOP
    end

    function get_tdop(pvt_sol::PVTSolution)
        return pvt_sol.DOP_val.TDOP
    end

    """
    Calculates elevation and azimuth of satellite

    $SIGNATURES
    user_pos_ecef: user position in ecef coordinates
    sat_pos_ecef: satellite position in ecef coordinates

    This function calculates the elevation and azimuth of the satellite in degree
    """

    function get_elevation_azimuth(user_pos_ecef::ECEF, sat_pos_ecef::ECEF)
        sat_enu = ENUfromECEF(user_pos_ecef, wgs84)(sat_pos_ecef)
        sat_enu_sph = SphericalFromCartesian()(sat_enu)
        az = rad2deg(sat_enu_sph.θ)
        el = rad2deg(sat_enu_sph.ϕ)
        return el, az
    end

    """
    Calculates ECEF position of user

    $SIGNATURES
    ´sat_state´: satellite state, combining decoded data, code- and carrierphase 
    
    ´min_elevation´: The minimum elevation for a satellite, that can be used for a positioning
                     Hint: The elevation starts here with 0° at the horizon and has its maximum with 90° in the zenith.
                     min_elevation = -100, deactivates the  elevation filter
    

    This function calculates the position of the user in ECEF coordinates
        incl. earth rotation correction

  
    The implementation follows IS-GPS-200K Table 20-IV.
    """
    function calc_PVT(
        satellite_states::AbstractVector{SatelliteState{Float64}},
        prev_pvt::PVTSolution = PVTSolution(),
        min_elevation::Int64 = -100 #elevation filter is normally turned off
        )

        c = satellite_states[1].decoder_state.constants.c
        usable_satellite_states = filter(x -> is_sat_healthy_and_decodable(x.decoder_state), satellite_states)
        length(usable_satellite_states) >= 4 || throw(SmallData("Not enough usable SV Data"))
        sv_positions = map(x -> sat_position_ECEF(x), usable_satellite_states)
        
        if get_num_used_sats(prev_pvt) == 0
            pseudoranges = pseudo_ranges(usable_satellite_states)
            pvt_sol = user_position([0,0,0], sv_positions, pseudoranges)
            prev_pos = pvt_sol.pos
        else
            prev_pos = prev_pvt.pos
        end

        # Remove Satellite with low elevation angle
        if min_elevation >= 0 && min_elevation < 90 # min elevaiton filter
            out = map( x-> get_elevation_azimuth(prev_pos, ECEF(x[1],x[2],x[3])), sv_positions)
            el, az = map(x->getindex.(out,x), 1:2)
            el_filter = el .>= min_elevation
            usable_satellite_states = usable_satellite_states[el_filter]
            el_save = el[el_filter]
            az_save = az[el_filter]
            PRN_save = map(x -> x.decoder_state.PRN, usable_satellite_states)
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
