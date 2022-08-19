module PositionVelocityTime
    using CoordinateTransformations,
        DocStringExtensions,
        Geodesy,
        GNSSDecoder,
        GNSSSignals,
        LinearAlgebra,
        Parameters,
        AstroTime,
        LsqFit,
        StaticArrays,
        Tracking

    using Unitful: s, Hz

    const SPEEDOFLIGHT = 299792458.0
    
    export  calc_pvt,
        PVTSolution,
        SatelliteState,
        get_LLA,
        get_num_used_sats,
        calc_satellite_position,
        get_sat_enu,
        get_gdop,
        get_pdop,
        get_hdop,
        get_vdop,
        get_tdop
            
    
    """
    Struct of decoder, code- and carrierphase of satellite
    """
    @with_kw struct SatelliteState{CP <: Real}
        decoder::GNSSDecoder.GNSSDecoderState
        system::AbstractGNSS
        code_phase::CP
        carrier_phase::CP = 0.0
    end

    function SatelliteState(
        decoder::GNSSDecoder.GNSSDecoderState,
        tracking_results::Tracking.TrackingResults
    )
        SatelliteState(
            decoder,
            get_system(tracking_results),
            get_code_phase(tracking_results),
            get_carrier_phase(tracking_results)
        )
    end

    """
    Dilution of Precision (DOP)
    """
    struct DOP
        GDOP::Float64
        PDOP::Float64
        VDOP::Float64
        HDOP::Float64
        TDOP::Float64
    end


    """
    PVT solution including DOP, used satellites and satellite
    positions.
    """
    @with_kw struct PVTSolution
        position::ECEF = ECEF(0, 0, 0)
        time_correction::Float64 = 0
        time::Union{TAIEpoch{Float64}, Nothing} = nothing
        dop::Union{DOP, Nothing} = nothing
        used_sats::Vector{Int64} = [] 
        sat_positions::Vector{ECEF} = []
    end

    function get_num_used_sats(pvt_solution::PVTSolution)
        length(pvt_solution.used_sats)
    end

    #Get methods for single DOP values
    function get_gdop(pvt_sol::PVTSolution)
        return pvt_sol.dop.GDOP
    end

    function get_pdop(pvt_sol::PVTSolution)
        return pvt_sol.dop.PDOP
    end

    function get_vdop(pvt_sol::PVTSolution)
        return pvt_sol.dop.VDOP
    end

    function get_hdop(pvt_sol::PVTSolution)
        return pvt_sol.dop.HDOP
    end

    function get_tdop(pvt_sol::PVTSolution)
        return pvt_sol.dop.TDOP
    end

    """
    Calculates East-North-Up (ENU) coordinates of satellite in spherical
    form (azimuth and elevation).

    $SIGNATURES
    user_pos_ecef: user position in ecef coordinates
    sat_pos_ecef: satellite position in ecef coordinates
    """
    function get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF)
        sat_enu = ENUfromECEF(user_pos_ecef, wgs84)(sat_pos_ecef)
        SphericalFromCartesian()(sat_enu)
    end

    """
    Calculates Position Velocity and Time (PVT).
    Note: Estimation of velocity still needs to be implemented.

    $SIGNATURES
    `system`: GNSS system
    `states`: Vector satellite states (SatelliteState)
    `prev_pvt` (optionally): Previous PVT solution to accelerate calculation of next
        PVT.
    """
    function calc_pvt(
        states::AbstractVector{<: SatelliteState},
        prev_pvt::PVTSolution = PVTSolution()
    )
        length(states) < 4 && throw(ArgumentError("You'll need at least 4 satellites to calculate PVT"))
        all(state -> state.system == states[1].system, states) || ArgumentError("For now all satellites need to be base on the same GNSS")
        healthy_states = filter(x -> is_sat_healthy(x.decoder), states)
        if length(healthy_states) < 4
            return prev_pvt
        end
        prev_ξ = [prev_pvt.position; prev_pvt.time_correction]
        healthy_prns = map(state -> state.decoder.prn, healthy_states)
        times = map(state -> calc_corrected_time(state), healthy_states)
        sat_positions = map((state, time) -> calc_satellite_position(state.decoder, time), healthy_states, times)
        pseudo_ranges, reference_time = calc_pseudo_ranges(times)
        ξ = user_position(sat_positions, pseudo_ranges, prev_ξ)
        position = ECEF(ξ[1], ξ[2], ξ[3])
        time_correction = ξ[4]
        corrected_reference_time = reference_time + time_correction / SPEEDOFLIGHT
    
        week = get_week(first(healthy_states).decoder)
        start_time = get_system_start_time(first(healthy_states).decoder)
        time = TAIEpoch(week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second, corrected_reference_time - floor(Int, corrected_reference_time))
    
        dop = calc_DOP(H(reduce(hcat, sat_positions), ξ))
        
        PVTSolution(position, time_correction, time, dop, healthy_prns, sat_positions)
    end

    function get_system_start_time(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data})
        TAIEpoch(1980, 1, 6, 0, 0, 19.0) # There were 19 leap seconds at 01/06/1999 compared to UTC
    end

    function get_system_start_time(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData})
        TAIEpoch(1999, 8, 22, 0, 0, (32 - 13.0)) # There were 32 leap seconds at 08/22/1999 compared to UTC
    end

    function get_week(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data})
        2048 + decoder.data.trans_week
    end
    
    function get_week(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData})
        decoder.data.WN
    end

    function get_LLA(pvt::PVTSolution)
        LLAfromECEF(wgs84)(pvt.position)
    end

    
    include("user_position.jl")
    include("sat_time.jl")
    include("sat_position.jl")
end
