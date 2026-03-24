module PositionVelocityTime
using CoordinateTransformations,
    DocStringExtensions,
    Geodesy,
    GNSSDecoder,
    GNSSSignals,
    LinearAlgebra,
    AstroTime,
    LsqFit,
    StaticArrays,
    Tracking,
    Unitful,
    Statistics

using Unitful: s, Hz

const SPEEDOFLIGHT = 299792458.0

export calc_pvt,
    PVTSolution,
    SatelliteState,
    get_LLA,
    get_num_used_sats,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    get_sat_enu,
    get_gdop,
    get_pdop,
    get_hdop,
    get_vdop,
    get_tdop,
    get_frequency_offset

"""
    SatelliteState{CP<:Real}

Combines the GNSS decoder state with code and carrier phase measurements for a single satellite.

# Fields
- `decoder::GNSSDecoderState`: GNSS decoder state containing decoded navigation data
- `system::AbstractGNSS`: GNSS system (e.g., `GPSL1()`, `GalileoE1B()`)
- `code_phase::CP`: Code phase measurement
- `carrier_doppler`: Carrier Doppler frequency in Hz
- `carrier_phase::CP`: Carrier phase measurement (default: `0.0`)

# Constructors
    SatelliteState(; decoder, system, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, system, sat_state::SatState)

The second constructor extracts code phase, carrier Doppler, and carrier phase from a
`Tracking.SatState`.
"""
@kwdef struct SatelliteState{CP<:Real}
    decoder::GNSSDecoder.GNSSDecoderState
    system::AbstractGNSS
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
end

function SatelliteState(
    decoder::GNSSDecoder.GNSSDecoderState,
    system::AbstractGNSS,
    sat_state::SatState,
)
    SatelliteState(
        decoder,
        system,
        get_code_phase(sat_state),
        get_carrier_doppler(sat_state),
        get_carrier_phase(sat_state),
    )
end

"""
    DOP

Dilution of Precision (DOP) values describing the geometric quality of the satellite
constellation used for a PVT solution.

# Fields
- `GDOP::Float64`: Geometric DOP (overall quality)
- `PDOP::Float64`: Position DOP (3D position quality)
- `VDOP::Float64`: Vertical DOP
- `HDOP::Float64`: Horizontal DOP
- `TDOP::Float64`: Time DOP
"""
struct DOP
    GDOP::Float64
    PDOP::Float64
    VDOP::Float64
    HDOP::Float64
    TDOP::Float64
end

struct SatInfo
    position::ECEF
    time::Float64
end

"""
    PVTSolution

Complete Position, Velocity, and Time solution from GNSS measurements.

# Fields
- `position::ECEF`: User position in ECEF coordinates (meters)
- `velocity::ECEF`: User velocity in ECEF coordinates (m/s)
- `time_correction::Float64`: Estimated receiver clock bias (meters)
- `time::Union{TAIEpoch{Float64}, Nothing}`: Estimated time as a TAI epoch
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `sats::Dict{Int, SatInfo}`: Dictionary mapping PRN to satellite info (position and time)
"""
@kwdef struct PVTSolution
    position::ECEF = ECEF(0, 0, 0)
    velocity::ECEF = ECEF(0, 0, 0)
    time_correction::Float64 = 0
    time::Union{TAIEpoch{Float64},Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    sats::Dict{Int,SatInfo} = Dict{Int,SatInfo}()
end

"""
    get_num_used_sats(pvt_solution::PVTSolution) -> Int

Return the number of satellites used in the PVT solution.
"""
function get_num_used_sats(pvt_solution::PVTSolution)
    length(pvt_solution.sats)
end

"""
    get_gdop(pvt_sol::PVTSolution) -> Float64

Return the Geometric Dilution of Precision (GDOP) from the PVT solution.
"""
function get_gdop(pvt_sol::PVTSolution)
    return pvt_sol.dop.GDOP
end

"""
    get_pdop(pvt_sol::PVTSolution) -> Float64

Return the Position Dilution of Precision (PDOP) from the PVT solution.
"""
function get_pdop(pvt_sol::PVTSolution)
    return pvt_sol.dop.PDOP
end

"""
    get_vdop(pvt_sol::PVTSolution) -> Float64

Return the Vertical Dilution of Precision (VDOP) from the PVT solution.
"""
function get_vdop(pvt_sol::PVTSolution)
    return pvt_sol.dop.VDOP
end

"""
    get_hdop(pvt_sol::PVTSolution) -> Float64

Return the Horizontal Dilution of Precision (HDOP) from the PVT solution.
"""
function get_hdop(pvt_sol::PVTSolution)
    return pvt_sol.dop.HDOP
end

"""
    get_tdop(pvt_sol::PVTSolution) -> Float64

Return the Time Dilution of Precision (TDOP) from the PVT solution.
"""
function get_tdop(pvt_sol::PVTSolution)
    return pvt_sol.dop.TDOP
end

"""
    get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF) -> Spherical

Convert satellite position to East-North-Up (ENU) spherical coordinates (azimuth and
elevation) relative to the user position.

# Arguments
- `user_pos_ecef::ECEF`: User position in ECEF coordinates
- `sat_pos_ecef::ECEF`: Satellite position in ECEF coordinates

# Returns
Spherical coordinates containing azimuth and elevation of the satellite as seen from
the user position.
"""
function get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF)
    sat_enu = ENUfromECEF(user_pos_ecef, wgs84)(sat_pos_ecef)
    SphericalFromCartesian()(sat_enu)
end

"""
    calc_pvt(states::AbstractVector{<:SatelliteState}, prev_pvt::PVTSolution=PVTSolution()) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS satellite measurements.

Requires at least 4 healthy satellites from the same GNSS system. Uses least-squares
estimation for position and time, and solves for velocity and clock drift from
carrier Doppler measurements.

# Arguments
- `states`: Vector of [`SatelliteState`](@ref) for observed satellites
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin)

# Returns
A [`PVTSolution`](@ref) containing position, velocity, time, DOP values, and
satellite information. Returns `prev_pvt` if fewer than 4 healthy satellites are
available or if the GDOP is negative.

# Throws
- `ArgumentError`: If fewer than 4 satellite states are provided
"""
function calc_pvt(
    states::AbstractVector{<:SatelliteState},
    prev_pvt::PVTSolution = PVTSolution(),
)
    length(states) < 4 &&
        throw(ArgumentError("You'll need at least 4 satellites to calculate PVT"))
    all(state -> state.system == states[1].system, states) ||
        ArgumentError("For now all satellites need to be base on the same GNSS")
    system = first(states).system
    healthy_states = filter(x -> is_sat_healthy(x.decoder), states)
    if length(healthy_states) < 4
        return prev_pvt
    end
    prev_ξ = [prev_pvt.position; prev_pvt.time_correction]
    healthy_prns = map(state -> state.decoder.prn, healthy_states)
    times = map(state -> calc_corrected_time(state), healthy_states)
    sat_positions_and_velocities = map(
        (state, time) -> calc_satellite_position_and_velocity(state.decoder, time),
        healthy_states,
        times,
    )
    sat_positions = map(get_sat_position, sat_positions_and_velocities)
    pseudo_ranges, reference_time = calc_pseudo_ranges(times)
    ξ, rmse = user_position(sat_positions, pseudo_ranges)
    user_velocity_and_clock_drift =
        calc_user_velocity_and_clock_drift(sat_positions_and_velocities, ξ, healthy_states)
    position = ECEF(ξ[1], ξ[2], ξ[3])
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    relative_clock_drift = user_velocity_and_clock_drift[4] / SPEEDOFLIGHT
    time_correction = ξ[4]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEEDOFLIGHT

    week = get_week(first(healthy_states).decoder)
    start_time = get_system_start_time(first(healthy_states).decoder)
    time = TAIEpoch(
        week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second,
        corrected_reference_time - floor(Int, corrected_reference_time),
    )

    sat_infos = SatInfo.(sat_positions, times)

    dop = calc_DOP(calc_H(reduce(hcat, sat_positions), ξ))
    if dop.GDOP < 0
        return prev_pvt
    end

    PVTSolution(
        position,
        velocity,
        time_correction,
        time,
        relative_clock_drift,
        dop,
        Dict(healthy_prns .=> sat_infos),
    )
end

"""
    get_frequency_offset(pvt::PVTSolution, base_frequency) -> typeof(base_frequency)

Calculate the receiver frequency offset from the relative clock drift and a base frequency.

# Arguments
- `pvt::PVTSolution`: PVT solution containing the relative clock drift
- `base_frequency`: Reference frequency (e.g., the carrier frequency of the GNSS signal)

# Returns
The frequency offset as `relative_clock_drift * base_frequency`.
"""
function get_frequency_offset(pvt::PVTSolution, base_frequency)
    pvt.relative_clock_drift * base_frequency
end

function get_system_start_time(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data},
)
    TAIEpoch(1980, 1, 6, 0, 0, 19.0) # There were 19 leap seconds at 01/06/1999 compared to UTC
end

function get_system_start_time(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
)
    TAIEpoch(1999, 8, 22, 0, 0, (32 - 13.0)) # There were 32 leap seconds at 08/22/1999 compared to UTC
end

function get_week(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data})
    2048 + decoder.data.trans_week
end

function get_week(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData})
    decoder.data.WN
end

"""
    get_LLA(pvt::PVTSolution) -> LLA

Convert the ECEF position in the PVT solution to geodetic coordinates
(latitude, longitude, altitude) using the WGS84 ellipsoid.
"""
function get_LLA(pvt::PVTSolution)
    LLAfromECEF(wgs84)(pvt.position)
end

include("user_position.jl")
include("sat_time.jl")
include("sat_position.jl")
end
