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
    Unitful,
    Statistics,
    Dates

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
    SatelliteState{CP<:Real,D<:GNSSDecoderState,S<:AbstractGNSSSignal}

Combines the GNSS decoder state with code and carrier phase measurements for a single satellite.

# Fields
- `decoder::GNSSDecoderState`: GNSS decoder state containing decoded navigation data
- `system::AbstractGNSSSignal`: GNSS system (e.g., `GPSL1CA()`, `GalileoE1B()`)
- `code_phase::CP`: Code phase measurement
- `carrier_doppler`: Carrier Doppler frequency in Hz
- `carrier_phase::CP`: Carrier phase measurement (default: `0.0`)

# Constructors
    SatelliteState(; decoder, system, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, system, sat_state)

The second constructor extracts code phase, carrier Doppler, and carrier phase from a
`Tracking` satellite state (`Tracking.SatState` for Tracking ≤ 1, `Tracking.TrackedSat` for
Tracking ≥ 2). It is provided by a package extension that is loaded automatically once
`Tracking` is available, so `Tracking` is only a weak dependency of this package.
"""
@kwdef struct SatelliteState{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    system::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
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
    calc_pvt(states::AbstractVector{<:SatelliteState},
             prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS satellite measurements.

Requires at least 4 healthy satellites from the same GNSS system. Uses least-squares
estimation for position and time, and solves for velocity and clock drift from
carrier Doppler measurements.

Unless disabled via `enable_ionospheric_correction`, the ionospheric delay is
corrected automatically using only the coefficients decoded from the navigation
messages. A single model is chosen for the whole solve and applied to every
satellite: NTCM-G if Galileo Effective Ionisation Level coefficients have been
decoded (the more accurate model), otherwise the Klobuchar model if GPS L1 α/β
have been decoded, otherwise no correction. See
[`select_ionospheric_correction`](@ref) and [`ionospheric_delay`](@ref).

Unless disabled via `enable_tropospheric_correction`, the tropospheric delay is
corrected with the blind Saastamoinen model (no broadcast coefficients needed).
See [`tropospheric_delay`](@ref).

# Arguments
- `states`: Vector of [`SatelliteState`](@ref) for observed satellites
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin)

# Keyword Arguments
- `approximate_year`: Calendar year of the observation, used to resolve the
  GPS L1 C/A 1024-week rollover ambiguity (legacy LNAV broadcasts only a
  10-bit week number, so the receiver needs external information to
  determine which 1024-week cycle the recording is in). Anything within
  ±9 years of the actual observation date works. Defaults to the current
  UTC year, which is correct for live signals; for processing archived
  recordings, pass the rough year of the recording.
- `enable_ionospheric_correction`: when `true` (default), apply the automatic
  ionospheric correction described above. Set to `false` to skip it entirely
  and solve from the raw pseudoranges (e.g. for diagnostics or when an external
  correction is applied elsewhere).
- `enable_tropospheric_correction`: when `true` (default), apply the Saastamoinen
  tropospheric correction. Set to `false` to skip it.

# Returns
A [`PVTSolution`](@ref) containing position, velocity, time, DOP values, and
satellite information. Returns `prev_pvt` if fewer than 4 healthy satellites are
available or if the GDOP is negative.

# Throws
- `ArgumentError`: If fewer than 4 satellite states are provided
"""
function calc_pvt(
    states::AbstractVector{<:SatelliteState},
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    length(states) < 4 &&
        throw(ArgumentError("You'll need at least 4 satellites to calculate PVT"))
    all(state -> state.system == states[1].system, states) ||
        ArgumentError("For now all satellites need to be base on the same GNSS")
    system = first(states).system
    healthy_indices = findall(x -> is_sat_healthy(x.decoder), states)
    length(healthy_indices) < 4 && return prev_pvt
    healthy_states = view(states, healthy_indices)
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
    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model is chosen for the whole solve from the
    # coefficients decoded across the constellation (NTCM-G if Galileo coefficients
    # are available, else Klobuchar if GPS α/β are available, else none; see
    # `select_ionospheric_correction`); the troposphere uses the blind Saastamoinen
    # model (see `tropospheric_delay`). A single corrected solve is enough: the
    # delays depend on position only through the satellite elevation/azimuth (and,
    # for the troposphere, the user height), and ∂delay/∂position is negligible over
    # the metre-level position uncertainty (a 15 m shift moves the elevation by
    # ~1e-5°), so delays predicted at a nearby position are accurate to well under a
    # millimetre — no iterate-to-convergence needed.
    correction =
        enable_ionospheric_correction ? select_ionospheric_correction(healthy_states) :
        nothing
    function predict_atmospheric_delays(position)
        user_pos = ECEF(position[1], position[2], position[3])
        # The user position is the same for every satellite, so the geodetic
        # coordinates and the ENU transform are computed once per epoch and reused.
        user_lla = LLAfromECEF(wgs84)(user_pos)
        enu_from_ecef = ENUfromECEF(user_pos, wgs84)
        map(healthy_states, sat_positions) do state, sat_pos
            elevation, azimuth = _elevation_azimuth(enu_from_ecef, sat_pos)
            iono = ionospheric_delay(
                correction,
                state.system,
                elevation,
                azimuth,
                user_lla,
                reference_time,
            )
            tropo =
                enable_tropospheric_correction ?
                tropospheric_delay(elevation, user_lla) : 0.0
            iono + tropo
        end
    end
    ξ, rmse = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, rmse_uncorrected =
            user_position(sat_positions, pseudo_ranges, prev_ξ)
        atmospheric_delays = predict_atmospheric_delays(ξ_uncorrected)
        any(!iszero, atmospheric_delays) ?
        user_position(sat_positions, pseudo_ranges .- atmospheric_delays, ξ_uncorrected) :
        (ξ_uncorrected, rmse_uncorrected)
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction.
        user_position(
            sat_positions,
            pseudo_ranges .- predict_atmospheric_delays(prev_ξ),
            prev_ξ,
        )
    end
    sat_positions_mat = reduce(hcat, sat_positions)
    H = calc_H(sat_positions_mat, ξ)
    user_velocity_and_clock_drift = calc_user_velocity_and_clock_drift(
        sat_positions_and_velocities, ξ, healthy_states, times, H)
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

    week = get_week(first(healthy_states).decoder; approximate_year)
    start_time = get_system_start_time(first(healthy_states).decoder)
    time = TAIEpoch(
        week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second,
        corrected_reference_time - floor(Int, corrected_reference_time),
    )

    sat_infos = SatInfo.(sat_positions, times)

    dop = calc_DOP(H)
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

"""
    get_week(decoder::GNSSDecoderState{<:GPSL1Data}; approximate_year)

Return the absolute GPS week number for a GPSL1 decoder, resolving the
1024-week rollover ambiguity using `approximate_year` as a calendar
anchor.

The legacy GPS L1 C/A LNAV message broadcasts only a 10-bit week number
(0–1023) modulo 1024, so the receiver cannot determine which 1024-week
cycle the recording is in from the data alone (IS-GPS-200, §20.3.3.3).
Each cycle is ~19.6 years, so any anchor within ±9 years of the true
observation date selects the correct cycle.

GPS week 0 is 1980-01-06; cycle boundaries fall on 1999-08-22,
2019-04-07, 2038-11-21, 2058-07-08, …

For Galileo, the broadcast WN is 12 bits and does not need this
treatment in any practical operational scenario.
"""
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1Data};
    approximate_year::Integer = year(now(UTC)),
)
    # GPS week 0 begins 1980-01-06. Compute the integer week count from
    # there to mid-`approximate_year`, then choose the cycle base such
    # that `cycle_base + trans_week` is closest to that anchor.
    days_at_anchor = Date(approximate_year, 6, 30) - Date(1980, 1, 6)
    weeks_at_anchor = Dates.value(days_at_anchor) ÷ 7
    n_cycles = round(Int, (weeks_at_anchor - decoder.data.trans_week) / 1024)
    return n_cycles * 1024 + decoder.data.trans_week
end

function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData};
    approximate_year::Integer = year(now(UTC)),
)
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
include("ionosphere.jl")
include("troposphere.jl")
end
