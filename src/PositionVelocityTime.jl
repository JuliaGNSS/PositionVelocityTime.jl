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
    Dates

using Unitful: s, Hz

const SPEEDOFLIGHT = 299792458.0

export calc_pvt,
    PVTSolution,
    SatInfo,
    SatelliteState,
    get_LLA,
    get_num_used_sats,
    get_sat_info,
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

"""
    SatInfo

Per-satellite information attached to a [`PVTSolution`](@ref) (one entry per
satellite used in the fix).

# Fields
- `position::ECEF`: Satellite ECEF position at transmit time (metres).
- `time::Float64`: Satellite transmit time (system time of week, seconds).
- `residual::Float64`: Post-fit least-squares pseudorange residual (metres) — the
  modeled minus the (atmosphere-corrected) measured pseudorange. A per-satellite
  fit-quality / outlier indicator.
"""
struct SatInfo
    position::ECEF
    time::Float64
    residual::Float64
end

"""
    PVTSolution

Complete Position, Velocity, and Time solution from GNSS measurements.

# Fields
- `position::ECEF`: User position in ECEF coordinates (meters)
- `velocity::ECEF`: User velocity in ECEF coordinates (m/s)
- `time_correction::Float64`: Estimated receiver clock bias of the reference GNSS
  (meters). For a multi-GNSS solution this is the bias of `reference_system`;
  other systems' biases are `time_correction + inter_system_biases[system]`.
- `time::Union{TAIEpoch{Float64}, Nothing}`: Estimated time as a TAI epoch
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `sats::Dict{Tuple{Symbol, Int}, SatInfo}`: Dictionary mapping `(GNSS time system, PRN)`
  to satellite info — position, transmit time, and post-fit residual (see [`SatInfo`](@ref)).
  The system tag (`:GPS`, `:Galileo`) is what keeps the same PRN apart across
  constellations (e.g. GPS PRN 5 vs Galileo E05).
- `reference_system::Union{Symbol, Nothing}`: GNSS time system (e.g. `:GPS`,
  `:Galileo`) that `time` and `time_correction` are referenced to.
- `inter_system_biases::Dict{Symbol, Float64}`: For each GNSS time system other
  than `reference_system`, its receiver clock bias relative to the reference
  (meters), i.e. the inter-system bias. Reference-independent: the difference of
  two entries is the offset between those two systems. Empty for a single-system
  solution. In a GGTO-aided solution the Galileo entry is the broadcast offset
  `c · Δt_systems`.
"""
@kwdef struct PVTSolution
    position::ECEF = ECEF(0, 0, 0)
    velocity::ECEF = ECEF(0, 0, 0)
    time_correction::Float64 = 0
    time::Union{TAIEpoch{Float64},Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    sats::Dict{Tuple{Symbol,Int},SatInfo} = Dict{Tuple{Symbol,Int},SatInfo}()
    reference_system::Union{Symbol,Nothing} = nothing
    inter_system_biases::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
end

"""
    get_num_used_sats(pvt_solution::PVTSolution) -> Int

Return the number of satellites used in the PVT solution.
"""
function get_num_used_sats(pvt_solution::PVTSolution)
    length(pvt_solution.sats)
end

"""
    get_sat_info(pvt_solution::PVTSolution, system::Symbol, prn::Integer) -> Union{SatInfo,Nothing}

Return the [`SatInfo`](@ref) (position, transmit time and post-fit residual) of the
satellite with the given `prn` in GNSS time `system` (`:GPS`, `:Galileo`), or
`nothing` if that satellite was not used in the fix. The system tag is required
because the same PRN can belong to different constellations; see the `sats` field
of [`PVTSolution`](@ref).
"""
function get_sat_info(pvt_solution::PVTSolution, system::Symbol, prn::Integer)
    get(pvt_solution.sats, (system, Int(prn)), nothing)
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
    get_sat_enu(enu_from_ecef::ENUfromECEF, sat_pos_ecef::ECEF) -> Spherical

Convert satellite position to East-North-Up (ENU) spherical coordinates (azimuth and
elevation) relative to the user position.

# Arguments
- `user_pos_ecef::ECEF`: User position in ECEF coordinates
- `enu_from_ecef::ENUfromECEF`: A precomputed `ENUfromECEF(user_pos_ecef, wgs84)`
  transform. Pass this form to reuse one transform across all satellites that
  share a user position (one geodetic conversion per epoch instead of per satellite).
- `sat_pos_ecef::ECEF`: Satellite position in ECEF coordinates

# Returns
Spherical coordinates containing azimuth and elevation of the satellite as seen from
the user position.
"""
function get_sat_enu(enu_from_ecef::ENUfromECEF, sat_pos_ecef::ECEF)
    SphericalFromCartesian()(enu_from_ecef(sat_pos_ecef))
end

get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF) =
    get_sat_enu(ENUfromECEF(user_pos_ecef, wgs84), sat_pos_ecef)

"""
    get_time_system(decoder) -> Symbol

Return the GNSS time system a decoder's measurements are referenced to
(`:GPS` or `:Galileo`). Satellites sharing a time system share one receiver
clock bias in the PVT estimation.
"""
get_time_system(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData}) = :GPS
get_time_system(::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData}) = :Galileo

"""
    decide_clock_bias_layout(states, systems, unique_systems, system_indices, times)

Decide how many receiver clock biases to estimate and which one applies to each
satellite. Returns `(clock_bias_indices, range_offsets, num_clock_biases)` where `clock_bias_indices[j]` is the clock
column (1…`num_clock_biases`) of satellite `j`, `range_offsets[j]` is a known per-satellite range
correction (meters), and `num_clock_biases` is the number of clock biases. Returns `nothing`
when the constellation cannot be solved.

With `num_systems` distinct time systems the independent inter-system-bias solve
needs `3 + num_systems` satellites. When that fails but GPS is present and at
least one Galileo satellite has decoded the GGTO, the Galileo clock bias is
collapsed onto GPS via the broadcast offset (carried in `range_offsets`), reducing the
unknowns to `3 + (num_systems − 1)`. The GGTO is a constellation-wide parameter,
so a single Galileo satellite carrying it suffices to convert every Galileo
measurement.
"""
function decide_clock_bias_layout(states, systems, unique_systems, system_indices, times)
    num_sats = length(states)
    num_systems = length(unique_systems)
    clock_bias_indices = [system_indices[sys] for sys in systems]
    range_offsets = zeros(num_sats)
    num_sats < 3 + num_systems || return clock_bias_indices, range_offsets, num_systems

    # Under-determined: try to collapse Galileo onto GPS using the broadcast GGTO.
    # The GGTO is the same constellation-wide offset whichever Galileo satellite
    # reports it, so one decoded copy is enough for every Galileo measurement.
    ggto_idx = findfirst(
        j -> systems[j] == :Galileo && ggto_available(states[j].decoder),
        1:num_sats,
    )
    (
        haskey(system_indices, :GPS) &&
        !isnothing(ggto_idx) &&
        num_sats >= 3 + (num_systems - 1)
    ) || return nothing

    # The Galileo satellites share GPS's clock column and carry the known
    # inter-system offset (the shared GGTO evaluated at each transmit time) as a
    # per-satellite range correction.
    ggto_decoder = states[ggto_idx].decoder
    for j in 1:num_sats
        systems[j] == :Galileo || continue
        range_offsets[j] = SPEEDOFLIGHT * calc_ggto_offset(ggto_decoder, times[j])
    end
    merged_systems = filter(!=(:Galileo), unique_systems)
    merged_indices = Dict(sys => i for (i, sys) in enumerate(merged_systems))
    merged_indices[:Galileo] = merged_indices[:GPS]
    return [merged_indices[sys] for sys in systems], range_offsets, length(merged_systems)
end

"""
    predict_atmospheric_delays(ξ, states, sat_positions, correction,
                               reference_time, enable_tropospheric_correction) -> Vector{Float64}

Per-satellite ionospheric + tropospheric delay (metres), to be subtracted from the
pseudoranges. The user position is the first three elements of the least-squares
state vector `ξ = [x, y, z, tc₁, …]` (ECEF, metres); the remaining clock-bias
components are ignored. `correction` is the constellation-wide ionospheric model
from [`select_ionospheric_correction`](@ref) (`nothing` skips the ionosphere); the
troposphere uses the blind Saastamoinen model (see [`tropospheric_delay`](@ref))
unless `enable_tropospheric_correction` is `false`.

A single corrected solve is enough: the delays depend on position only through the
satellite elevation/azimuth (and, for the troposphere, the user height), and
∂delay/∂position is negligible over the metre-level position uncertainty (a 15 m
shift moves the elevation by ~1e-5°), so delays predicted at a nearby position are
accurate to well under a millimetre — no iterate-to-convergence needed. The user
geodetic coordinates and the ENU transform depend only on `ξ`, so they are
computed once and reused across satellites.
"""
function predict_atmospheric_delays(
    ξ,
    states,
    sat_positions,
    correction,
    reference_time,
    enable_tropospheric_correction,
)
    user_pos = ECEF(ξ[1], ξ[2], ξ[3])
    user_lla = LLAfromECEF(wgs84)(user_pos)
    enu_from_ecef = ENUfromECEF(user_pos, wgs84)
    map(states, sat_positions) do state, sat_pos
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
            enable_tropospheric_correction ? tropospheric_delay(elevation, user_lla) : 0.0
        iono + tropo
    end
end

"""
    calc_pvt(states::AbstractVector{<:SatelliteState},
             prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS satellite measurements.

Satellites from different GNSS may be combined. Because each constellation
references its broadcasts to its own system time, one receiver clock bias is
estimated per GNSS time system, so the state vector is
`[x, y, z, tc₁, …, tc_M]` for `M` distinct systems. Position and time are found
by least squares; velocity and clock drift are solved from carrier Doppler.

A solution requires `n ≥ 3 + M` healthy satellites (each system needs at least
one satellite, and a system contributing a single satellite spends it entirely
on that system's clock bias). When that condition fails but both GPS and Galileo
are present and the Galileo message carries the GGTO (Galileo–GPS Time Offset),
the Galileo clock bias is collapsed onto GPS using the broadcast offset, which
makes a 4-satellite GPS+Galileo fix possible. Estimating an independent bias is
preferred whenever the geometry allows it, since it avoids the GGTO's own
broadcast error.

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
satellite information. Returns `prev_pvt` if too few healthy satellites are
available to solve the constellation (including the GGTO fallback) or if the
GDOP is negative.

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
    healthy_indices = findall(x -> is_sat_healthy(x.decoder), states)
    length(healthy_indices) < 4 && return prev_pvt
    healthy_states = view(states, healthy_indices)
    num_sats = length(healthy_states)

    healthy_sat_keys =
        map(state -> (get_time_system(state.decoder), state.decoder.prn), healthy_states)
    times = map(calc_corrected_time, healthy_states)

    # Group satellites by GNSS time system (ordered by first appearance, so the
    # column layout is deterministic), then decide the clock-bias layout: one
    # column per system (independent inter-system bias estimation), or a
    # GGTO-aided collapse when the independent solve would be under-determined.
    systems = map(state -> get_time_system(state.decoder), healthy_states)
    unique_systems = unique(systems)
    system_indices = Dict(sys => i for (i, sys) in enumerate(unique_systems))
    clock_bias_layout =
        decide_clock_bias_layout(healthy_states, systems, unique_systems, system_indices, times)
    isnothing(clock_bias_layout) && return prev_pvt
    clock_bias_indices, range_offsets, num_clock_biases = clock_bias_layout

    # Propagating the ephemerides is only worthwhile once the constellation is
    # known to be solvable (otherwise we'd return prev_pvt above).
    sat_positions_and_velocities = map(
        (state, time) -> calc_satellite_position_and_velocity(state.decoder, time),
        healthy_states,
        times,
    )
    sat_positions = map(get_sat_position, sat_positions_and_velocities)
    # Built once here and reused by every least-squares solve below and by calc_H,
    # instead of being rebuilt inside each user_position call.
    sat_positions_mat = reduce(hcat, sat_positions)

    # Primary system — its clock bias, reference time, week and start epoch
    # define the reported time. A GGTO collapse (signalled by fewer clock biases
    # than systems) is anchored on GPS, so GPS must be primary there; otherwise
    # pick the system with the most satellites (best-conditioned reported time),
    # breaking ties by first appearance.
    primary_system =
        num_clock_biases < length(unique_systems) ? :GPS :
        unique_systems[argmax([count(==(sys), systems) for sys in unique_systems])]
    primary_clock_index = clock_bias_indices[findfirst(==(primary_system), systems)]

    # The common reference cancels out of the reported time (the primary clock
    # bias absorbs it), so any latest-transmit-time reference works.
    pseudo_ranges, reference_time = calc_pseudo_ranges(times)

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias.
    prev_abs_bias(sys) =
        sys == prev_pvt.reference_system ? prev_pvt.time_correction :
        prev_pvt.time_correction + get(prev_pvt.inter_system_biases, sys, 0.0)
    prev_ξ = zeros(3 + num_clock_biases)
    prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
    for j in 1:num_sats
        prev_ξ[3+clock_bias_indices[j]] = prev_abs_bias(systems[j])
    end

    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model is chosen for the whole solve from the
    # coefficients decoded across the constellation (NTCM-G if Galileo coefficients
    # are available, else Klobuchar if GPS α/β are available, else none; see
    # `select_ionospheric_correction`). The prediction is the top-level
    # `predict_atmospheric_delays`, a function barrier that specialises on the
    # concrete type of the `Union`-typed `ionospheric_correction`.
    ionospheric_correction =
        enable_ionospheric_correction ? select_ionospheric_correction(healthy_states) :
        nothing
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, resid_uncorrected = user_position(
            sat_positions_mat, pseudo_ranges, clock_bias_indices, num_clock_biases,
            range_offsets, prev_ξ)
        atmospheric_delays = predict_atmospheric_delays(
            ξ_uncorrected, healthy_states, sat_positions, ionospheric_correction,
            reference_time, enable_tropospheric_correction)
        any(!iszero, atmospheric_delays) ?
        user_position(
            sat_positions_mat, pseudo_ranges .- atmospheric_delays, clock_bias_indices,
            num_clock_biases, range_offsets, ξ_uncorrected) :
        (ξ_uncorrected, resid_uncorrected)
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction.
        user_position(
            sat_positions_mat,
            pseudo_ranges .- predict_atmospheric_delays(
                prev_ξ, healthy_states, sat_positions, ionospheric_correction,
                reference_time, enable_tropospheric_correction),
            clock_bias_indices,
            num_clock_biases,
            range_offsets,
            prev_ξ,
        )
    end
    H = calc_H(sat_positions_mat, ξ, clock_bias_indices, num_clock_biases)
    # H feeds both the per-system DOP and the velocity/clock-drift solve, which
    # reuses its line-of-sight block and collapses its clock columns to a single
    # common receiver clock drift.
    user_velocity_and_clock_drift = calc_user_velocity_and_clock_drift(
        sat_positions_and_velocities, healthy_states, times, H)
    position = ECEF(ξ[1], ξ[2], ξ[3])
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    relative_clock_drift = user_velocity_and_clock_drift[4] / SPEEDOFLIGHT
    time_correction = ξ[3+primary_clock_index]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEEDOFLIGHT

    primary_decoder = healthy_states[findfirst(==(primary_system), systems)].decoder
    week = get_week(primary_decoder; approximate_year)
    start_time = get_system_start_time(primary_decoder)
    time = TAIEpoch(
        week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second,
        corrected_reference_time - floor(Int, corrected_reference_time),
    )

    sat_infos = SatInfo.(sat_positions, times, residuals)

    dop = calc_DOP(H, position, primary_clock_index)
    if dop.GDOP < 0
        return prev_pvt
    end

    # Inter-system biases relative to the reference (primary) system's clock, in
    # meters. The reference is omitted (its bias is `time_correction`); for a
    # GGTO-collapsed system this is the broadcast offset c·Δt_systems.
    inter_system_biases = Dict{Symbol,Float64}()
    for j in 1:num_sats
        systems[j] == primary_system && continue
        inter_system_biases[systems[j]] =
            ξ[3+clock_bias_indices[j]] + range_offsets[j] - time_correction
    end

    PVTSolution(
        position,
        velocity,
        time_correction,
        time,
        relative_clock_drift,
        dop,
        Dict(healthy_sat_keys .=> sat_infos),
        primary_system,
        inter_system_biases,
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
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData},
)
    TAIEpoch(1980, 1, 6, 0, 0, 19.0) # There were 19 leap seconds at 01/06/1999 compared to UTC
end

function get_system_start_time(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GalileoE1BData},
)
    TAIEpoch(1999, 8, 22, 0, 0, (32 - 13.0)) # There were 32 leap seconds at 08/22/1999 compared to UTC
end

"""
    get_week(decoder::GNSSDecoderState{<:GPSL1CAData}; approximate_year)

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
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData};
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
