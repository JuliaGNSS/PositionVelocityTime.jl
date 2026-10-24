# The pseudorange and pseudorange-rate model: predicted pseudoranges, their design matrix
# and curvature, the Earth-rotation correction, the rate model, and the elevation and
# azimuth of a satellite.

# Satellite positions come as a `3 × N` matrix or as a vector of `SVector{3,Float64}`
# (what the solver keeps, as it resizes in place).
satellite_count(satellite_positions::AbstractMatrix) = size(satellite_positions, 2)
satellite_count(satellite_positions::AbstractVector) = length(satellite_positions)
get_satellite_position(satellite_positions::AbstractMatrix, j) = SVector{3}(
    satellite_positions[1, j], satellite_positions[2, j], satellite_positions[3, j])
get_satellite_position(satellite_positions::AbstractVector, j) =
    SVector{3,Float64}(satellite_positions[j])

"""
    predict_pseudoranges!(pseudoranges, satellite_positions, ξ, bias_columns::BiasColumns)
        -> pseudoranges

The predicted pseudoranges ρ̂ at the state `ξ`, written into `pseudoranges`: per
satellite, the range to the Earth-rotation-corrected satellite plus its clock and
inter-frequency biases, then one modelled clock difference per [`ClockConstraint`](@ref)
([`num_observations`](@ref) entries in all).

- `satellite_positions`: a `3 × N` matrix or a vector of `SVector{3}`s.
- `ξ`: `[x, y, z, tc₁, …, ifb₁, …]`, laid out by `bias_columns` (see
  [`BiasColumns`](@ref)).

Atmospheric delays are not modelled here; [`calc_pvt`](@ref) removes them from the
measured pseudoranges.
"""
function predict_pseudoranges!(
    pseudoranges,
    satellite_positions,
    ξ,
    bias_columns::BiasColumns,
)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    num_sats = satellite_count(satellite_positions)
    for j in 1:num_sats
        sat_pos = get_satellite_position(satellite_positions, j)
        ifb_column = bias_columns.inter_frequency_bias_columns[j]
        ifb = ifb_column == 0 ? 0.0 : ξ[3+bias_columns.num_clock_biases+ifb_column]
        pseudoranges[j] = norm(satellite_position_at_reception(rₙ, sat_pos) - rₙ) +
               ξ[3+bias_columns.clock_columns[j]] + ifb
    end
    for (i, constraint) in enumerate(bias_columns.constraints)
        pseudoranges[num_sats+i] = ξ[3+constraint.column] - ξ[3+constraint.target_column]
    end
    return pseudoranges
end

"""
    num_observations(satellite_positions, bias_columns::BiasColumns) -> Int

Rows of the least-squares problem: one per satellite and one per
[`ClockConstraint`](@ref).
"""
num_observations(satellite_positions, bias_columns::BiasColumns) =
    satellite_count(satellite_positions) + length(bias_columns.constraints)

# The rotation carrying a vector from the ECEF frame of transmission into that of
# reception, `Δt` later.
@inline function earth_rotation(Δt)
    α = EARTH_ROTATION_RATE * Δt
    @SMatrix [
        cos(α) sin(α) 0
        -sin(α) cos(α) 0
        0 0 1
    ]
end

# The Earth's turn over the travel time `|s − u|/c`. Every part of the model — the
# pseudorange, its partials, its curvature and its rate — takes the rotation from here,
# so they describe one and the same satellite.
@inline reception_frame_rotation(user_pos, sat_pos) =
    earth_rotation(norm(sat_pos - user_pos) / SPEED_OF_LIGHT)

# The satellite at `sat_pos` in the ECEF frame of reception at `user_pos`.
@inline satellite_position_at_reception(user_pos, sat_pos) =
    reception_frame_rotation(user_pos, sat_pos) * sat_pos

# The satellite's position and velocity in the frame of reception: `s' = R(τ)·s`, moving
# at `R(τ)·ṡ + ω_e·τ̇·(s'_y, −s'_x, 0)`, the second term being the rate of the rotation
# itself. `τ̇` is taken for a resting receiver, `û·v_sat/c`; the user-velocity share is
# `travel_time_velocity_partial`.
@inline function satellite_state_at_reception(user_pos, sat_pos, sat_vel)
    to_sat = sat_pos - user_pos
    travel_time_rate = dot(to_sat, sat_vel) / (norm(to_sat) * SPEED_OF_LIGHT)
    rotation = reception_frame_rotation(user_pos, sat_pos)
    position = rotation * sat_pos
    velocity =
        rotation * sat_vel +
        EARTH_ROTATION_RATE * travel_time_rate * SVector(position[2], -position[1], 0.0)
    return position, velocity
end

# The user-velocity partial of the rotation-rate term: `−ω_e·(e · J·s')·û/c`, from the
# `−û·v_user/c` in `τ̇`. It is ~2e-6 of the line of sight, so it is taken on the unrotated
# geometry (`e ≈ û`, `s' ≈ s`, so `e · J·s' ≈ −u·J·s/d`), an error of ~1e-11 of the
# partial.
@inline function travel_time_velocity_partial(user_pos, sat_pos)
    to_sat = sat_pos - user_pos
    # `−u·J·s`; the two `1/d` combine with `to_sat` into `1/d²`.
    sagnac = user_pos[2] * sat_pos[1] - user_pos[1] * sat_pos[2]
    scale = -EARTH_ROTATION_RATE * sagnac / (SPEED_OF_LIGHT * dot(to_sat, to_sat))
    return scale * to_sat
end

"""
    calc_pseudorange_rate(user_position, user_velocity, receiver_clock_drift,
        satellite_position, satellite_velocity, satellite_clock_drift) -> Float64

The modelled pseudorange rate (m/s) of one satellite: the exact time derivative, per
second of receiver clock, of what [`predict_pseudoranges!`](@ref) models. It is the model
[`estimate_velocity`](@ref) fits, for an estimator that predicts the rate from its own
state. Positive while the satellite recedes (see [`SatelliteInfo`](@ref)).

Positions and velocities are ECEF `SVector{3}`s in m and m/s, the satellite's at
transmission. `receiver_clock_drift` is in m/s, as `PVTSolution.clock_drift`;
`satellite_clock_drift` in s/s, as [`calc_satellite_clock_drift`](@ref) returns it.

Beyond the line-of-sight relative velocity it carries two terms of order `v/c`, each up
to ~1 mm/s: the rate of the Earth-rotation angle, `ω_e·τ̇` with `τ̇ = û·(v_sat − v_user)/c`,
and the light-time factor `k = 1 − ρ̇/c` (a second of receiver clock spans that many of
the satellite's).

Allocation-free and `@inline`.
"""
@inline calc_pseudorange_rate(user_position, user_velocity, receiver_clock_drift,
    satellite_position, satellite_velocity, satellite_clock_drift) =
    first(pseudorange_rate_and_gradient(user_position, user_velocity, receiver_clock_drift,
        satellite_position, satellite_velocity, satellite_clock_drift))

# `calc_pseudorange_rate` and its gradient in `[v; ċ]`, the velocity solve's design row:
#
#     ρ̇ = k·s + r,   s = rate + q · v,   r = ċ − e · v,   k = 1 − (s + r)/c,
#
# with `e` the line of sight to the rotated satellite, `rate` the satellite's share for a
# resting receiver (rotated velocity, rotation rate and clock drift), `q` its change with
# the receiver velocity (`travel_time_velocity_partial`), and `k` the light-time factor.
# The gradient is
#
#     ∂ρ̇/∂v = k·q − e − s·(q − e)/c,   ∂ρ̇/∂ċ = 1 − s/c.
@inline function pseudorange_rate_and_gradient(user_pos, user_vel, receiver_clock_drift,
    sat_pos, sat_vel, sat_clock_drift)
    rotated_sat_pos, rotated_sat_vel =
        satellite_state_at_reception(user_pos, sat_pos, sat_vel)
    to_rotated_sat = rotated_sat_pos - user_pos
    e = to_rotated_sat / norm(to_rotated_sat)
    q = travel_time_velocity_partial(user_pos, sat_pos)
    s = dot(e, rotated_sat_vel) - sat_clock_drift * SPEED_OF_LIGHT + dot(q, user_vel)
    r = receiver_clock_drift - dot(e, user_vel)
    k = 1 - (s + r) / SPEED_OF_LIGHT
    ∂v = k * q - e - (s / SPEED_OF_LIGHT) * (q - e)
    return k * s + r, SVector(∂v[1], ∂v[2], ∂v[3], 1 - s / SPEED_OF_LIGHT)
end

# The unit line of sight from the user to the rotated satellite; the design matrix's
# position partials are its negative.
function calc_line_of_sight(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    to_sat = satellite_position_at_reception(rₙ, sat_pos) - rₙ
    return to_sat / norm(to_sat)
end

"""
    calc_design_matrix!(H, satellite_positions, ξ, bias_columns::BiasColumns) -> H

The Jacobian of [`predict_pseudoranges!`](@ref) at `ξ`, written into `H`, with
[`num_parameters`](@ref)`(bias_columns)` columns. A satellite row holds the negated line
of sight and a `1.0` in its clock column and, off the reference band, in its
inter-frequency-bias column. A [`ClockConstraint`](@ref) row holds `+1` in its `column`
and `−1` in its `target_column`.
"""
function calc_design_matrix!(H, satellite_positions, ξ, bias_columns::BiasColumns)
    fill!(H, 0.0)
    for j in 1:satellite_count(satellite_positions)
        sat_pos = get_satellite_position(satellite_positions, j)
        e = calc_line_of_sight(sat_pos, ξ)
        H[j, 1] = -e[1]
        H[j, 2] = -e[2]
        H[j, 3] = -e[3]
        H[j, 3+bias_columns.clock_columns[j]] = 1.0
        ifb_column = bias_columns.inter_frequency_bias_columns[j]
        if ifb_column != 0
            H[j, 3+bias_columns.num_clock_biases+ifb_column] = 1.0
        end
    end
    num_sats = satellite_count(satellite_positions)
    for (i, constraint) in enumerate(bias_columns.constraints)
        H[num_sats+i, 3+constraint.column] = 1.0
        H[num_sats+i, 3+constraint.target_column] = -1.0
    end
    return H
end

"""
    calc_design_matrix(satellite_positions, ξ, bias_columns::BiasColumns) -> Matrix{Float64}

Allocating form of [`calc_design_matrix!`](@ref PositionVelocityTime.calc_design_matrix!).
"""
calc_design_matrix(satellite_positions, ξ, bias_columns::BiasColumns) = calc_design_matrix!(
    Matrix{Float64}(undef, num_observations(satellite_positions, bias_columns),
        num_parameters(bias_columns)),
    satellite_positions, ξ, bias_columns)

"""
    calc_pseudorange_curvature!(curvature, satellite_positions, ξ, v) -> curvature

The second directional derivative of `predict_pseudoranges!` along the step `v`, for the
Levenberg–Marquardt geodesic acceleration: `(‖v_r‖² − (û · v_r)²) / d_j` per satellite,
treating the rotated satellite as fixed (its dependence on `ξ` is scaled by
`ω_e/c ≈ 2e-13`). The clock-constraint rows are linear, so theirs is zero.
"""
function calc_pseudorange_curvature!(curvature, satellite_positions, ξ, v)
    fill!(curvature, 0.0)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    v_r = SVector{3}(v[1], v[2], v[3])
    v_r_sq = dot(v_r, v_r)
    for j in 1:satellite_count(satellite_positions)
        sat_pos = get_satellite_position(satellite_positions, j)
        u = satellite_position_at_reception(rₙ, sat_pos) - rₙ
        d = norm(u)
        û = u / d
        curvature[j] = (v_r_sq - dot(û, v_r)^2) / d
    end
    return curvature
end

"""
    calc_elevation_and_azimuth(user_position::ECEF, satellite_position)
        -> (elevation, azimuth)
    calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, satellite_position)
        -> (elevation, azimuth)

Elevation and azimuth (radians) of the satellite at `satellite_position` (ECEF), seen
from the user position in its local East-North-Up frame; the azimuth runs clockwise from
North in `[0, 2π)`. Pass `enu_from_ecef = ENUfromECEF(user_position, wgs84)` to share one
geodetic conversion across the satellites of an epoch.
"""
function calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, satellite_position)
    enu = SphericalFromCartesian()(enu_from_ecef(ECEF(satellite_position)))
    elevation = enu.ϕ
    # `SphericalFromCartesian` measures θ counter-clockwise from East (the ENU +x
    # axis); azimuth is measured clockwise from North, i.e. π/2 − θ.
    azimuth = mod2pi(π / 2 - enu.θ)
    return elevation, azimuth
end

calc_elevation_and_azimuth(user_position::ECEF, satellite_position) =
    calc_elevation_and_azimuth(ENUfromECEF(user_position, wgs84), satellite_position)
