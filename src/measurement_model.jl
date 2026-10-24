# The pseudorange and pseudorange-rate model: predicted pseudoranges, their design matrix
# and curvature, the Earth-rotation correction, the rate model, and the elevation and
# azimuth of a satellite.

# The satellite positions the model and its estimators take come in either of two layouts: a
# `3 × N` matrix (one column per satellite), or a vector of `SVector{3,Float64}` — what
# the solver itself keeps, since a satellite position is a fixed-size quantity and a
# vector of them is resized in place where a matrix cannot be.
satellite_count(satellite_positions::AbstractMatrix) = size(satellite_positions, 2)
satellite_count(satellite_positions::AbstractVector) = length(satellite_positions)
get_satellite_position(satellite_positions::AbstractMatrix, j) = SVector{3}(
    satellite_positions[1, j], satellite_positions[2, j], satellite_positions[3, j])
get_satellite_position(satellite_positions::AbstractVector, j) =
    SVector{3,Float64}(satellite_positions[j])

"""
    predict_pseudoranges!(pseudoranges, satellite_positions, ξ, bias_columns::BiasColumns)
        -> pseudoranges

The predicted pseudoranges ρ̂, written into `pseudoranges`: per satellite, the distance
between the satellite and the estimated user position, plus the receiver clock and
inter-frequency biases.

- `satellite_positions`: Satellite positions, a `3 × N` matrix or a vector of `SVector{3}`s.
- `ξ`: Estimated user position, one receiver clock bias per GNSS time system, and one
  receiver inter-frequency bias per frequency band beyond the reference band, i.e.
  `[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_inter_frequency_biases]`.
- `bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock column and
  inter-frequency-bias column of each satellite). Known range corrections (atmosphere)
  are not modelled here: [`calc_pvt`](@ref) folds them into the measured pseudoranges.

After the satellites, `pseudoranges` holds one modelled clock difference per
[`ClockConstraint`](@ref) of `bias_columns`, so it has
[`num_observations`](@ref)`(satellite_positions, bias_columns)` entries.
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

# The ECEF frame turns by `EARTH_ROTATION_RATE · Δt` while a signal is in flight for
# `Δt`; this rotation carries a vector given in the ECEF frame of transmission into the
# one of reception.
@inline function earth_rotation(Δt)
    α = EARTH_ROTATION_RATE * Δt
    @SMatrix [
        cos(α) sin(α) 0
        -sin(α) cos(α) 0
        0 0 1
    ]
end

# The rotation that carries the satellite at `sat_pos`, in the ECEF frame of its
# transmission, into the ECEF frame of reception at `user_pos`: the Earth's turn over the
# signal's travel time `|s − u|/c`. The pseudorange model (`predict_pseudoranges!`), its
# position partials (`calc_line_of_sight` for `calc_design_matrix!`), its curvature
# (`calc_pseudorange_curvature!`) and its rate (`pseudorange_rate_and_gradient`) all take
# the Earth-rotation correction from here, so they describe one and the same satellite.
@inline reception_frame_rotation(user_pos, sat_pos) =
    earth_rotation(norm(sat_pos - user_pos) / SPEED_OF_LIGHT)

# The satellite at `sat_pos` in the ECEF frame of reception at `user_pos`.
@inline satellite_position_at_reception(user_pos, sat_pos) =
    reception_frame_rotation(user_pos, sat_pos) * sat_pos

# The range-rate model: the exact time derivative of the pseudorange `predict_pseudoranges!`
# models, which the velocity solve fits (`estimate_velocity`) and an estimator of its own
# predicts with (`calc_pseudorange_rate`).

# The satellite's position and velocity in the ECEF frame of reception, as
# `predict_pseudoranges!` models the position: `s' = R(τ)·s`, rotated by the Earth's turn
# over the travel time. The rotation angle `ω_e·τ` changes with the travel time, so `s'`
# moves at `R(τ)·ṡ + ω_e·τ̇·(s'_y, −s'_x, 0)` — the rotated velocity plus the rate of the
# rotation itself (`dR/dα = R·J` with `J` the quarter turn about z, which commutes with
# `R`). `τ̇` is taken for a resting receiver, `û·v_sat/c`; its user-velocity share is linear
# and is carried separately (`travel_time_velocity_partial`).
@inline function satellite_state_at_reception(user_pos, sat_pos, sat_vel)
    to_sat = sat_pos - user_pos
    travel_time_rate = dot(to_sat, sat_vel) / (norm(to_sat) * SPEED_OF_LIGHT)
    # One rotation for both, the one `predict_pseudoranges!` turns the satellite by.
    rotation = reception_frame_rotation(user_pos, sat_pos)
    position = rotation * sat_pos
    velocity =
        rotation * sat_vel +
        EARTH_ROTATION_RATE * travel_time_rate * SVector(position[2], -position[1], 0.0)
    return position, velocity
end

# The user-velocity partial of the travel-time rate's share of the range rate: the
# Earth-rotation-rate term `ω_e·τ̇·(e · J·s')` (with `J·s' = (s'_y, −s'_x, 0)`) takes
# `−ω_e·(e · J·s')·û/c` from the `−û·v_user/c` in `τ̇ = û·(v_sat − v_user)/c`. It is ~2e-6
# relative to the line of sight, so it is evaluated on the unrotated geometry: with
# `e ≈ û` and `s' ≈ s`, and `s·J·s = 0`, `e · J·s' ≈ −u·J·s/d`. The neglected rotation is
# ~5e-6 relative of this term, ~1e-11 of the partial.
@inline function travel_time_velocity_partial(user_pos, sat_pos)
    to_sat = sat_pos - user_pos
    # `−u·J·s`; the `1/d` of `e · J·s'` and of `û` combine with `to_sat` into `1/d²`.
    sagnac = user_pos[2] * sat_pos[1] - user_pos[1] * sat_pos[2]
    scale = -EARTH_ROTATION_RATE * sagnac / (SPEED_OF_LIGHT * dot(to_sat, to_sat))
    return scale * to_sat
end

"""
    calc_pseudorange_rate(user_position, user_velocity, receiver_clock_drift,
        satellite_position, satellite_velocity, satellite_clock_drift) -> Float64

The modelled pseudorange rate (m/s) of one satellite: the exact time derivative, per second
of receiver clock, of the pseudorange [`predict_pseudoranges!`](@ref) models — the satellite
rotated by the Earth's turn over the travel time into the ECEF frame of reception, plus the
receiver clock bias, minus the satellite clock offset. It is the model the velocity solve
([`estimate_velocity`](@ref)) fits, and its rate residuals are the measured rate minus this
at the solved state; an estimator that predicts the rate from its own state — a
vector-tracking navigation filter, say — calls it to predict exactly what the solve fits.

Positive while the range grows (the satellite recedes), the sense of the velocity solve
and of RTKLIB's `resdop`; the Doppler-signed rate a tracking loop works in,
`λ · carrier_doppler`, is its negative.

`user_position`, `user_velocity`, `satellite_position` and `satellite_velocity` are ECEF
`SVector{3}`s in m and m/s, the satellite's taken at transmission. `receiver_clock_drift` is
the receiver clock drift in m/s, as `PVTSolution`'s `clock_drift` reports it, and
`satellite_clock_drift` the satellite clock drift in s/s, as
[`calc_satellite_clock_drift`](@ref) returns it: the receiver's clock is estimated in
range units, the satellite's is broadcast in time units.

Beyond the line-of-sight projection of the relative velocity it carries two terms of
order `v/c`, each up to about a millimetre per second:

- the rate of the Earth-rotation angle, `ω_e·τ̇` with `τ̇ = û·(v_sat − v_user)/c`
  ([`EARTH_ROTATION_RATE`](@ref)), and
- the light-time factor `k = 1 − ρ̇/c`: one second of receiver clock spans that many
  seconds of the satellite's motion and clock, `ρ̇` being the rate of the range including
  the receiver clock drift.

Allocation-free and `@inline`.
"""
@inline calc_pseudorange_rate(user_position, user_velocity, receiver_clock_drift,
    satellite_position, satellite_velocity, satellite_clock_drift) =
    first(pseudorange_rate_and_gradient(user_position, user_velocity, receiver_clock_drift,
        satellite_position, satellite_velocity, satellite_clock_drift))

# `calc_pseudorange_rate` and its gradient in `[v; ċ]`, the velocity solve's design row.
# With the receiver velocity `v` and clock drift `ċ` (m/s) the pseudorange rate is
#
#     ρ̇ = k·s + r,   s = rate + q · v,   r = ċ − e · v,   k = 1 − (s + r)/c,
#
# where
# - `e` is the unit receiver→satellite vector to the Earth-rotation-corrected satellite,
#   as `calc_line_of_sight` gives it;
# - `rate` is the satellite's share of the range rate for a resting receiver (m/s),
#   `e · d(s')/dt − c·ḃ_sat` — its rotated velocity, the rate of the rotation itself, and
#   its clock drift;
# - `q` is how that share changes with the receiver velocity, through the travel-time
#   rate the rotation rate is scaled by (~2e-6 against the line of sight);
# - `k` is the light-time factor.
#
# `s` and `r` are affine in `[v; ċ]`, so `ρ̇` is quadratic, with the gradient
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

# The unit line of sight from the user to the Earth-rotation-corrected satellite — the
# textbook receiver→satellite direction. The pseudorange's position partial is its
# *negative*, which is what the design-matrix rows of `calc_design_matrix!` carry. `ξ` holds
# the user position in its first three entries.
function calc_line_of_sight(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    to_sat = satellite_position_at_reception(rₙ, sat_pos) - rₙ
    return to_sat / norm(to_sat)
end

"""
    calc_design_matrix!(H, satellite_positions, ξ, bias_columns::BiasColumns) -> H

The least-squares design matrix at the state `ξ`, written into `H`: the Jacobian of
[`predict_pseudoranges!`](@ref), whose arguments it shares.

The matrix has [`num_parameters`](@ref)`(bias_columns)` columns: three for the position
partials, one per GNSS time system, and one per frequency band beyond the reference band.
Each row carries a `1.0` in its system's clock column (`3 + clock_columns[j]`) and,
unless the satellite is on the reference band, a `1.0` in its band's inter-frequency-bias
column (`3 + num_clock_biases + inter_frequency_bias_columns[j]`); the remaining bias
columns are zero. Each [`ClockConstraint`](@ref) adds a row after the satellites, `+1` in
its `column` and `−1` in its `target_column`.
"""
function calc_design_matrix!(H, satellite_positions, ξ, bias_columns::BiasColumns)
    fill!(H, 0.0)
    for j in 1:satellite_count(satellite_positions)
        sat_pos = get_satellite_position(satellite_positions, j)
        # The position partial ∂ρ/∂r is the negative of the receiver→satellite
        # line of sight.
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

Allocating form of [`calc_design_matrix!`](@ref PositionVelocityTime.calc_design_matrix!):
the least-squares design matrix at the state `ξ`, one row per satellite and per
[`ClockConstraint`](@ref), with
[`num_parameters`](@ref PositionVelocityTime.num_parameters)`(bias_columns)` columns.
"""
calc_design_matrix(satellite_positions, ξ, bias_columns::BiasColumns) = calc_design_matrix!(
    Matrix{Float64}(undef, num_observations(satellite_positions, bias_columns),
        num_parameters(bias_columns)),
    satellite_positions, ξ, bias_columns)

"""
    calc_pseudorange_curvature!(curvature, satellite_positions, ξ, v) -> curvature

The curvature of the predicted pseudoranges along the step `v`: the second directional
derivative of `predict_pseudoranges!`, written into `curvature`, used by the
Levenberg–Marquardt geodesic acceleration.

For each satellite j the residual is `r_j(ξ) = ‖s_j' - r_n‖ + tc - ρ_j`,
where s_j' is the Earth-rotation-corrected satellite position. Treating
s_j' as constant w.r.t. ξ (the rotation depends on ξ via travel time,
but ω_e/c ≈ 2e-13 makes that contribution negligible here), the Hessian
is block-diagonal: the position block is `(I - û û^T) / d_j`, the time
block is zero. So `v^T H_j v = (‖v_r‖² - (û · v_r)²) / d_j`. The clock-constraint rows
after the satellites are linear, so their entries are zero.
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

Elevation and azimuth (radians) of the satellite at `satellite_position` (ECEF) as seen
from the user position, in the local East-North-Up frame there. Azimuth is measured
clockwise from North and wrapped to `[0, 2π)`, the GNSS convention.

Pass a precomputed `enu_from_ecef = ENUfromECEF(user_position, wgs84)` to reuse one
transform across all satellites that share a user position (one geodetic conversion per
epoch instead of per satellite).
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
