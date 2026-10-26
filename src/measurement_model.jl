# The pseudorange and pseudorange-rate model: predicted pseudoranges, their design matrix
# and curvature, the Earth-rotation correction, the rate model, and the elevation and
# azimuth of a satellite.

# The satellite positions the model and its estimators take come in either of two layouts: a
# `3 × N` matrix (one column per satellite), or a vector of `SVector{3,Float64}` — what
# the solver itself keeps, since a satellite position is a fixed-size quantity and a
# vector of them is resized in place where a matrix cannot be.
satellite_count(sat_positions::AbstractMatrix) = size(sat_positions, 2)
satellite_count(sat_positions::AbstractVector) = length(sat_positions)
get_satellite_position(sat_positions::AbstractMatrix, j) = SVector{3}(
    sat_positions[1, j], sat_positions[2, j], sat_positions[3, j])
get_satellite_position(sat_positions::AbstractVector, j) =
    SVector{3,Float64}(sat_positions[j])

"""
    predict_pseudoranges!(pseudoranges, sat_positions, ξ, bias_columns::BiasColumns)
        -> pseudoranges

The predicted pseudoranges ρ̂ at the state `ξ`, written into `pseudoranges`: per
satellite, the range to the Earth-rotation-corrected satellite plus its clock and
inter-frequency biases.

- `sat_positions`: a `3 × N` matrix or a vector of `SVector{3}`s.
- `ξ`: `[x, y, z, tc₁, …, ifb₁, …]`, laid out by `bias_columns` (see
  [`BiasColumns`](@ref)).

Atmospheric delays and the broadcast offsets onto the hub are not modelled here;
[`calc_pvt`](@ref) removes them from the measured pseudoranges.
"""
function predict_pseudoranges!(
    pseudoranges,
    sat_positions,
    ξ,
    bias_columns::BiasColumns,
)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    for j in 1:satellite_count(sat_positions)
        sat_pos = get_satellite_position(sat_positions, j)
        ifb_column = bias_columns.ifb_columns[j]
        ifb = ifb_column == 0 ? 0.0 : ξ[3+bias_columns.num_clock_biases+ifb_column]
        pseudoranges[j] = norm(satellite_position_at_reception(rₙ, sat_pos) - rₙ) +
               ξ[3+bias_columns.clock_columns[j]] + ifb
    end
    return pseudoranges
end

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
# signal's travel time `|s − u|/c`. The pseudorange model (`predict_pseudoranges!`), its position
# partials (`calc_line_of_sight` for `calc_design_matrix!`, `calc_pseudorange_curvature!`) and its rate
# (`pseudorange_rate_and_gradient`) all take the Earth-rotation correction from here, so
# they describe one and the same satellite.
@inline reception_frame_rotation(user_pos, sat_pos) =
    earth_rotation(norm(sat_pos - user_pos) / SPEED_OF_LIGHT)

# The satellite at `sat_pos` in the ECEF frame of reception at `user_pos`.
@inline satellite_position_at_reception(user_pos, sat_pos) =
    reception_frame_rotation(user_pos, sat_pos) * sat_pos

# The range-rate model: the exact time derivative of the pseudorange `predict_pseudoranges!` models,
# which the velocity solve fits (`estimate_velocity`) and an estimator of
# its own predicts with (`calc_pseudorange_rate`).

# The satellite's position and velocity in the ECEF frame of reception, as `predict_pseudoranges!`
# models the position: `s' = R(τ)·s`, rotated by the Earth's turn over the travel time.
# The rotation angle `ω_e·τ` changes with the travel time, so `s'` moves at
# `R(τ)·ṡ + ω_e·τ̇·(s'_y, −s'_x, 0)` — the rotated velocity plus the rate of the rotation
# itself (`dR/dα = R·J` with `J` the quarter turn about z, which commutes with `R`). `τ̇`
# is taken for a resting receiver, `û·v_sat/c`; its user-velocity share is linear and is
# carried separately (`travel_time_velocity_partial`).
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
    calc_pseudorange_rate(user_pos, user_vel, user_clock_drift, sat_pos, sat_vel,
        sat_clock_drift) -> Float64

The modelled pseudorange rate (m/s) of one satellite: the exact time derivative, per
second of receiver clock, of the pseudorange [`predict_pseudoranges!`](@ref) models — the satellite
rotated by the Earth's turn over the travel time into the ECEF frame of reception, plus
the receiver clock bias, minus the satellite clock offset. It is the model the velocity
solve ([`estimate_velocity`](@ref)) fits, and its rate residuals are the
measured rate minus this at the solved state; an estimator that predicts the rate from
its own state — a vector-tracking navigation filter, say — calls it to predict exactly
what the solve fits.

Positive while the range grows (the satellite recedes), the sense of the velocity solve
and of RTKLIB's `resdop`; the Doppler-signed rate a tracking loop works in,
`λ · carrier_doppler`, is its negative.

`user_pos`, `user_vel`, `sat_pos` and `sat_vel` are ECEF `SVector{3}`s in m and m/s, the
satellite's taken at transmission. `user_clock_drift` is the receiver clock drift in m/s,
but `sat_clock_drift` the satellite clock drift as a rate in s/s, as
[`calc_satellite_clock_drift`](@ref) returns it.

Beyond the line-of-sight projection of the relative velocity it carries two terms of
order `v/c`, each up to about a millimetre per second:

- the rate of the Earth-rotation angle, `ω_e·τ̇` with `τ̇ = û·(v_sat − v_user)/c`
  ([`EARTH_ROTATION_RATE`](@ref)), and
- the light-time factor `k = 1 − ρ̇/c`: one second of receiver clock spans that many
  seconds of the satellite's motion and clock, `ρ̇` being the rate of the range including
  the receiver clock drift.

Allocation-free and `@inline`.
"""
@inline calc_pseudorange_rate(user_pos, user_vel, user_clock_drift, sat_pos, sat_vel,
    sat_clock_drift) = first(pseudorange_rate_and_gradient(user_pos, user_vel,
    user_clock_drift, sat_pos, sat_vel, sat_clock_drift))

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
@inline function pseudorange_rate_and_gradient(user_pos, user_vel, user_clock_drift,
    sat_pos, sat_vel, sat_clock_drift)
    rotated_sat_pos, rotated_sat_vel =
        satellite_state_at_reception(user_pos, sat_pos, sat_vel)
    to_rotated_sat = rotated_sat_pos - user_pos
    e = to_rotated_sat / norm(to_rotated_sat)
    q = travel_time_velocity_partial(user_pos, sat_pos)
    s = dot(e, rotated_sat_vel) - sat_clock_drift * SPEED_OF_LIGHT + dot(q, user_vel)
    r = user_clock_drift - dot(e, user_vel)
    k = 1 - (s + r) / SPEED_OF_LIGHT
    ∂v = k * q - e - (s / SPEED_OF_LIGHT) * (q - e)
    return k * s + r, SVector(∂v[1], ∂v[2], ∂v[3], 1 - s / SPEED_OF_LIGHT)
end

# The unit line of sight from the user to the Earth-rotation-corrected satellite — the
# textbook receiver→satellite direction. The pseudorange's position partial is its
# *negative*, which is what the design-matrix rows of `calc_design_matrix!` carry. `ξ` holds the user
# position in its first three entries.
function calc_line_of_sight(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    to_sat = satellite_position_at_reception(rₙ, sat_pos) - rₙ
    return to_sat / norm(to_sat)
end

"""
    calc_design_matrix!(H, sat_positions, ξ, bias_columns::BiasColumns) -> H

The Jacobian of [`predict_pseudoranges!`](@ref) at `ξ`, written into `H`, with
[`num_parameters`](@ref)`(bias_columns)` columns: three for the position partials, one
per GNSS time system, and one per frequency band beyond the reference band. Each row
holds the negated line of sight and a `1.0` in its clock column and, off the reference
band, in its inter-frequency-bias column.
"""
function calc_design_matrix!(H, sat_positions, ξ, bias_columns::BiasColumns)
    fill!(H, 0.0)
    for j in 1:satellite_count(sat_positions)
        sat_pos = get_satellite_position(sat_positions, j)
        e = calc_line_of_sight(sat_pos, ξ)
        H[j, 1] = -e[1]
        H[j, 2] = -e[2]
        H[j, 3] = -e[3]
        H[j, 3+bias_columns.clock_columns[j]] = 1.0
        ifb_column = bias_columns.ifb_columns[j]
        if ifb_column != 0
            H[j, 3+bias_columns.num_clock_biases+ifb_column] = 1.0
        end
    end
    return H
end

"""
    calc_design_matrix(sat_positions, ξ, bias_columns::BiasColumns) -> Matrix{Float64}

Allocating form of [`calc_design_matrix!`](@ref PositionVelocityTime.calc_design_matrix!).
"""
calc_design_matrix(sat_positions, ξ, bias_columns::BiasColumns) = calc_design_matrix!(
    Matrix{Float64}(undef, satellite_count(sat_positions),
        num_parameters(bias_columns)),
    sat_positions, ξ, bias_columns)

"""
    calc_pseudorange_curvature!(curvature, sat_positions, ξ, v) -> curvature

The second directional derivative of `predict_pseudoranges!` along the step `v`, for the
Levenberg–Marquardt geodesic acceleration: `(‖v_r‖² − (û · v_r)²) / d_j` per satellite,
treating the rotated satellite as fixed (its dependence on `ξ` is scaled by
`ω_e/c ≈ 2e-13`).
"""
function calc_pseudorange_curvature!(curvature, sat_positions, ξ, v)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    v_r = SVector{3}(v[1], v[2], v[3])
    v_r_sq = dot(v_r, v_r)
    for j in 1:satellite_count(sat_positions)
        sat_pos = get_satellite_position(sat_positions, j)
        u = satellite_position_at_reception(rₙ, sat_pos) - rₙ
        d = norm(u)
        û = u / d
        curvature[j] = (v_r_sq - dot(û, v_r)^2) / d
    end
    return curvature
end

"""
    calc_elevation_and_azimuth(user_pos::ECEF, sat_pos)
        -> (elevation, azimuth)
    calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, sat_pos)
        -> (elevation, azimuth)

Elevation and azimuth (radians) of the satellite at `sat_pos` (ECEF), seen
from the user position in its local East-North-Up frame; the azimuth runs clockwise from
North in `[0, 2π)`. Pass `enu_from_ecef = ENUfromECEF(user_pos, wgs84)` to share one
geodetic conversion across the satellites of an epoch.
"""
function calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, sat_pos)
    enu = SphericalFromCartesian()(enu_from_ecef(ECEF(sat_pos)))
    elevation = enu.ϕ
    # `SphericalFromCartesian` measures θ counter-clockwise from East (the ENU +x
    # axis); azimuth is measured clockwise from North, i.e. π/2 − θ.
    azimuth = mod2pi(π / 2 - enu.θ)
    return elevation, azimuth
end

calc_elevation_and_azimuth(user_pos::ECEF, sat_pos) =
    calc_elevation_and_azimuth(ENUfromECEF(user_pos, wgs84), sat_pos)
