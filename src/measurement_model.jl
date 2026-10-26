# The pseudorange and pseudorange-rate model: predicted pseudoranges, their design matrix
# and curvature, the Earth-rotation correction, the rate model, and the elevation and
# azimuth of a satellite.

# The satellite positions the model and its estimators take come in either of two layouts: a
# `3 × N` matrix (one column per satellite), or a vector of `SVector{3,Float64}` — what
# the solver itself keeps, since a satellite position is a fixed-size quantity and a
# vector of them is resized in place where a matrix cannot be.
satellite_count(sat_positions::AbstractMatrix) = size(sat_positions, 2)
satellite_count(sat_positions::AbstractVector) = length(sat_positions)
satellite_position(sat_positions::AbstractMatrix, j) =
    SVector{3}(sat_positions[1, j], sat_positions[2, j], sat_positions[3, j])
satellite_position(sat_positions::AbstractVector, j) = SVector{3,Float64}(sat_positions[j])

"""
Computes ̂ρ, the distance between the satellite and the estimated user position

$SIGNATURES
`ξ`: Estimated user position, one clock correction per GNSS time system, and one
     receiver inter-frequency bias per frequency band beyond the reference band, i.e.
     `[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifb]`
`sat_positions`: Satellite positions, a `3 × N` matrix or a vector of `SVector{3}`s
`bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock column and
     inter-frequency-bias column of each satellite). Known range corrections (atmosphere,
     GGTO offset) are already folded into the measured `ρ` by [`calc_pvt`](@ref).
"""
function calc_ρ_hat!(ρ, sat_positions, ξ, bias_columns::BiasColumns)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    for j in 1:satellite_count(sat_positions)
        sat_pos = satellite_position(sat_positions, j)
        ifb =
            bias_columns.ifb_indices[j] == 0 ? 0.0 :
            ξ[3+bias_columns.num_clock_biases+bias_columns.ifb_indices[j]]
        ρ[j] = norm(satellite_position_at_reception(rₙ, sat_pos) - rₙ) +
               ξ[3+bias_columns.clock_bias_indices[j]] + ifb
    end
    return ρ
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
# signal's travel time `|s − u|/c`. The pseudorange model (`calc_ρ_hat!`), its position
# partials (`calc_line_of_sight` for `calc_H!`, `calc_Avv!`) and its rate
# (`pseudorange_rate_and_gradient`) all take the Earth-rotation correction from here, so
# they describe one and the same satellite.
@inline reception_frame_rotation(user_pos, sat_pos) =
    earth_rotation(norm(sat_pos - user_pos) / SPEED_OF_LIGHT)

# The satellite at `sat_pos` in the ECEF frame of reception at `user_pos`.
@inline satellite_position_at_reception(user_pos, sat_pos) =
    reception_frame_rotation(user_pos, sat_pos) * sat_pos

# The range-rate model: the exact time derivative of the pseudorange `calc_ρ_hat!` models,
# which the velocity solve fits (`calc_user_velocity_and_clock_drift`) and an estimator of
# its own predicts with (`calc_pseudorange_rate`).

# The satellite's position and velocity in the ECEF frame of reception, as `calc_ρ_hat!`
# models the position: `s' = R(τ)·s`, rotated by the Earth's turn over the travel time.
# The rotation angle `ω_e·τ` changes with the travel time, so `s'` moves at
# `R(τ)·ṡ + ω_e·τ̇·(s'_y, −s'_x, 0)` — the rotated velocity plus the rate of the rotation
# itself (`dR/dα = R·J` with `J` the quarter turn about z, which commutes with `R`). `τ̇`
# is taken for a resting receiver, `û·v_sat/c`; its user-velocity share is linear and is
# carried separately (`travel_time_velocity_partial`).
@inline function satellite_state_at_reception(user_pos, sat_pos, sat_vel)
    to_sat = sat_pos - user_pos
    travel_time_rate = dot(to_sat, sat_vel) / (norm(to_sat) * SPEED_OF_LIGHT)
    # One rotation for both, the one `calc_ρ_hat!` turns the satellite by.
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
second of receiver clock, of the pseudorange [`calc_ρ_hat!`](@ref) models — the satellite
rotated by the Earth's turn over the travel time into the ECEF frame of reception, plus
the receiver clock bias, minus the satellite clock offset. It is the model the velocity
solve ([`calc_user_velocity_and_clock_drift`](@ref)) fits, and its rate residuals are the
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
# *negative*, which is what the design-matrix rows of `calc_H!` carry. `ξ` holds the user
# position in its first three entries.
function calc_line_of_sight(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    to_sat = satellite_position_at_reception(rₙ, sat_pos) - rₙ
    return to_sat / norm(to_sat)
end

"""
Computes Geometry Matrix H

The matrix has [`num_lsq_params`](@ref)`(layout)` columns: three for the position
partials, one per GNSS time system, and one per frequency band beyond the reference
band. Each row carries a `1.0` in its system's clock column (`3 + clock_bias_indices[j]`)
and, unless the satellite is on the reference band, a `1.0` in its band's
inter-frequency-bias column (`3 + num_clock_biases + ifb_indices[j]`); the remaining
bias columns are zero.

$SIGNATURES
`ξ`: Estimated user position, per-system clock corrections, and per-band inter-frequency biases
`sat_positions`: Satellite positions, a `3 × N` matrix or a vector of `SVector{3}`s
`bias_columns`: The per-satellite [`BiasColumns`](@ref)
"""
function calc_H!(H, sat_positions, ξ, bias_columns::BiasColumns)
    fill!(H, 0.0)
    for j in 1:satellite_count(sat_positions)
        sat_pos = satellite_position(sat_positions, j)
        # The position partial ∂ρ/∂r is the negative of the receiver→satellite
        # line of sight.
        e = calc_line_of_sight(sat_pos, ξ)
        H[j, 1] = -e[1]
        H[j, 2] = -e[2]
        H[j, 3] = -e[3]
        H[j, 3+bias_columns.clock_bias_indices[j]] = 1.0
        if bias_columns.ifb_indices[j] != 0
            H[j, 3+bias_columns.num_clock_biases+bias_columns.ifb_indices[j]] = 1.0
        end
    end
    return H
end

"""
    calc_H(sat_positions, ξ, bias_columns::BiasColumns) -> Matrix{Float64}

Allocating form of [`calc_H!`](@ref PositionVelocityTime.calc_H!): the
least-squares design matrix at the state `ξ`, one row per satellite, with
[`num_lsq_params`](@ref PositionVelocityTime.num_lsq_params)`(bias_columns)`
columns.
"""
calc_H(sat_positions, ξ, bias_columns::BiasColumns) =
    calc_H!(Matrix{Float64}(undef, satellite_count(sat_positions), num_lsq_params(bias_columns)),
        sat_positions, ξ, bias_columns)

"""
Computes the directional second derivative of `calc_ρ_hat` along `v`,
used by the Levenberg-Marquardt geodesic acceleration.

For each satellite j the residual is `r_j(ξ) = ‖s_j' - r_n‖ + tc - ρ_j`,
where s_j' is the Earth-rotation-corrected satellite position. Treating
s_j' as constant w.r.t. ξ (the rotation depends on ξ via travel time,
but ω_e/c ≈ 2e-13 makes that contribution negligible here), the Hessian
is block-diagonal: the position block is `(I - û û^T) / d_j`, the time
block is zero. So `v^T H_j v = (‖v_r‖² - (û · v_r)²) / d_j`.

$SIGNATURES
"""
function calc_Avv!(dir_deriv, sat_positions, ξ, v)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    v_r = SVector{3}(v[1], v[2], v[3])
    v_r_sq = dot(v_r, v_r)
    for j in 1:satellite_count(sat_positions)
        sat_pos = satellite_position(sat_positions, j)
        u = satellite_position_at_reception(rₙ, sat_pos) - rₙ
        d = norm(u)
        û = u / d
        dir_deriv[j] = (v_r_sq - dot(û, v_r)^2) / d
    end
    return dir_deriv
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
