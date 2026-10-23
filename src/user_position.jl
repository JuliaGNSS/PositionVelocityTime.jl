
# The satellite positions every function below takes come in either of two layouts: a
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

# The rank test of `factorize_normal_matrix!` on the extreme pivots of the Cholesky factor of a
# normal-equations matrix `HᵀH`, which is positive definite exactly when `H` has full
# column rank — so a failure means the geometry left some parameter unobservable, and
# the epoch cannot be solved.
#
# A successful factorisation alone is not a sufficient test: a rank-deficient design
# usually does not produce an exactly zero pivot at all. Rounding leaves a tiny positive
# one instead, and the factorisation "succeeds" with a solution made of rounding noise
# (velocities of 1e8 m/s and the like). For a normal-equations matrix the Cholesky
# diagonal is on the scale of `H`'s singular values, so its smallest-to-largest ratio is
# ~`1/cond(H)`. A rank-deficient design leaves the ratio at rounding level (~1e-8 for
# `Float64`, or exactly 0), whereas even a barely usable GNSS geometry stays above
# ~1e-3; `cbrt(eps)` ≈ 6e-6 sits between the two with orders of magnitude of margin
# either way.
passes_rank_tolerance(min_pivot, max_pivot) =
    min_pivot > cbrt(eps(typeof(max_pivot))) * max_pivot

# `HᵀH` for the design `H`, Cholesky-factorised in place into the upper triangle of
# `normal_matrix` (`n × n` for the `n` columns of `H`). `false` when `H` does not have
# full column rank to within `passes_rank_tolerance`, and the factor is then garbage.
function factorize_normal_matrix!(normal_matrix, H)
    mul!(normal_matrix, transpose(H), H)
    _, info = LAPACK.potrf!('U', normal_matrix)
    info == 0 || return false
    passes_rank_tolerance(extrema(i -> normal_matrix[i, i], 1:size(H, 2))...)
end

"""
    calc_dop(H, user_position::ECEF, reference_clock_column = 1) -> DOP

The dilution of precision of the design matrix `H` ([`calc_design_matrix`](@ref)) at the
user ECEF position `user_position`.

`H` has `3 + num_clock_biases + num_inter_frequency_biases` columns (three position
partials, one per GNSS time system, and one per frequency band beyond the reference).
`D = (HᵀH)⁻¹` is the parameter covariance in the units/frame of `H`. Because
`calc_design_matrix!` builds `H` from ECEF unit vectors (the negated receiver→satellite
lines of sight), the position block of `D` is in ECEF; it is rotated into the local ENU
frame at `user_position` before the horizontal/vertical split, so `HDOP`/`VDOP` are taken in
the user's tangent plane (`PDOP`/`GDOP` are trace-invariant and unaffected by the rotation).
`GDOP` spans all parameters; `TDOP` reports the clock variance of the reference system, the
clock column `reference_clock_column` (1…`num_clock_biases`) — see [`PVTSolution`](@ref) —
while the other systems' clock (inter-system-bias) and the inter-frequency-bias variances
enter `GDOP` only.

A rank-deficient geometry makes `HᵀH` singular (not positive definite, or only to
within rounding); the factorisation detects this and the function returns the sentinel
`DOP(-1, …)` instead of erroring.
"""
calc_dop(H, user_position::ECEF, reference_clock_column = 1) = calc_dop!(
    Matrix{Float64}(undef, size(H, 2), size(H, 2)), H, user_position,
    reference_clock_column)

"""
    calc_dop!(normal_matrix, H, user_position::ECEF, reference_clock_column = 1) -> DOP

[`calc_dop`](@ref) with the normal-equations matrix `HᵀH` formed, factorised and
inverted in `normal_matrix` (`n × n` for the `n` columns of `H`, overwritten), so
that it allocates nothing. The factorisation is LAPACK's Cholesky `potrf`, followed by
a relative rank test on its pivots, and the inverse is `potri` on that factor — what
`inv` of a `Cholesky` computes.
"""
function calc_dop!(normal_matrix, H, user_position::ECEF, reference_clock_column = 1)
    # HᵀH is symmetric positive definite iff H has full column rank, so a
    # rank-deficient (singular) geometry fails gracefully here instead of throwing —
    # see `factorize_normal_matrix!`. The inverse of an SPD matrix is itself SPD, so
    # the DOP variances on the diagonal are then guaranteed non-negative.
    n = size(H, 2)
    factorize_normal_matrix!(normal_matrix, H) || return DOP(-1, -1, -1, -1, -1)
    # The inverse, in the upper triangle only; `D(i, j)` reads it symmetrically.
    LAPACK.potri!('U', normal_matrix)
    D(i, j) = i <= j ? normal_matrix[i, j] : normal_matrix[j, i]

    # Rotate the ECEF position covariance into the local ENU (East, North, Up)
    # frame so the horizontal/vertical split is taken in the user's tangent plane.
    # `R` is the ECEF→ENU rotation at the user position; it matches Geodesy's
    # `ENUfromECEF` convention (verified equal), built explicitly here as it is
    # marginally cheaper than recovering it from the transform.
    lla = LLAfromECEF(wgs84)(user_position)
    sφ, cφ = sincosd(lla.lat)
    sλ, cλ = sincosd(lla.lon)
    R = @SMatrix [
        -sλ     cλ      0.0
        -sφ*cλ  -sφ*sλ  cφ
        cφ*cλ   cφ*sλ   sφ
    ]
    D_position = SMatrix{3,3}(D(1, 1), D(2, 1), D(3, 1), D(1, 2), D(2, 2), D(3, 2),
        D(1, 3), D(2, 3), D(3, 3))
    D_enu = R * D_position * R'

    HDOP = sqrt(D_enu[1, 1] + D_enu[2, 2])   # horizontal dop (East² + North²)
    VDOP = sqrt(D_enu[3, 3])                 # vertical dop (Up)
    PDOP = sqrt(D(1, 1) + D(2, 2) + D(3, 3)) # position dop (trace-invariant)
    # temporal dop (reference system)
    TDOP = sqrt(D(3 + reference_clock_column, 3 + reference_clock_column))
    GDOP = sqrt(sum(i -> D(i, i), 1:n))      # geometrical dop (all parameters)

    return DOP(GDOP, PDOP, VDOP, HDOP, TDOP)
end

"""
    estimate_position(satellite_positions, pseudoranges, bias_columns::BiasColumns,
                      prev_ξ = zeros(num_parameters(bias_columns))) -> (ξ, residuals)

The user position, receiver clock biases and inter-frequency biases that fit the
pseudoranges by least squares, starting from `prev_ξ`. The algorithm is based on the
common reception method.

- `satellite_positions`: Satellite positions, a `3 × N` matrix or a vector of `SVector{3}`s.
- `pseudoranges`: The satellites' pseudoranges followed by the known `offset` of each
  [`ClockConstraint`](@ref) of `bias_columns`; any other length is an `ArgumentError`.
- `bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock and
  inter-frequency-bias columns).

Returns `(ξ, residuals)`: the solved state vector
`ξ = [x, y, z, tc₁, …, ifb₁, …]` and the per-satellite post-fit residual vector
(measured minus modelled pseudorange, metres), in the same satellite order as
`pseudoranges`; the constraint rows get none.

`curve_fit` (LsqFit's, mirrored by [`LevenbergMarquardt`](@ref)) reports its own residual
as `model - data`, so the returned vector negates it.
Measured − modelled ("observed minus computed") is how GNSS software reports observation
residuals — RTKLIB's `rescode`, and GNSS-SDR and PocketSDR through it — and the negation
is the whole of the difference: it is applied to the converged fit, so the solve itself
is untouched.

The allocating form of [`estimate_position!`](@ref), on fresh buffers.
"""
estimate_position(satellite_positions, pseudoranges, bias_columns::BiasColumns,
    prev_ξ = zeros(num_parameters(bias_columns))) = estimate_position!(
    LMWorkspace(), Vector{Float64}(undef, satellite_count(satellite_positions)),
    satellite_positions, pseudoranges, bias_columns, prev_ξ)

"""
    estimate_position!(workspace::LMWorkspace, residuals, satellite_positions,
                       pseudoranges, bias_columns, prev_ξ) -> (ξ, residuals)

[`estimate_position`](@ref) on the buffers of `workspace`, writing the post-fit residuals
into `residuals` (resized to the satellite count). The returned `ξ` **is** `workspace.x`, so
it is overwritten by the next solve on the same workspace; `prev_ξ` may be that same vector,
to restart a solve from the previous one's solution.
"""
function estimate_position!(
    workspace::LMWorkspace,
    residuals,
    satellite_positions,
    pseudoranges,
    bias_columns::BiasColumns,
    prev_ξ,
)
    length(pseudoranges) == num_observations(satellite_positions, bias_columns) ||
        throw(ArgumentError("pseudoranges has $(length(pseudoranges)) entries, but the " *
            "$(satellite_count(satellite_positions)) satellites and " *
            "$(length(bias_columns.constraints)) clock constraints need " *
            "$(num_observations(satellite_positions, bias_columns))"))
    model! = (out, x, par) -> predict_pseudoranges!(out, x, par, bias_columns)
    jacobian! = (J, x, par) -> calc_design_matrix!(J, x, par, bias_columns)

    # Two departures from LsqFit's Levenberg–Marquardt defaults, both about the scale of
    # this particular problem:
    #
    #  - `x_tol` is a *relative* step tolerance — the iteration stops once a step falls
    #    below `x_tol·(x_tol + ‖ξ‖)` — and `‖ξ‖` here is dominated not by the position but
    #    by the clock bias, which carries the ~2e7 m of common range the pseudoranges are
    #    referred to. The 1e-8 default therefore calls it converged at a step of ~0.2 m,
    #    so a solve that starts a metre from the optimum — an ordinary warm start from the
    #    previous epoch — stops most of a metre short of it. 1e-13 puts the threshold at a
    #    ~2 µm step, still ~1e3 × the rounding resolution of `ξ`.
    #  - `lambda` is the initial (inverse) trust-region radius, and the default 10 damps
    #    the first steps to a fraction of the Gauss–Newton step — which is what makes the
    #    tolerance above bite so early. The pseudorange model is only mildly nonlinear, so
    #    the full step is a good step here: 1e-8 leaves the iteration effectively
    #    Gauss–Newton — measured at 2 to 5 iterations to nanometre level from 1 m, 100 km
    #    or a cold start, against 8 to 16 with the default or with damping switched off
    #    entirely — while keeping LM's machinery available, since a step that fails to
    #    improve still grows `lambda` from there.
    #
    # Geodesic acceleration additionally helps when starting far from the optimum (cold
    # start from origin, ~6e6 m away) by trading per-iteration work for fewer iterations.
    # When prev_ξ is already near-converged, the extra curvature evaluations are pure
    # overhead. Detect cold by checking the default zeros sentinel (origin position).
    ξ_fit_ols = if iszero(prev_ξ)
        curve_fit!(
            workspace, model!, jacobian!, satellite_positions, pseudoranges, prev_ξ;
            inplace = true,
            avv! = (curvature, par, v) ->
                calc_pseudorange_curvature!(curvature, satellite_positions, par, v),
            lambda = 1e-8,
            min_step_quality = 0.0,
            x_tol = 1e-13,
        )
    else
        curve_fit!(
            workspace, model!, jacobian!, satellite_positions, pseudoranges, prev_ξ;
            inplace = true,
            lambda = 1e-8,
            x_tol = 1e-13,
        )
    end
    # Negated into the caller's buffer rather than in place: `resid` is the workspace's
    # own residual buffer, which the next fit starts from.
    num_sats = satellite_count(satellite_positions)
    residuals = resize!(residuals, num_sats)
    residuals .= .-view(ξ_fit_ols.resid, 1:num_sats)
    return ξ_fit_ols.param, residuals
end

"""
    estimate_velocity(measurements, user_position,
                      prev_velocity_and_drift = zero(SVector{4,Float64}))
        -> (velocity_and_drift::SVector{4}, rate_residuals)

The user velocity and a single receiver clock drift, returned together with
the per-satellite post-fit range-rate residuals as
`([vx, vy, vz, ċ], rate_residuals)`. `ċ` is the clock drift in m/s, the rate of the clock
bias in metres, as [`PVTSolution`](@ref)'s `clock_drift` reports it.
Unlike the position solve — which estimates one clock *bias*
per GNSS time system (the inter-system bias / GGTO is metre-level and must be
resolved) — a single clock *drift* is shared by all satellites regardless of
GNSS: the receiver has one oscillator, and the only per-system difference is the
drift of the inter-system time offset (e.g. the GGTO rate `A_1G`), which is
~1e-6 m/s — far below the Doppler velocity resolution. Using one common drift lets
every satellite constrain the four unknowns instead of spending a column per
system.

The measured range rate is fitted to [`calc_pseudorange_rate`](@ref), the exact time
derivative of the pseudorange [`predict_pseudoranges!`](@ref) models (see
`calc_pseudorange_rate` for its two terms of order `v/c`), by one Gauss–Newton step from
`prev_velocity_and_drift` (`[vx, vy, vz, ċ]`, m/s; the resting receiver by default). The
model is so nearly linear that one step is the least-squares fit: within ~1e-11 m/s from a
start a few m/s off, such as the previous epoch's solution, and within ~1e-7 m/s from rest
even for an aircraft's speed and 30 km/s of clock drift — far below any Doppler's noise
either way. That model needs the satellite's travel time, its rate and its
Earth-rotation-corrected position, hence the solved `user_position` (ECEF, m — the position
block of the converged `ξ`).

The residuals are `measured − modelled` range rate (m/s): the measured rate minus
[`calc_pseudorange_rate`](@ref) at the solved state, to the same precision, the same
orientation as the pseudorange residuals of [`estimate_position`](@ref) and in the same
satellite order as `measurements`. They are the range-rate analogue of the post-fit
pseudorange residual: a per-satellite Doppler-consistency / outlier indicator. Measured
and modelled are both taken in that model's sense, in which a *receding* satellite reads
positive — the same quantity and sign as RTKLIB's `resdop` residual, and hence the
negative of the Doppler-signed range rate a tracking loop works in (see the note at the
residual loop).

Needs a geometry whose position design at `user_position` ([`calc_design_matrix!`](@ref))
has full column rank, as [`calc_dop`](@ref) establishes. Where the velocity design itself is
not of full rank to within the tolerance `calc_dop` applies, the solve does not throw but
returns `NaN` for the state and every residual.

The allocating form of [`estimate_velocity!`](@ref).
"""
estimate_velocity(measurements, user_position,
    prev_velocity_and_drift = zero(SVector{4,Float64})) = estimate_velocity!(
    Vector{Float64}(undef, length(measurements)), SVector{4,Float64}[], measurements,
    user_position, prev_velocity_and_drift)

"""
    estimate_velocity!(rate_residuals, design, measurements, user_position,
                       prev_velocity_and_drift = zero(SVector{4,Float64}))
        -> (velocity_and_drift::SVector{4}, rate_residuals)

[`estimate_velocity`](@ref) writing the range-rate residuals into
`rate_residuals`, with `design` — a `Vector{SVector{4,Float64}}` — as scratch for the
design rows; both are resized to the satellite count. The solve itself is a fixed `4 × 4`
static one and allocates nothing.
"""
function estimate_velocity!(rate_residuals, design, measurements,
    user_position, prev_velocity_and_drift = zero(SVector{4,Float64}))
    position = SVector{3,Float64}(user_position[1], user_position[2], user_position[3])
    prev_velocity = SVector{3,Float64}(
        prev_velocity_and_drift[1], prev_velocity_and_drift[2], prev_velocity_and_drift[3])
    rate_residuals = resize!(rate_residuals, length(measurements))
    design = resize!(design, length(measurements))
    # One Gauss–Newton step on `measured − calc_pseudorange_rate` from the start, in
    # normal-equations form: accumulating AᵀA (4×4) and Aᵀr (length 4) row by row keeps
    # the solve a fixed 4×4 regardless of satellite count — no per-count recompilation and
    # no per-epoch heap allocation. Every quantity the loop needs — the Doppler, the
    # carrier, the satellite clock drift, position and velocity — is a field of the flat
    # [`SatelliteMeasurement`](@ref) row, so no decoder is touched and this compiles once
    # for every constellation mix.
    #
    # The design row is the gradient of `calc_pseudorange_rate` at the start: the negated
    # rotated line of sight plus a ~3e-6 relative share from the travel-time rate and the
    # light-time factor, and a clock-drift entry of `1 − s/c` (the single common drift
    # column; H's per-system clock columns, unit indicators, collapsed). The model is
    # linear in `ċ` and quadratic in `v` only through the light-time factor times the
    # travel-time rate's share: curvatures of `~2|q|/c ≈ 1e-14` per (m/s)² in `v` and
    # `|q|/c` in `v·ċ`. One step from a start `Δ` away therefore leaves an error of order
    # `1e-14·Δ²`, amplified by the geometry: on the fixtures ≤ 1e-11 m/s for a start a few
    # m/s off, ≤ 4e-9 m/s from rest for up to 780 m/s and 3 km/s of clock drift, and
    # ≤ 6e-8 m/s for 780 m/s and 30 km/s. A second step would take that up and nothing
    # more, so the solve takes one.
    AtA = zero(SMatrix{4,4,Float64})
    Atr = zero(SVector{4,Float64})
    for j in eachindex(measurements)
        measurement = measurements[j]
        modelled, a = pseudorange_rate_and_gradient(position, prev_velocity,
            prev_velocity_and_drift[4], measurement.position, measurement.velocity,
            measurement.satellite_clock_drift)
        rⱼ = measured_pseudorange_rate(measurement) - modelled
        rate_residuals[j] = rⱼ
        design[j] = a
        AtA += a * a'
        Atr += a * rⱼ
    end
    # The design row is H's row times `T` — `T` keeps H's three position columns, sums its
    # clock columns into the drift column and drops the IFB columns — plus a perturbation:
    # `k·q − s·(q − e)/c` in the position entries and `−s/c` in the drift entry, with
    # `|q| ≲ 2e-6` and `|s|/c ≲ 3e-6` for a satellite's range rate of up to ~900 m/s, so a
    # few 1e-6 relative. `T` has orthogonal columns of norms 1, 1, 1, √M, so a full-rank `H`
    # gives a full-rank `H·T`, with `σ_min(H·T)/σ_max(H·T) ≥ σ_min(H)/(√M·σ_max(H))`. The
    # DOP check `calc_pvt` runs first (`passes_rank_tolerance`) only asks the Cholesky
    # pivots of `HᵀH` for a ratio above `cbrt(eps)` ≈ 6e-6, and the pivots bound the
    # singular values from one side only (the smallest pivot is at least `σ_min`). At that
    # margin the perturbation is of the same order, so a geometry the DOP check just passes
    # is not guaranteed a positive definite `AtA`. Any usable geometry is far from it — a
    # pivot ratio of 1e-3 already means a GDOP of ~1e3 — but rather than throw from inside
    # `calc_pvt` for one that is not, the solve applies the same rank test to its own
    # factor and returns `NaN`, which `calc_pvt` treats as an unsolvable epoch.
    factor = cholesky(Symmetric(AtA); check = false)
    if !(issuccess(factor) && passes_rank_tolerance(extrema(diag(factor.U))...))
        fill!(rate_residuals, NaN)
        return SVector{4,Float64}(NaN, NaN, NaN, NaN), rate_residuals
    end
    step = factor \ Atr

    # Post-fit residual, measured minus modelled: the residual at the start moved on by the
    # step along its design row — `calc_pseudorange_rate` at the solved state to the step's
    # precision above, without a second rotation per satellite. It is orthogonal to the
    # design, as a least-squares residual is.
    #
    # The measured rate carries `-doppler * λ`, so it — and this residual with it — runs
    # in the geometric range-rate sense (positive while the satellite recedes), not the
    # Doppler-signed sense of a tracking loop's `λ · carrier_doppler`. A consumer that
    # forms its own rate residual from a loop's Doppler and compares it with this one has
    # to negate one of the two; magnitudes agree either way. This matches RTKLIB's
    # `resdop`, which likewise residuates `-lam * D` against a receding-positive model.
    for j in eachindex(measurements)
        rate_residuals[j] -= dot(design[j], step)
    end
    return prev_velocity_and_drift + step, rate_residuals
end

# The measured range rate, receding-positive as the model runs.
@inline measured_pseudorange_rate(measurement) =
    -measurement.carrier_doppler * SPEED_OF_LIGHT / measurement.center_frequency
