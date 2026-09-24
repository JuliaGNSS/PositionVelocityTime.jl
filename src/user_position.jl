
# The satellite positions every function below takes come in either of two layouts: a
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
        travel_time = norm(sat_pos - rₙ) / SPEED_OF_LIGHT
        rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
        ifb =
            bias_columns.ifb_indices[j] == 0 ? 0.0 :
            ξ[3+bias_columns.num_clock_biases+bias_columns.ifb_indices[j]]
        ρ[j] = norm(rotated_sat_pos - rₙ) + ξ[3+bias_columns.clock_bias_indices[j]] + ifb
    end
    return ρ
end

function rotate_by_earth_rotation(sat_pos, Δt)
    ω_e = 7.2921151467e-5
    α = ω_e * Δt
    Rz = @SMatrix [
        cos(α) sin(α) 0
        -sin(α) cos(α) 0
        0 0 1
    ]
    Rz * sat_pos
end

"""
Computes the unit line-of-sight vector from the user to the
(Earth-rotation-corrected) satellite position — the textbook receiver→satellite
direction. The pseudorange's position partial is its *negative*, which is what
the design-matrix rows of [`calc_H!`](@ref) carry.

$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_pos`: Single satellite positions
"""
function calc_line_of_sight(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    travel_time = norm(sat_pos - rₙ) / SPEED_OF_LIGHT
    rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
    (rotated_sat_pos - rₙ) / norm(rotated_sat_pos - rₙ)
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
        travel_time = norm(sat_pos - rₙ) / SPEED_OF_LIGHT
        rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
        u = rotated_sat_pos - rₙ
        d = norm(u)
        û = u / d
        dir_deriv[j] = (v_r_sq - dot(û, v_r)^2) / d
    end
    return dir_deriv
end

"""
    positive_definite_cholesky(A::Symmetric) -> Union{Cholesky,Nothing}

Cholesky factorization of `A`, or `nothing` when `A` is not positive definite to within
the rank tolerance below. Used to solve and invert the normal-equations matrices of this
package (`HᵀH` for a design matrix `H`), which are positive definite exactly when `H`
has full column rank — so `nothing` means the geometry left some parameter
unobservable, and the epoch cannot be solved.

Cholesky is the right factorization for a symmetric positive definite matrix (≈2×
cheaper than a general or symmetric-indefinite one, and numerically faithful), and
`check = false` reports the non-positive-definite case instead of throwing. `issuccess`
alone is not a sufficient test, for two reasons:

  - StaticArrays accepts a pivot of exactly zero (its check is `pivot ≥ 0`) and returns
    a "successful" factor whose diagonal contains that zero — the triangular solve that
    follows then throws the very `SingularException` this is meant to avoid.
  - A rank-deficient design usually does not produce an exactly zero pivot at all:
    rounding leaves a tiny positive one instead, and the factorization "succeeds" with a
    solution made of rounding noise (velocities of 1e8 m/s and the like).

Both are caught by a relative rank tolerance on the Cholesky diagonal, which for a
normal-equations matrix is on the scale of `H`'s singular values — so its
smallest-to-largest ratio is ~`1/cond(H)`. A rank-deficient design leaves the ratio at
rounding level (~1e-8 for `Float64`, or exactly 0), whereas even a barely usable GNSS
geometry stays above ~1e-3; `cbrt(eps)` ≈ 6e-6 sits between the two with orders of
magnitude of margin either way. The test costs a few comparisons and does not allocate.

$SIGNATURES
"""
function positive_definite_cholesky(A::Symmetric)
    F = cholesky(A; check = false)
    issuccess(F) || return nothing
    pivots = diag(F.U)
    passes_rank_tolerance(minimum(pivots), maximum(pivots)) ? F : nothing
end

# The rank test of `positive_definite_cholesky` on the extreme pivots of a factor, shared
# with the in-place factorisation of `calc_DOP!`.
passes_rank_tolerance(min_pivot, max_pivot) =
    min_pivot > cbrt(eps(typeof(max_pivot))) * max_pivot

"""
Calculates the dilution of precision for a given geometry matrix H

`H_GEO` has `3 + num_clock_biases + num_ifb` columns (three position partials, one per
GNSS time system, and one per frequency band beyond the reference). `D = (HᵀH)⁻¹` is
the parameter covariance in the units/frame of `H`. Because `calc_H!` builds `H` from
ECEF unit vectors (the negated receiver→satellite lines of sight), the position
block of `D` is in ECEF; it is rotated into
the local ENU frame at `user_pos` before the horizontal/vertical split, so `HDOP`/`VDOP`
are taken in the user's tangent plane (`PDOP`/`GDOP` are trace-invariant and unaffected
by the rotation). `GDOP` spans all parameters; `TDOP` reports the clock variance of the
primary (reference) system — see [`PVTSolution`](@ref) — while the other systems' clock
(inter-system-bias) and the inter-frequency-bias variances enter `GDOP` only.

A rank-deficient geometry makes `HᵀH` singular (not positive definite);
[`positive_definite_cholesky`](@ref) detects this and the function returns the sentinel
`DOP(-1, …)` instead of erroring.

$SIGNATURES
`H_GEO`: Geometry matrix
`user_pos`: User ECEF position, for the ECEF→ENU rotation of the horizontal/vertical DOPs
`primary_clock_index`: Index (1…num_clock_biases) of the clock column whose variance is reported as TDOP
"""
calc_DOP(H_GEO, user_pos::ECEF, primary_clock_index = 1) = calc_DOP!(
    Matrix{Float64}(undef, size(H_GEO, 2), size(H_GEO, 2)), H_GEO, user_pos,
    primary_clock_index)

"""
    calc_DOP!(normal_matrix, H_GEO, user_pos::ECEF, primary_clock_index = 1) -> DOP

[`calc_DOP`](@ref) with the normal-equations matrix `HᵀH` formed, factorised and
inverted in `normal_matrix` (`n × n` for the `n` columns of `H_GEO`, overwritten), so
that it allocates nothing. The factorisation and the rank test are those of
[`positive_definite_cholesky`](@ref), and the inverse is LAPACK's `potri` on that factor
— what `inv` of a `Cholesky` computes — so the DOP is the same.
"""
function calc_DOP!(normal_matrix, H_GEO, user_pos::ECEF, primary_clock_index = 1)
    # HᵀH is symmetric positive definite iff H has full column rank, so a
    # rank-deficient (singular) geometry fails gracefully here instead of throwing —
    # see `positive_definite_cholesky`. The inverse of an SPD matrix is itself SPD, so
    # the DOP variances on the diagonal are then guaranteed non-negative.
    n = size(H_GEO, 2)
    mul!(normal_matrix, transpose(H_GEO), H_GEO)
    _, info = LAPACK.potrf!('U', normal_matrix)
    info == 0 || return DOP(-1, -1, -1, -1, -1)
    min_pivot, max_pivot = extrema(i -> normal_matrix[i, i], 1:n)
    passes_rank_tolerance(min_pivot, max_pivot) || return DOP(-1, -1, -1, -1, -1)
    # The inverse, in the upper triangle only; `D(i, j)` reads it symmetrically.
    LAPACK.potri!('U', normal_matrix)
    D(i, j) = i <= j ? normal_matrix[i, j] : normal_matrix[j, i]

    # Rotate the ECEF position covariance into the local ENU (East, North, Up)
    # frame so the horizontal/vertical split is taken in the user's tangent plane.
    # `R` is the ECEF→ENU rotation at the user position; it matches Geodesy's
    # `ENUfromECEF` convention (verified equal), built explicitly here as it is
    # marginally cheaper than recovering it from the transform.
    lla = LLAfromECEF(wgs84)(user_pos)
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
    TDOP = sqrt(D(3 + primary_clock_index, 3 + primary_clock_index)) # temporal dop (reference system)
    GDOP = sqrt(sum(i -> D(i, i), 1:n))      # geometrical dop (all parameters)

    return DOP(GDOP, PDOP, VDOP, HDOP, TDOP)
end

"""
Computes user position

$SIGNATURES
`sat_positions_mat`: Satellite positions, as a `(3, N)` matrix (xyz per satellite) or a
vector of `SVector{3}`s.
`ρ`: Array of pseudo ranges

Calculates the user position by least squares method. The algorithm is based on the common reception method.

`bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock and inter-frequency-bias columns).

Returns `(ξ, residuals)`: the solved state vector
`ξ = [x, y, z, tc₁, …, ifb₁, …]` and the per-satellite post-fit residual vector
(measured minus modeled pseudorange, metres), in the same satellite order as `ρ`.

`curve_fit` (LsqFit's, mirrored by [`LevenbergMarquardt`](@ref)) reports its own residual
as `model - data`, so the returned vector negates it.
Measured − modeled ("observed minus computed") is how GNSS software reports observation
residuals — RTKLIB's `rescode`, and GNSS-SDR and PocketSDR through it — and the negation
is the whole of the difference: it is applied to the converged fit, so the solve itself
is untouched.

The allocating form of [`user_position!`](@ref), on fresh buffers.
"""
user_position(sat_positions_mat, ρ, bias_columns::BiasColumns,
    prev_ξ = zeros(num_lsq_params(bias_columns))) = user_position!(
    LMWorkspace(), Vector{Float64}(undef, length(ρ)), sat_positions_mat, ρ, bias_columns,
    prev_ξ)

"""
    user_position!(workspace::LMWorkspace, residuals, sat_positions, ρ, bias_columns, prev_ξ)
        -> (ξ, residuals)

[`user_position`](@ref) on the buffers of `workspace`, writing the post-fit residuals into
`residuals` (resized to the satellite count). The returned `ξ` **is** `workspace.x`, so it
is overwritten by the next solve on the same workspace; `prev_ξ` may be that same vector,
to restart a solve from the previous one's solution.
"""
function user_position!(
    workspace::LMWorkspace,
    residuals,
    sat_positions_mat,
    ρ,
    bias_columns::BiasColumns,
    prev_ξ,
)
    model! = (out, x, par) -> calc_ρ_hat!(out, x, par, bias_columns)
    jacobian! = (J, x, par) -> calc_H!(J, x, par, bias_columns)

    # Two departures from LsqFit's Levenberg-Marquardt defaults, both about the scale of
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
    #    the first steps to a fraction of the Gauss-Newton step — which is what makes the
    #    tolerance above bite so early. The pseudorange model is only mildly nonlinear, so
    #    the full step is a good step here: 1e-8 leaves the iteration effectively
    #    Gauss-Newton — measured at 2 to 5 iterations to nanometre level from 1 m, 100 km
    #    or a cold start, against 8 to 16 with the default or with damping switched off
    #    entirely — while keeping LM's machinery available, since a step that fails to
    #    improve still grows `lambda` from there.
    #
    # Geodesic acceleration additionally helps when starting far from the optimum (cold
    # start from origin, ~6e6 m away) by trading per-iteration work for fewer iterations.
    # When prev_ξ is already near-converged, the extra Avv! evals are pure overhead.
    # Detect cold by checking the default zeros sentinel (origin position).
    ξ_fit_ols = if iszero(prev_ξ)
        curve_fit!(
            workspace, model!, jacobian!, sat_positions_mat, ρ, prev_ξ;
            inplace = true,
            avv! = (dir_deriv, par, v) -> calc_Avv!(dir_deriv, sat_positions_mat, par, v),
            lambda = 1e-8,
            min_step_quality = 0.0,
            x_tol = 1e-13,
        )
    else
        curve_fit!(
            workspace, model!, jacobian!, sat_positions_mat, ρ, prev_ξ;
            inplace = true,
            lambda = 1e-8,
            x_tol = 1e-13,
        )
    end
    #    wt = 1 ./ (ξ_fit_ols.resid .^ 2)
    #    ξ_fit_wls = curve_fit(ρ_hat, H, sat_positions_mat, ρ, wt, collect(prev_ξ))
    # Negated into the caller's buffer rather than in place: `resid` is the workspace's
    # own residual buffer, which the next fit starts from.
    residuals = resize!(residuals, length(ρ))
    residuals .= .-ξ_fit_ols.resid
    return ξ_fit_ols.param, residuals
end

"""
Computes user velocity

$SIGNATURES

Calculates the user velocity and a single receiver clock drift, returned together with
the per-satellite post-fit range-rate residuals as `([vx, vy, vz, ċ], rate_residuals)`.
Unlike the position solve — which estimates one clock *bias*
per GNSS time system (the inter-system bias / GGTO is metre-level and must be
resolved) — a single clock *drift* is shared by all satellites regardless of
GNSS: the receiver has one oscillator, and the only per-system difference is the
drift of the inter-system time offset (e.g. the GGTO rate `A_1G`), which is
~1e-6 m/s — far below the Doppler velocity resolution. Using one common drift lets
every satellite constrain the four unknowns instead of spending a column per
system.

The residuals are `measured − modeled` range rate (m/s), the same orientation as the
pseudorange residuals of [`user_position`](@ref) and in the same satellite order as
`measurements`. They are the range-rate analogue of the post-fit pseudorange residual: a
per-satellite Doppler-consistency / outlier indicator. Measured and modeled are both
taken in `yⱼ`'s sense below, in which a *receding* satellite reads positive — the same
quantity and sign as RTKLIB's `resdop` residual, and hence the negative of the
Doppler-signed range rate a tracking loop works in (see the note at the residual loop).

Requires a geometry whose position design `H` has full column rank — the caller
establishes that with [`calc_DOP`](@ref) before calling this; see the comment at the solve.

The allocating form of [`calc_user_velocity_and_clock_drift!`](@ref).
"""
calc_user_velocity_and_clock_drift(measurements, H) = calc_user_velocity_and_clock_drift!(
    Vector{Float64}(undef, length(measurements)), measurements, H)

"""
    calc_user_velocity_and_clock_drift!(rate_residuals, measurements, H)
        -> (velocity_and_drift::SVector{4}, rate_residuals)

[`calc_user_velocity_and_clock_drift`](@ref) writing the range-rate residuals into
`rate_residuals` (resized to the satellite count). The solve itself is a fixed `4 × 4`
static one and allocates nothing.
"""
function calc_user_velocity_and_clock_drift!(rate_residuals, measurements, H)
    num_sats = length(measurements)
    # Normal-equations form of the 4-unknown velocity + clock-drift least squares.
    # The velocity design row is [eₓ e_y e_z 1]: the pseudorange's position partial
    # (H's first three columns, the negated receiver→satellite line of sight) plus the
    # single common clock-drift column (H's per-system clock columns, unit indicators,
    # collapsed to a 1). Accumulating HᵀH (4×4) and
    # Hᵀy (length 4) row by row keeps the (num_sats × 4) design matrix unmaterialised
    # and the solve a fixed 4×4 regardless of satellite count — no per-count
    # recompilation and no per-epoch heap allocation. The Doppler wavelength is
    # evaluated per satellite from its own carrier frequency. Every quantity the loop
    # needs — the Doppler, the carrier, the satellite clock drift and velocity — is a
    # field of the flat [`SatelliteMeasurement`](@ref) row, so no decoder is touched
    # and this compiles once for every constellation mix.
    HtH = zero(SMatrix{4,4,Float64})
    Hty = zero(SVector{4,Float64})
    # The normal-equations form does not keep the design rows, so the measurements are
    # kept here instead and turned into post-fit residuals in place after the solve —
    # cheaper than a second pass that recomputes each `yⱼ` from the row.
    rate_residuals = resize!(rate_residuals, num_sats)
    for j in 1:num_sats
        measurement = measurements[j]
        λ = SPEED_OF_LIGHT / measurement.center_frequency
        doppler = measurement.carrier_doppler
        clock_drift = measurement.clock_drift
        # The pseudorange's position partial — the negative of the
        # receiver→satellite line of sight — already computed for the position
        # solve and stored in H's first three columns (calc_H).
        e = SVector{3}(view(H, j, 1:3))
        a = SVector(e[1], e[2], e[3], 1.0)
        yⱼ = -(doppler * λ - clock_drift * SPEED_OF_LIGHT - dot(e, measurement.velocity))
        rate_residuals[j] = yⱼ
        HtH += a * a'
        Hty += a * yⱼ
    end
    # `HtH` is positive definite whenever the position design `H` has full column rank,
    # which `calc_pvt` has established via `calc_DOP` before calling this: writing the
    # velocity design as `A = H·T` — `T` keeps H's three position columns, sums its clock
    # columns into the drift column and drops the IFB columns — `T` has orthogonal columns
    # of norms 1, 1, 1, √M, so `null(T) = {0}` and a full-rank `H` gives a full-rank `A`
    # (with `cond(A) ≤ √M·cond(H)`). Hence a plain Cholesky solve, no rank test of its own.
    #
    # This is why `calc_pvt` must keep the DOP check *ahead* of this call: with the order
    # reversed, a rank-deficient geometry reaches the solve below and throws
    # `SingularException` — the failure this arrangement exists to prevent.
    velocity_and_drift = cholesky(Symmetric(HtH)) \ Hty

    # Post-fit residual, measured minus modeled: the measurement stored in the
    # accumulation loop, minus the design row `[e 1]` (rebuilt from H's
    # position-partial columns, as above) applied to the solved state.
    #
    # `yⱼ` carries `-doppler * λ`, so it — and this residual with it — runs in the
    # geometric range-rate sense (positive while the satellite recedes), not the
    # Doppler-signed sense of a tracking loop's `λ · carrier_doppler`. A consumer that
    # forms its own rate residual from a loop's Doppler and compares it with this one
    # has to negate one of the two; magnitudes agree either way. This matches RTKLIB's
    # `resdop`, which likewise residuates `-lam * D` against a receding-positive model.
    velocity = SVector{3}(
        velocity_and_drift[1], velocity_and_drift[2], velocity_and_drift[3])
    for j in 1:num_sats
        e = SVector{3}(view(H, j, 1:3))
        rate_residuals[j] -= dot(e, velocity) + velocity_and_drift[4]
    end
    return velocity_and_drift, rate_residuals
end

"""
    get_sat_position(sat_pv) -> SVector{3,Float64}

Satellite ECEF position (m) of one [`calc_satellite_position_and_velocity`](@ref) result.
"""
get_sat_position(x) = x.position

"""
    get_sat_velocity(sat_pv) -> SVector{3,Float64}

Satellite ECEF velocity (m/s) of one [`calc_satellite_position_and_velocity`](@ref) result.
"""
get_sat_velocity(x) = x.velocity
