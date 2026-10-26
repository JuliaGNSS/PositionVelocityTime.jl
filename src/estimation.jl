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

A rank-deficient geometry makes `HᵀH` singular (not positive definite, or only to
within rounding); the factorisation detects this and the function returns the sentinel
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
that it allocates nothing. The factorisation is LAPACK's Cholesky `potrf`, followed by
a relative rank test on its pivots, and the inverse is `potri` on that factor — what
`inv` of a `Cholesky` computes.
"""
function calc_DOP!(normal_matrix, H_GEO, user_pos::ECEF, primary_clock_index = 1)
    # HᵀH is symmetric positive definite iff H has full column rank, so a
    # rank-deficient (singular) geometry fails gracefully here instead of throwing —
    # see `factorize_normal_matrix!`. The inverse of an SPD matrix is itself SPD, so
    # the DOP variances on the diagonal are then guaranteed non-negative.
    n = size(H_GEO, 2)
    factorize_normal_matrix!(normal_matrix, H_GEO) || return DOP(-1, -1, -1, -1, -1)
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
`ρ`: Array of pseudoranges

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

The measured range rate is fitted to [`calc_pseudorange_rate`](@ref), the exact time
derivative of the pseudorange [`calc_ρ_hat!`](@ref) models (see `calc_pseudorange_rate`
for its two terms of order `v/c`), by one Gauss–Newton step from `start`
(`[vx, vy, vz, ċ]`, m/s; the resting receiver by default). The model is so nearly linear
that one step is the least-squares fit: within ~1e-11 m/s from a start a few m/s off,
such as the previous epoch's solution, and within ~1e-7 m/s from rest even for an
aircraft's speed and 30 km/s of clock drift — far below any Doppler's noise either way.
That model needs the satellite's travel time, its rate and its Earth-rotation-corrected
position, hence the solved `user_pos` (ECEF, m — the position block of the converged
`ξ`).

The residuals are `measured − modeled` range rate (m/s): the measured rate minus
[`calc_pseudorange_rate`](@ref) at the solved state, to the same precision, the same
orientation as the pseudorange residuals of [`user_position`](@ref) and in the same
satellite order as `measurements`. They are the range-rate analogue of the post-fit
pseudorange residual: a per-satellite Doppler-consistency / outlier indicator. Measured
and modeled are both taken in that model's sense, in which a *receding* satellite reads
positive — the same quantity and sign as RTKLIB's `resdop` residual, and hence the
negative of the Doppler-signed range rate a tracking loop works in (see the note at the
residual loop).

Needs a geometry whose position design at `user_pos` ([`calc_H!`](@ref)) has full
column rank, as [`calc_DOP`](@ref) establishes. Where the velocity design itself is not of
full rank to within the tolerance `calc_DOP` applies, the solve does not throw but returns
`NaN` for the state and every residual.

The allocating form of [`calc_user_velocity_and_clock_drift!`](@ref).
"""
calc_user_velocity_and_clock_drift(measurements, user_pos,
    start = zero(SVector{4,Float64})) = calc_user_velocity_and_clock_drift!(
    Vector{Float64}(undef, length(measurements)), SVector{4,Float64}[], measurements,
    user_pos, start)

"""
    calc_user_velocity_and_clock_drift!(rate_residuals, design, measurements, user_pos,
        start = zero(SVector{4,Float64})) -> (velocity_and_drift::SVector{4}, rate_residuals)

[`calc_user_velocity_and_clock_drift`](@ref) writing the range-rate residuals into
`rate_residuals`, with `design` — a `Vector{SVector{4,Float64}}` — as scratch for the
design rows; both are resized to the satellite count. The solve itself is a fixed `4 × 4`
static one and allocates nothing.
"""
function calc_user_velocity_and_clock_drift!(rate_residuals, design, measurements,
    user_pos, start = zero(SVector{4,Float64}))
    user_position = SVector{3,Float64}(user_pos[1], user_pos[2], user_pos[3])
    start_vel = SVector{3,Float64}(start[1], start[2], start[3])
    rate_residuals = resize!(rate_residuals, length(measurements))
    design = resize!(design, length(measurements))
    # One Gauss–Newton step on `measured − calc_pseudorange_rate` from `start`, in
    # normal-equations form: accumulating AᵀA (4×4) and Aᵀr (length 4) row by row keeps
    # the solve a fixed 4×4 regardless of satellite count — no per-count recompilation and
    # no per-epoch heap allocation. Every quantity the loop needs — the Doppler, the
    # carrier, the satellite clock drift, position and velocity — is a field of the flat
    # [`SatelliteMeasurement`](@ref) row, so no decoder is touched and this compiles once
    # for every constellation mix.
    #
    # The design row is the gradient of `calc_pseudorange_rate` at `start`: the negated
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
        modelled, a = pseudorange_rate_and_gradient(user_position, start_vel, start[4],
            measurement.position, measurement.velocity, measurement.clock_drift)
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

    # Post-fit residual, measured minus modeled: the residual at `start` moved on by the
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
    return start + step, rate_residuals
end

# The measured range rate, receding-positive as the model runs.
@inline measured_pseudorange_rate(measurement) =
    -measurement.carrier_doppler * SPEED_OF_LIGHT / measurement.center_frequency

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
