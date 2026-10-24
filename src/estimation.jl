# The rank test on the extreme Cholesky pivots of a normal matrix `HᵀH`. A successful
# factorisation is not enough: a rank-deficient `H` usually leaves a tiny positive pivot
# rather than a zero one, and a solution of rounding noise. The pivots are on the scale of
# `H`'s singular values, so their ratio is ~`1/cond(H)`: rounding level (~1e-8) when rank
# deficient, above ~1e-3 for even a barely usable geometry. `cbrt(eps)` ≈ 6e-6 sits
# between the two.
passes_rank_tolerance(min_pivot, max_pivot) =
    min_pivot > cbrt(eps(typeof(max_pivot))) * max_pivot

# `HᵀH`, Cholesky-factorised in place into the upper triangle of `normal_matrix`. `false`
# when `H` fails `passes_rank_tolerance`, and the factor is then garbage.
function factorize_normal_matrix!(normal_matrix, H)
    mul!(normal_matrix, transpose(H), H)
    _, info = LAPACK.potrf!('U', normal_matrix)
    info == 0 || return false
    passes_rank_tolerance(extrema(i -> normal_matrix[i, i], 1:size(H, 2))...)
end

"""
    calc_dop(H, user_position::ECEF, reference_clock_column = 1) -> DOP

The dilution of precision of the design matrix `H` ([`calc_design_matrix`](@ref)) at the
user position. The position covariance is rotated from ECEF into the local ENU frame for
`HDOP`/`VDOP`. `GDOP` spans all parameters; `TDOP` is the variance of the clock column
`reference_clock_column`, the reference system's (see [`PVTSolution`](@ref)).

A rank-deficient geometry (see `passes_rank_tolerance`) returns `NaN` for every DOP
instead of throwing, as [`estimate_velocity`](@ref) does for its state.
"""
calc_dop(H, user_position::ECEF, reference_clock_column = 1) = calc_dop!(
    Matrix{Float64}(undef, size(H, 2), size(H, 2)), H, user_position,
    reference_clock_column)

"""
    calc_dop!(normal_matrix, H, user_position::ECEF, reference_clock_column = 1) -> DOP

[`calc_dop`](@ref) working in `normal_matrix` (`n × n` for the `n` columns of `H`,
overwritten). Allocates nothing.
"""
function calc_dop!(normal_matrix, H, user_position::ECEF, reference_clock_column = 1)
    n = size(H, 2)
    factorize_normal_matrix!(normal_matrix, H) || return DOP(NaN, NaN, NaN, NaN, NaN)
    # The inverse, in the upper triangle only; `D(i, j)` reads it symmetrically.
    LAPACK.potri!('U', normal_matrix)
    D(i, j) = i <= j ? normal_matrix[i, j] : normal_matrix[j, i]

    # The ECEF→ENU rotation at the user position, as Geodesy's `ENUfromECEF` builds it.
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

The user position, clock biases and inter-frequency biases that fit `pseudoranges` by
least squares from `prev_ξ`, and the per-satellite post-fit residuals (m, measured minus
modelled; see [`SatelliteInfo`](@ref)).

- `satellite_positions`: a `3 × N` matrix or a vector of `SVector{3}`s.
- `pseudoranges`: the satellites' pseudoranges followed by each [`ClockConstraint`](@ref)'s
  `offset`; any other length is an `ArgumentError`.
- `bias_columns`: the [`BiasColumns`](@ref).

The allocating form of [`estimate_position!`](@ref).
"""
estimate_position(satellite_positions, pseudoranges, bias_columns::BiasColumns,
    prev_ξ = zeros(num_parameters(bias_columns))) = estimate_position!(
    LMWorkspace(), Vector{Float64}(undef, satellite_count(satellite_positions)),
    satellite_positions, pseudoranges, bias_columns, prev_ξ)

"""
    estimate_position!(workspace::LMWorkspace, residuals, satellite_positions,
                       pseudoranges, bias_columns, prev_ξ) -> (ξ, residuals)

[`estimate_position`](@ref) on the buffers of `workspace`. The returned `ξ` **is**
`workspace.x`, overwritten by the next solve; `prev_ξ` may be that same vector.
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

    # Two departures from LsqFit's defaults, both for this problem's scale:
    #
    #  - `x_tol` is relative to `‖ξ‖`, which the clock bias dominates (~2e7 m), so the
    #    1e-8 default stops at a ~0.2 m step — most of a metre short on a warm start.
    #    1e-13 stops at ~2 µm, still ~1e3 × the rounding of `ξ`.
    #  - `lambda = 10` damps the first steps well short of Gauss–Newton. The model is only
    #    mildly nonlinear, so 1e-8 makes it effectively Gauss–Newton: 2 to 5 iterations
    #    from 1 m, 100 km or a cold start, against 8 to 16, while a failing step still
    #    grows `lambda`.
    #
    # Geodesic acceleration pays off only from a cold start (the zero seed, ~6e6 m away).
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
    # `curve_fit` reports `model − data`; negated into the caller's buffer, as `resid` is
    # the workspace's own.
    num_sats = satellite_count(satellite_positions)
    residuals = resize!(residuals, num_sats)
    residuals .= .-view(ξ_fit_ols.resid, 1:num_sats)
    return ξ_fit_ols.param, residuals
end

"""
    estimate_velocity(measurements, user_position,
                      prev_velocity_and_drift = zero(SVector{4,Float64}))
        -> (velocity_and_drift::SVector{4}, rate_residuals)

The user velocity and receiver clock drift `[vx, vy, vz, ċ]` (m/s; `ċ` as
`PVTSolution.clock_drift`), and the post-fit range-rate residuals (m/s, measured minus
modelled; see [`SatelliteInfo`](@ref)), in the order of `measurements`.

One drift is shared by every system: the receiver has one oscillator, and the drift
between system time scales (e.g. the GGTO's `A_1G`, ~1e-6 m/s) is far below the Doppler
resolution.

The measured range rates are fitted to [`calc_pseudorange_rate`](@ref) by one
Gauss–Newton step from `prev_velocity_and_drift`, which is the least-squares fit to well
below any Doppler noise (~1e-11 m/s from a few m/s off, ~1e-7 m/s from rest).
`user_position` is the solved position (ECEF, m).

Where the velocity design fails the rank test `calc_dop` applies, returns `NaN` for the
state and every residual instead of throwing.

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

[`estimate_velocity`](@ref) into `rate_residuals`, with `design`, a
`Vector{SVector{4,Float64}}`, as scratch; both are resized to the satellite count.
Allocates nothing.
"""
function estimate_velocity!(rate_residuals, design, measurements,
    user_position, prev_velocity_and_drift = zero(SVector{4,Float64}))
    position = SVector{3,Float64}(user_position[1], user_position[2], user_position[3])
    prev_velocity = SVector{3,Float64}(
        prev_velocity_and_drift[1], prev_velocity_and_drift[2], prev_velocity_and_drift[3])
    rate_residuals = resize!(rate_residuals, length(measurements))
    design = resize!(design, length(measurements))
    # One Gauss–Newton step in normal-equations form, accumulated row by row, so the solve
    # stays a static 4×4 whatever the satellite count.
    #
    # The model is quadratic in `v` only through the light-time factor times `q`: a
    # curvature of ~1e-14 per (m/s)². One step from a start `Δ` away leaves ~`1e-14·Δ²`,
    # times the geometry: on the fixtures ≤ 1e-11 m/s from a few m/s off, ≤ 4e-9 m/s from
    # rest at 780 m/s and 3 km/s of drift, ≤ 6e-8 m/s at 30 km/s.
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
    # The design row is `H`'s row times `T` (position columns kept, clock columns summed,
    # IFB columns dropped) plus a perturbation of a few 1e-6. `T`'s column norms are
    # 1, 1, 1, √M, so a full-rank `H` gives a full-rank `H·T`. But the DOP check only
    # bounds `H`'s condition to `cbrt(eps)`, the same order as the perturbation, so a
    # geometry that just passes it is not guaranteed a definite `AtA`: the same test is
    # applied here, returning `NaN` rather than throwing from inside `calc_pvt`.
    factor = cholesky(Symmetric(AtA); check = false)
    if !(issuccess(factor) && passes_rank_tolerance(extrema(diag(factor.U))...))
        fill!(rate_residuals, NaN)
        return SVector{4,Float64}(NaN, NaN, NaN, NaN), rate_residuals
    end
    step = factor \ Atr

    # The post-fit residual: the start's residual moved on by the step along its design
    # row, which is the model at the solved state to the precision above.
    for j in eachindex(measurements)
        rate_residuals[j] -= dot(design[j], step)
    end
    return prev_velocity_and_drift + step, rate_residuals
end

# The measured range rate, receding-positive like the model (the negative of a tracking
# loop's `λ · carrier_doppler`).
@inline measured_pseudorange_rate(measurement) =
    -measurement.carrier_doppler * SPEED_OF_LIGHT / measurement.center_frequency
