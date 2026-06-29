
"""
Computes ̂ρ, the distance between the satellite and the estimated user position

$SIGNATURES
`ξ`: Estimated user position, one clock correction per GNSS time system, and one
     receiver inter-frequency bias per frequency band beyond the reference band, i.e.
     `[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifb]`
`sat_positions`: Satellite Positions
`bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock column and
     inter-frequency-bias column of each satellite). Known range corrections (atmosphere,
     GGTO offset) are already folded into the measured `ρ` by [`calc_pvt`](@ref).
"""
function calc_ρ_hat!(ρ, sat_positions, ξ, bias_columns::BiasColumns)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    num_sats = size(sat_positions, 2)
    for j in 1:num_sats
        sat_pos = SVector{3}(sat_positions[1, j], sat_positions[2, j], sat_positions[3, j])
        travel_time = norm(sat_pos - rₙ) / SPEEDOFLIGHT
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
Computes e, the direction vector from satellite to user

$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_pos`: Single satellite positions
"""
function calc_e(sat_pos, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    travel_time = norm(sat_pos - rₙ) / SPEEDOFLIGHT
    rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
    (rₙ - rotated_sat_pos) / norm(rₙ - rotated_sat_pos)
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
`sat_positions`: Matrix of satellite positions
`bias_columns`: The per-satellite [`BiasColumns`](@ref)
"""
function calc_H!(H, sat_positions, ξ, bias_columns::BiasColumns)
    num_sats = size(sat_positions, 2)
    fill!(H, 0.0)
    for j in 1:num_sats
        sat_pos = SVector{3}(view(sat_positions, :, j))
        e = calc_e(sat_pos, ξ)
        H[j, 1] = e[1]
        H[j, 2] = e[2]
        H[j, 3] = e[3]
        H[j, 3+bias_columns.clock_bias_indices[j]] = 1.0
        if bias_columns.ifb_indices[j] != 0
            H[j, 3+bias_columns.num_clock_biases+bias_columns.ifb_indices[j]] = 1.0
        end
    end
    return H
end

calc_H(sat_positions, ξ, bias_columns::BiasColumns) =
    calc_H!(Matrix{Float64}(undef, size(sat_positions, 2), num_lsq_params(bias_columns)),
        sat_positions, ξ, bias_columns)

"""
Computes the directional second derivative of `calc_ρ_hat` along `v`,
used by LsqFit's geodesic acceleration.

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
    num_sats = size(sat_positions, 2)
    for j in 1:num_sats
        sat_pos = SVector{3}(sat_positions[1, j], sat_positions[2, j], sat_positions[3, j])
        travel_time = norm(sat_pos - rₙ) / SPEEDOFLIGHT
        rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
        u = rotated_sat_pos - rₙ
        d = norm(u)
        û = u / d
        dir_deriv[j] = (v_r_sq - dot(û, v_r)^2) / d
    end
    return dir_deriv
end

"""
Calculates the dilution of precision for a given geometry matrix H

`H_GEO` has `3 + num_clock_biases + num_ifb` columns (three position partials, one per
GNSS time system, and one per frequency band beyond the reference). `D = (HᵀH)⁻¹` is
the parameter covariance in the units/frame of `H`. Because `calc_e` builds `H` from
ECEF line-of-sight vectors, the position block of `D` is in ECEF; it is rotated into
the local ENU frame at `user_pos` before the horizontal/vertical split, so `HDOP`/`VDOP`
are taken in the user's tangent plane (`PDOP`/`GDOP` are trace-invariant and unaffected
by the rotation). `GDOP` spans all parameters; `TDOP` reports the clock variance of the
primary (reference) system — see [`PVTSolution`](@ref) — while the other systems' clock
(inter-system-bias) and the inter-frequency-bias variances enter `GDOP` only.

A rank-deficient geometry makes `HᵀH` singular (not positive definite); a
non-throwing Cholesky detects this and the function returns the sentinel
`DOP(-1, …)` instead of erroring.

$SIGNATURES
`H_GEO`: Geometry matrix
`user_pos`: User ECEF position, for the ECEF→ENU rotation of the horizontal/vertical DOPs
`primary_clock_index`: Index (1…num_clock_biases) of the clock column whose variance is reported as TDOP
"""
function calc_DOP(H_GEO, user_pos::ECEF, primary_clock_index = 1)
    # HᵀH is symmetric positive definite iff H has full column rank. Cholesky is
    # the right factorization for an SPD matrix (≈2× cheaper than a general or
    # symmetric-indefinite inverse, and numerically faithful), and `check = false`
    # lets a rank-deficient (singular) geometry fail gracefully via `issuccess`
    # instead of throwing. The inverse of an SPD matrix is itself SPD, so the DOP
    # variances on the diagonal are then guaranteed non-negative.
    F = cholesky(Symmetric(H_GEO' * H_GEO); check = false)
    issuccess(F) || return DOP(-1, -1, -1, -1, -1)
    D = inv(F)

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
    D_enu = R * SMatrix{3,3}(@view D[1:3, 1:3]) * R'

    HDOP = sqrt(D_enu[1, 1] + D_enu[2, 2])   # horizontal dop (East² + North²)
    VDOP = sqrt(D_enu[3, 3])                 # vertical dop (Up)
    PDOP = sqrt(D[1, 1] + D[2, 2] + D[3, 3]) # position dop (trace-invariant)
    TDOP = sqrt(D[3+primary_clock_index, 3+primary_clock_index]) # temporal dop (reference system)
    GDOP = sqrt(tr(D))                       # geometrical dop (all parameters)

    return DOP(GDOP, PDOP, VDOP, HDOP, TDOP)
end

"""
Computes user position

$SIGNATURES
`sat_positions_mat`: Satellite positions as a `(3, N)` matrix (xyz per satellite).
`ρ`: Array of pseudo ranges

Calculates the user position by least squares method. The algorithm is based on the common reception method.

`bias_columns`: The per-satellite [`BiasColumns`](@ref) (clock and inter-frequency-bias columns).

Returns `(ξ, residuals)`: the solved state vector
`ξ = [x, y, z, tc₁, …, ifb₁, …]` and the per-satellite post-fit residual vector
(modeled minus measured pseudorange, metres), in the same satellite order as `ρ`.
"""
function user_position(sat_positions_mat, ρ, bias_columns::BiasColumns, prev_ξ = zeros(num_lsq_params(bias_columns)))
    model! = (out, x, par) -> calc_ρ_hat!(out, x, par, bias_columns)
    jacobian! = (J, x, par) -> calc_H!(J, x, par, bias_columns)

    # Geodesic acceleration helps when starting far from the optimum (cold start
    # from origin, ~6e6 m away) by trading per-iteration work for fewer iterations.
    # When prev_ξ is already near-converged, the extra Avv! evals are pure overhead.
    # Detect cold by checking the default zeros sentinel (origin position).
    ξ_fit_ols = if iszero(prev_ξ)
        curve_fit(
            model!, jacobian!, sat_positions_mat, ρ, collect(prev_ξ);
            inplace = true,
            avv! = (dir_deriv, par, v) -> calc_Avv!(dir_deriv, sat_positions_mat, par, v),
            lambda = 0.0,
            min_step_quality = 0.0,
        )
    else
        curve_fit(
            model!, jacobian!, sat_positions_mat, ρ, collect(prev_ξ);
            inplace = true,
        )
    end
    #    wt = 1 ./ (ξ_fit_ols.resid .^ 2)
    #    ξ_fit_wls = curve_fit(ρ_hat, H, sat_positions_mat, ρ, wt, collect(prev_ξ))
    return ξ_fit_ols.param, ξ_fit_ols.resid
end

"""
Computes user velocity

$SIGNATURES

Calculates the user velocity and a single receiver clock drift, returned as
`[vx, vy, vz, ċ]`. Unlike the position solve — which estimates one clock *bias*
per GNSS time system (the inter-system bias / GGTO is metre-level and must be
resolved) — a single clock *drift* is shared by all satellites regardless of
GNSS: the receiver has one oscillator, and the only per-system difference is the
drift of the inter-system time offset (e.g. the GGTO rate `A_1G`), which is
~1e-6 m/s — far below the Doppler velocity resolution. Using one common drift lets
every satellite constrain the four unknowns instead of spending a column per
system.

The velocity geometry is built from the position geometry matrix `H` rather than
recomputed: its first three columns are the line-of-sight unit vectors (which `H`
already holds), and the single clock-drift column is the row-collapse of `H`'s
per-system clock columns — these are unit indicators, so they sum to `1` in every
row. The Doppler wavelength is evaluated per satellite from its own carrier
frequency.
"""
function calc_user_velocity_and_clock_drift(sat_positions_and_velocities, states, times, H)
    num_sats = length(states)
    # Normal-equations form of the 4-unknown velocity + clock-drift least squares.
    # The velocity design row is [eₓ e_y e_z 1]: the line-of-sight unit vector (H's
    # first three columns) plus the single common clock-drift column (H's per-system
    # clock columns, unit indicators, collapsed to a 1). Accumulating HᵀH (4×4) and
    # Hᵀy (length 4) row by row keeps the (num_sats × 4) design matrix unmaterialised
    # and the solve a fixed 4×4 regardless of satellite count — no per-count
    # recompilation and no per-epoch heap allocation. (HᵀH is well-conditioned for
    # the position geometry, the same matrix calc_DOP already inverts.) The Doppler
    # wavelength is evaluated per satellite from its own carrier frequency.
    HtH = zero(SMatrix{4,4,Float64})
    Hty = zero(SVector{4,Float64})
    for j in 1:num_sats
        state = states[j]
        sat_pv = sat_positions_and_velocities[j]
        λ = SPEEDOFLIGHT / upreferred(get_center_frequency(state.system) / Hz)
        doppler = upreferred(state.carrier_doppler / Hz)
        clock_drift = calc_satellite_clock_drift(state.decoder, times[j])
        # Line-of-sight unit vector, already computed for the position solve and
        # stored in H's first three columns (calc_H) — no need to recompute calc_e.
        e = SVector{3}(view(H, j, 1:3))
        a = SVector(e[1], e[2], e[3], 1.0)
        yⱼ = -(doppler * λ - clock_drift * SPEEDOFLIGHT - dot(e, get_sat_velocity(sat_pv)))
        HtH += a * a'
        Hty += a * yⱼ
    end
    HtH \ Hty
end

get_sat_position(x) = x.position
get_sat_velocity(x) = x.velocity
