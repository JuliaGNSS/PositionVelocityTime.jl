
"""
Computes ̂ρ, the distance between the satellite and the estimated user position 

$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_positions`: Satellite Positions
"""
function calc_ρ_hat!(ρ, sat_positions, ξ)
    rₙ = SVector{3}(ξ[1], ξ[2], ξ[3])
    tc = ξ[4]
    n = size(sat_positions, 2)
    for j in 1:n
        sat_pos = SVector{3}(sat_positions[1, j], sat_positions[2, j], sat_positions[3, j])
        travel_time = norm(sat_pos - rₙ) / SPEEDOFLIGHT
        rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
        ρ[j] = norm(rotated_sat_pos - rₙ) + tc
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


$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_positions`: Matrix of satellite positions
"""
function calc_H!(H, sat_positions, ξ)
    n = size(sat_positions, 2)
    for j in 1:n
        sat_pos = SVector{3}(view(sat_positions, :, j))
        e = calc_e(sat_pos, ξ)
        H[j, 1] = e[1]
        H[j, 2] = e[2]
        H[j, 3] = e[3]
        H[j, 4] = 1.0
    end
    return H
end

calc_H(sat_positions, ξ) =
    calc_H!(Matrix{Float64}(undef, size(sat_positions, 2), 4), sat_positions, ξ)

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
    n = size(sat_positions, 2)
    for j in 1:n
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

$SIGNATURES
`H_GEO`: Geometry matrix
"""
function calc_DOP(H_GEO)
    D = inv(SMatrix{4,4}(H_GEO' * H_GEO))
    if D[4, 4] < 0 ||
        D[3, 3] < 0 ||
        D[1, 1] + D[2, 2] < 0 ||
        D[1, 1] + D[2, 2] + D[3, 3] < 0 ||
        sum(diag(D)) < 0
        # Something has gone wrong
        # This could probably be detected somewhere else
        # more efficiently.
        return DOP(-1, -1, -1, -1, -1)
    end
    TDOP = sqrt(D[4, 4]) # temporal dop
    VDOP = sqrt(D[3, 3]) # vertical dop
    HDOP = sqrt(D[1, 1] + D[2, 2]) # horizontal dop
    PDOP = sqrt(D[1, 1] + D[2, 2] + D[3, 3]) # position dop
    GDOP = sqrt(sum(diag(D))) # geometrical dop

    return DOP(GDOP, PDOP, VDOP, HDOP, TDOP)
end

"""
Computes user position

$SIGNATURES
sat_positions: Array of satellite positions. Needs 3 values per satellite (xyz), size must be (3, N)")
`ρ`: Array of pseudo ranges 

Calculates the user position by least squares method. The algorithm is based on the common reception method. 
"""
function user_position(sat_positions, ρ, prev_ξ = zeros(4))
    sat_positions_mat = reduce(hcat, sat_positions)

    # Geodesic acceleration helps when starting far from the optimum (cold start
    # from origin, ~6e6 m away) by trading per-iteration work for fewer iterations.
    # When prev_ξ is already near-converged, the extra Avv! evals are pure overhead.
    # Detect cold by checking the default zeros(4) sentinel.
    ξ_fit_ols = if iszero(prev_ξ)
        curve_fit(
            calc_ρ_hat!, calc_H!, sat_positions_mat, ρ, collect(prev_ξ);
            inplace = true,
            avv! = (dir_deriv, p, v) -> calc_Avv!(dir_deriv, sat_positions_mat, p, v),
            lambda = 0.0,
            min_step_quality = 0.0,
        )
    else
        curve_fit(
            calc_ρ_hat!, calc_H!, sat_positions_mat, ρ, collect(prev_ξ);
            inplace = true,
        )
    end
    #    wt = 1 ./ (ξ_fit_ols.resid .^ 2)
    #    ξ_fit_wls = curve_fit(ρ_hat, H, sat_positions_mat, ρ, wt, collect(prev_ξ))
    rmse = sqrt(mean(ξ_fit_ols.resid .^ 2))
    return ξ_fit_ols.param, rmse
end

"""
Computes user velocity

$SIGNATURES

Calculates the user velocity and clock drift
"""
function calc_user_velocity_and_clock_drift(sat_positions_and_velocities, ξ, states, times, H)
    center_frequency = get_center_frequency(first(states).system)
    λ = SPEEDOFLIGHT / upreferred(center_frequency / Hz)
    n = length(states)
    y = Vector{Float64}(undef, n)
    for j in 1:n
        state = states[j]
        sat_pv = sat_positions_and_velocities[j]
        doppler = upreferred(state.carrier_doppler / Hz)
        clock_drift = calc_satellite_clock_drift(state.decoder, times[j])
        e = calc_e(get_sat_position(sat_pv), ξ)
        y[j] = -(doppler * λ - clock_drift * SPEEDOFLIGHT - dot(e, get_sat_velocity(sat_pv)))
    end
    _solve_velocity_qr(Val(n), H, y)
end

@generated function _solve_velocity_qr(::Val{n}, H, y) where {n}
    quote
        Hs = SMatrix{$n, 4}(H)
        ys = SVector{$n}(y)
        F = qr(Hs)
        SMatrix{4,4}(F.R) \ (F.Q' * ys)
    end
end

get_sat_position(x) = x.position
get_sat_velocity(x) = x.velocity
