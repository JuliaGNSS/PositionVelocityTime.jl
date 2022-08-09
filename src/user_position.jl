
"""
Computes ̂ρ, the distance between the satellite and the estimated user position 

$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_positions`: Satellite Positions
"""
function ρ_hat(sat_positions, ξ)
    rₙ = ξ[1:3]
    tc = ξ[4]
    
    map(eachcol(sat_positions)) do sat_pos
        travel_time = norm(sat_pos - rₙ) / SPEEDOFLIGHT
        rotated_sat_pos = rotate_by_earth_rotation(sat_pos, travel_time)
        norm(rotated_sat_pos - rₙ) + tc
    end
end


function rotate_by_earth_rotation(sat_pos, Δt)
    ω_e = 7.2921151467e-5
    α = ω_e * Δt
    Rz = @SMatrix [  cos(α) sin(α)   0;
                    -sin(α) cos(α)   0;
                    0       0        1];
    Rz * sat_pos
end

"""
Computes e, the direction vector from satellite to user

$SIGNATURES
`ξ`: Combination of estimated user position and time correction
`sat_pos`: Single satellite positions
"""
function e(sat_pos, ξ)
    rₙ = ξ[1:3]
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
function H(sat_positions, ξ)
    mapreduce(sat_pos -> [transpose(e(sat_pos, ξ)) 1], vcat, eachcol(sat_positions))
end



"""
Calculates the dilution of precision for a given geometry matrix H

$SIGNATURES
`H_GEO`: Geometry matrix
"""
function calc_DOP(H_GEO)
    D = inv(H_GEO' * H_GEO)
    TDOP = sqrt(D[4,4]) # temporal dop
    VDOP = sqrt(D[3,3]) # vertical dop
    HDOP = sqrt(D[1,1] + D[2,2]) # horizontal dop
    PDOP = sqrt(D[1,1] + D[2,2] + D[3, 3]) # position dop
    GDOP = sqrt(sum(diag(D))) # geometrical dop

    return DOP(GDOP,PDOP,VDOP,HDOP,TDOP)
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

    ξ_fit_ols = curve_fit(ρ_hat, H, sat_positions_mat, ρ, collect(prev_ξ))
#    wt = 1 ./ (ξ_fit_ols.resid .^ 2)
#    ξ_fit_wls = curve_fit(ρ_hat, H, sat_positions_mat, ρ, wt, collect(prev_ξ))

    return ξ_fit_ols.param
end