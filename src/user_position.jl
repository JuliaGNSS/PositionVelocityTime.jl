"""
Creates empty PVTSolution
"""
function null_PVT()
    pvt = PVTSolution(ECEF([0,0,0]), 0, nothing, nothing, nothing, nothing, nothing, 0, nothing, nothing, nothing)
end

"""
Computes ̂ρ , the distance between the satellite and the assumed user position 

$SIGNATURES
´ξ´: Combination of estimated user position and time correction
´rSat´: Single satellite Positions
"""
function ρ_hat(ξ, rSats)
    rₙ = ξ[1:3]
    tc = ξ[4]
    
    ρ_hat = map(i -> norm(rSats[i] - rₙ) + tc, 1:length(rSats))
end

"""
Computes e, the direction vector from satellite to user

$SIGNATURES
´ξ´: Combination of estimated user position and time correction
´rSat´: Single satellite Positions
"""
function e(ξ, rSat)
    rₙ = ξ[1:3]
    e = (rₙ - rSat) / norm(rₙ - rSat)
end

"""
Computes Geometry Matrix H


$SIGNATURES
´ξ´: Combination of estimated user position and time correction
´rSats´: Matrix of satellite Positions
"""
function H(ξ, rSats)
    H = mapreduce(rSat -> [transpose(e(ξ, rSat)) 1], vcat, rSats)
end



"""
Calculates the dilution of precision for a given geometry matrix H

$SIGNATURES
´H_GEO´: Geometry matrix
"""
function calc_DOP(H_GEO)
    D = inv(H_GEO' * H_GEO)
    TDOP = sqrt(D[4,4]) # temporal dop
    VDOP = sqrt(D[3,3]) # vertical dop
    HDOP = sqrt(D[1,1] + D[2,2]) # horizontal dop
    PDOP = sqrt(D[1,1] + D[2,2] + D[3, 3]) # position dop
    GDOP = sqrt(sum(diag(D))) # geometrical dop
    return GDOP, PDOP, VDOP, HDOP, TDOP
end

"""
Computes user position

$SIGNATURES
´rSats´: Array of Satellite positions. needs 3 values per satellite (xyz), size must be (3, N)")
´ρ´: Array of pseudo ranges 

Calculates the user position by least squares method. The algorithm is based on the common reception method. 
"""
function user_position(prev_pos ,rSats, ρ, accuracy::Float64 = 0.02)

    length(rSats) >= 4 || throw(SmallData("Not enough usable SV Data in user_position estimation."))
        
    # First Guesses of Position (Center of WGS84, time error = 0)
    #r₀ = [0.0, 0.0, 0.0]
    r₀ = prev_pos
    t = 0.0
    ξ = [r₀; t]  

    Δξ = 100
    while norm(Δξ) > accuracy
        Δρ = ρ - ρ_hat(ξ, rSats)
        Δξ = H(ξ, rSats) \ Δρ #(transpose(H(ξ)) * H(ξ)) \ (transpose(H(ξ)) * Δρ) #H(ξ) \ Δρ
        ξ = ξ + Δξ # ξₙ₊₁ = ξₙ + Δξ
    end

    pos = ECEF(ξ[1:3])
    dt_receiver = ξ[4]
    GDOP, PDOP, VDOP, HDOP, TDOP = calc_DOP(H(ξ, rSats))
    PVT = PVTSolution(pos, dt_receiver, GDOP, PDOP, VDOP, HDOP, TDOP, length(ρ), nothing, nothing, nothing)
    return PVT
end