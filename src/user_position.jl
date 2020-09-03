ITERATIONS = 20

"""
Computes user position

$SIGNATURES
´pos´: Array of Satellite positions
´p_ranges´: Array of pseudo ranges 

Calculates the user position by least squares method. The algorithm is based on the common reception method. 
"""
function user_position(sv_pos, p_ranges)
        
        
        sizes = [length(sv_pos), length(p_ranges)]
        nSat = minimum(sizes)[1]
        
        nSat >= 4 || throw(SmallData("At least 4 Satellites needed for position computing"))
        (length(sv_pos) == length(p_ranges)) || throw(IncompatibleData("Length of Input Arrays must be equal"))
        size(sv_pos)[2] == 3 || throw(InvalidData("Position Matrix needs 3 values per satellite (xyz), size must be (3, N)"))
        
        ξ = [0, 0, 0, 0]
        H = zeros(nSat, 4) # initialization of Geometry Matrix
        
        Δpos(nr) = [ξ[1] - sv_pos[nr][1], ξ[2] - sv_pos[nr][2], ξ[3] - sv_pos[nr][3]]  # Computes the Distance between estimated user position and Satellite
        e(nr) = Δpos(nr) ./ norm(Δpos(nr)) 
        
        
        for i in 1:ITERATIONS
            # Geometry Matrix
            for i = 1:nSat
                H[i, 1:3] = e(i)
                H[i, 4] = 1
            end

            # Calculated pseudoranges
            calc_pr = map(i->norm(Δpos(i)), 1:nSat)


            Δξ = (transpose(H) * H)^-1 * transpose(H) * ((p_ranges[1:nSat] .- ξ[4]) - calc_pr)

            ξ = ξ + Δξ
    end
    
    return ξ, calc_DOP(H)
end


"""
Computes user position

$SIGNATURES
´x_sat´: Array of x positions of sat positions
´y_sat´: Array of y positions of sat positions
´z_sat´: Array of z positions of sat positions
´p_ranges´: Array of pseudo ranges 

Calculates the user position by least squares method. 
The algorithm is based on the common reception method. 
"""
function user_position(x_sat, y_sat, z_sat, p_ranges)
        
        
    sizes = [length(x_sat), length(x_sat), length(x_sat), length(p_ranges)]
    nSat = minimum(sizes)[1]
    

    (size(x_sat) == size(y_sat) == size(z_sat) == size(p_ranges)) || throw(IncompatibleData("Size of Input Arrays must be equal ((1 , N), N = Nmax)"))
    nSat >= 4 || throw(SmallData("At least 4 Satellites needed for position computing"))    
    
    ξ = [0, 0, 0, 0]
    H = zeros(nSat, 4) # initialization of Geometry Matrix
    
    Δpos(nr) = [ξ[1] - x_sat[nr], ξ[2] - y_sat[nr], ξ[3] - z_sat[nr]]  # Computes the Distance between estimated user position and Satellite
    e(nr) = Δpos(nr) ./ norm(Δpos(nr)) 
    
    
    for i in 1:ITERATIONS
        # Geometry Matrix
        for i = 1:nSat
            H[i, 1:3] = e(i)
            H[i, 4] = 1
        end

        # Calculated pseudoranges
        calc_pr = map(i->norm(Δpos(i)), 1:nSat)


        Δξ = (transpose(H) * H)^-1 * transpose(H) * ((p_ranges[1:nSat] .- ξ[4]) - calc_pr)

        ξ = ξ + Δξ
    end
    return ξ, calc_DOP(H)
end



function calc_DOP(H)
    # Calculates the dilution of precision for a given geometry matrix H
    D = inv(H' * H)
    TDOP = sqrt(D[4,4]) # temporal dop
    VDOP = sqrt(D[3,3]) # vertical dop
    HDOP = sqrt(D[1,1] + D[2,2]) # horizontal dop
    PDOP = sqrt(D[1,1] + D[2,2] + D[3,3]) # position dop
    GDOP = sqrt(sum(diag(D))) # geometrical dop
    return GDOP
end