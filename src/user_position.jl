ITERATIONS = 20

function user_position(pos, p_ranges)
        
        
        sizes = [length(pos), length(p_ranges)]
        nSat = minimum(sizes)[1]
        
        nSat >= 4 || throw(SmallData("At least 4 Satellites needed for position computing"))
        (length(pos) == length(p_ranges)) || throw(IncompatibleData("Length of Input Arrays must be equal"))
        size(pos)[2] == 3 || throw(InvalidData("Position Matrix needs 3 values per satellite (xyz), size must be (3, N)"))
        
        ξ = [0, 0, 0, 0]
        H = zeros(nSat, 4) # initialization of Geometry Matrix
        
        Δpos(nr) = [ξ[1] - pos[nr][1], ξ[2] - pos[nr][2], ξ[3] - pos[nr][3]]  # Computes the Distance between estimated user position and Satellite
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
    
    return ξ ,calc_DOP(H)
end



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
    return ξ
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