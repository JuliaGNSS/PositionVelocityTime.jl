#MAP POSITION ESTIMATOR

using LinearAlgebra

function userPosition(xSat,ySat,zSat, prMat)
#input: raw pseudoranges, satellite ephemeris information
#output: user position and user clock bias

# Satellite position
rSat = [[xSat[1],ySat[1],zSat[1]] [xSat[2],ySat[2],zSat[2]] [xSat[3],ySat[3],zSat[3]] [xSat[4],ySat[4],zSat[4]]]
nSat = size(rSat,2) #number of satellites

# Initial position estimator and time error
rEst = [-50, -50, -50] #user position
dT = 0 #user clock bias
c = 299792458

#Pseudorange estimation
rhoT = prMat

# Gauss-Newton based position estimation

for k = 1:10 #number of iterations
    # Geometry matrix
    e = [[0.0,0.0,0.0] [0.0,0.0,0.0] [0.0,0.0,0.0] [0.0,0.0,0.0]]
    for n = 1:nSat
    e[:,n] = (rEst-rSat[:,n]) ./ norm(rEst - rSat[:,n])
    end

    H = transpose(e)
    H= [H [1,1,1,1]]

    # Space-time vector
    xi = [rEst; dT]
    rhoXi = [0.0 0.0 0.0 0.0]
    for n = 1:nSat
    rhoXi[n] = norm(xi[1:end-1] - rSat[:,n])
    end

    dXi = (transpose(H) * H)^-1 * transpose(H) * transpose((rhoT' - [dT dT dT dT] - rhoXi)) 
    # Update
    rEst = rEst + dXi[1:end-1]
    dT = dT + dXi[end]
    println("iteration number: $k")
    println("estimated position: ")
    println(rEst)
    println("estimated time error: ")
    println(dT)

end
    dop = calcDOP(H)
    return rEst, dop
end

function calcDOP(H)
    #Calculates the dilution of precision for a given geometry matrix H
    D = inv(H'*H)
    TDOP = sqrt(D[4,4]) #temporal dop
    VDOP = sqrt(D[3,3]) #vertical dop
    HDOP = sqrt(D[1,1] + D[2,2]) #horizontal dop
    PDOP = sqrt(D[1,1] + D[2,2] +D[3,3]) #position dop
    GDOP = sqrt(sum(diag(D))) #geometrical dop
end
