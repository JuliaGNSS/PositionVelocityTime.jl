function userPosition(xSat,ySat,zSat,prMeas)
    #Estimates the user position using a least-squares approach
    c = 299792458
    uPos = [0;0;0]
    tRx = 0
    dPos = [Inf; Inf; Inf; Inf]

    while sum(map(abs,dPos))>1
        prCalc = sqrt((xSat - uPos[1]).^2 + (ySat - uPos[2]).^2 + (zSat - uPos[3]).^2) + c.*tRx

        #Difference between calculated and measured pseudorange
        dpr = prCalc - prMeas
        r  = sqrt((xSat-uPos[1]).^2 + (ySat-uPos[2]).^2 + (zSat-uPos[3]).^2)

        ax = (xSat-uPos[1])./r
        ay = (ySat-uPos[2])./r
        az = (zSat-uPos[3])./r

        H = [ax ay az ones(length(ax),1)] #dice que esta matriz es singular
        dPos = (H'*H)\H'*dpr
        uPos = uPos + dPos[1:3]
        tRx  = tRx  - dPos[4]/c
    end
    dop = calcDOP(H)
    return uPos,dop
end

include("calcDOP")
