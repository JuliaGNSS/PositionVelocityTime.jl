function calcDOP(H)
    #Calculates the dilution of precision for a given geometry matrix H
    D = inv(H'*H)
    TDOP = sqrt(D[4,4]) #temporal dop
    VDOP = sqrt(D[3,3]) #vertical dop
    HDOP = sqrt(D[1,1] + D[2,2]) #horizontal dop
    PDOP = sqrt(D[1,1] + D[2,2] +D[3,3]) #position dop
    GDOP = sqrt(sum(diag(D))) #geometrical dop
end
