function correctSatTime(tTX,eph,Ek)
#Corrects the time transmitted by satellites with the data from the nav message
    #Parameters
    mu = 3.986005e14
    c  = 299792458
    F  = -2*sqrt(mu)/c^2
    #calculate satellite clock offset dtTX
    dtTX = tTX - eph.toc
    #relativistic correction term
    dtr = F * eph.e * eph.sqrt_A * sin(Ek)
    #satellite clock error
    tTXerror = eph.af0 + eph.af1*dtTX + eph.af2*dtTX^2 + dtr - eph.groupDelayDifferential
   return totc = tTX - tTXerror
end
