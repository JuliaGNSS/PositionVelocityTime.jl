module PVT

    using DocStringExtensions, Parameters, FixedPointNumbers, GNSSDecoder

    #export init_decode,

    function calcSinglePosition(decRes,tRXref)
        #decRes: decoding results for navigation data. Hay que pasarle decode.data
        #decRes es un array. Cada elemento del array es del tipo ::GPSData
        #tRXref: receiver time

        #preallocation
        xSat = [NaN for i=1:length(decRes)]
        ySat = [NaN for i=1:length(decRes)]
        zSat = [NaN for i=1:length(decRes)]
        prMat = [NaN for i=1:length(decRes)]

        #calculate satellite position and pseudorange
        satInd = 0
        for decInd = 1:length(decRes)
        #calculate satellite position
                tTX = (decRes[decInd].tow - 1)*6 #seconds since the transmission of the previous week
                #calculate Ek with uncorrected satellite clock
                a,b,c,Ek = satPosition(decRes[decInd], tTX)
                #correct satellite clock with Ek
                tTXcorr = correctSatTime(tTX, decRes[decInd], Ek)
                #calculate ECEF coordinates
                xSatTmp,ySatTmp,zSatTmp,d = satPosition(decRes[decInd], tTXcorr)

                #Pseudorange Estimation
                pseudorange = ((deltaT)*299792458) / 1e3

                #store variables
                satInd = satInd + 1
                xSat[satInd,1] = xSatTmp
                ySat[satInd,1] = ySatTmp
                zSat[satInd,1] = zSatTmp
                prMat[satInd,1] = pseudorange
        end

        #calculate user position
        uPos, dop = userPosition(xSat, ySat, zSat, prMat)
    end

    include("satPosition.jl")
    include("correctSatTime.jl")
    include("userPosition.jl")

end #end module
