module PVT

    using DocStringExtensions, Parameters, FixedPointNumbers, GNSSDecoder

    export calcSinglePosition,
    satPosition

    function calcSinglePosition(decRes)
        #decRes: decoding results for navigation data.
        #tRXref: receiver time

        #preallocation
        xSat = [NaN for i=1:length(decRes)]
        ySat = [NaN for i=1:length(decRes)]
        zSat = [NaN for i=1:length(decRes)]
        prMat = [NaN for i=1:length(decRes)]

        #calculate satellite position and pseudorange
        for decInd = 1:length(decRes)
        #calculate satellite position
                tTX = (decRes[decInd].tow - 1)*6 #seconds since the transmission of the previous week
                #calculate Ek with uncorrected satellite clock
                a,b,c,Ek = satPosition(decRes[decInd], tTX)
                #correct satellite clock with Ek
                tTXcorr = correctSatTime(tTX, decRes[decInd], Ek)
                #calculate ECEF coordinates
                xSatTmp,ySatTmp,zSatTmp,d = satPosition(decRes[decInd], tTXcorr)

#TODO
                #Pseudorange Estimation

                #store variables

                xSat[decInd,1] = xSatTmp
                ySat[decInd,1] = ySatTmp
                zSat[decInd,1] = zSatTmp
                #prMat[decInd,1] = pseudorange
        end

        #calculate user position
        uPos, dop = userPosition(xSat, ySat, zSat, prMat)
    end

    include("satPosition.jl")
    include("correctSatTime.jl")
    include("userPosition.jl")

end #end module
