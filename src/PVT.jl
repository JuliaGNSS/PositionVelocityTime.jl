module PVT

    using DocStringExtensions, Parameters, FixedPointNumbers

    #export init_decode,

    abstract type GNSSData end

    @with_kw mutable struct GPSData <: GNSSData
        IODC::Union{Nothing, String} = nothing
    end

    function calcSinglePosition(decRes,selPRN,tRXref,msInterpol)
        #decRes: decoding results for navigation data
        #selPRN: struct containing the ids of the selected satellite
        #tRXref: receiver reference time
        #msInterpol: interpolation time for the PVT

        #preallocation
        xSat = [NaN for i=1:length(selPRN)]
        ySat = [NaN for i=1:length(selPRN)]
        zSat = [NaN for i=1:length(selPRN)]
        prMat = [NaN for i=1:length(selPRN)]

        #calculate satellite position and pseudorange
        satInd = 0
        for decInd = 1:length(decRes)
        #calculate satellite position
                tTX = (decRes{decInd}.navigationData.tow - 1)*6 + 1e-3*msInterpol
                #calculate Ek with uncorrected satellite clock
                a,b,c,Ek = satPosition(decRes{decInd}.data, tTX)
                #correct satellite clock with Ek
                tTXcorr = correctSatTime(tTX, decRes{decInd}.data, Ek)
                #calculate ECEF coordinates
                xSatTmp,ySatTmp,zSatTmp,d = satPosition(decRes{decInd}.data, tTXcorr)

                #Pseudorange Estimation
                #tRX = calcUserTime(decRes{decInd}, msInterpol)

                #use first preamble time for first estimation of pseudorange
                pseudorange = (tRX - tRXref - tTXcorr)*299792458

                #store variables
                satInd = satInd + 1
                xSat[satInd,1] = xSatTmp
                ySat[satInd,1] = ySatTmp
                zSat[satInd,1] = zSatTmp
                prMat[satInd,1] = pseudorange
        end

        ecef.xSat = xSat'
        ecef.ySat = ySat'
        ecef.zSat = zSat'
        ecef.pseudorange = prMat'

        #calculate user position
        ecef.uPos, ecef.dop = userPosition(xSat, ySat, zSat, prMat)
    end

    include("satPosition.jl")
    include("calcDOP")
    include("correctSatTime.jl")
    include("userPosition.jl")

end #end module
