module PVT

    using DocStringExtensions, Parameters, FixedPointNumbers, GNSSDecoder

    export calcSinglePosition,
    satPosition

    function calcSinglePosition(decRes,code_phase)
        #decRes: decoding results for navigation data.
        #code_phase: code phases from tracking

        preamble_offset = [NaN for i=1:length(decRes)]
        for decInd = 1:length(decRes)
            preamble_offset[decInd,1] = decRes[decInd].found_preambles.preamble_pos
        end
#       min_offset = minimum(preamble_offset)
#       pos_min = findall(x -> x == min_offset, preamble_offset)

        #preallocation
        xSat = [NaN for i=1:length(decRes)]
        ySat = [NaN for i=1:length(decRes)]
        zSat = [NaN for i=1:length(decRes)]
        prMat = [NaN for i=1:length(decRes)]

        #calculate satellite position and pseudorange
        for decInd = 1:length(decRes)
        #calculate satellite position
                tTX = (decRes[decInd].data.tow - 1)*6 #seconds since the transmission of the previous week
                #calculate Ek with uncorrected satellite clock
                a,b,c,Ek = satPosition(decRes[decInd].data, tTX)
                #correct satellite clock with Ek
                tTXcorr = correctSatTime(tTX, decRes[decInd].data, Ek)
                #calculate ECEF coordinates
                xSatTmp,ySatTmp,zSatTmp,d = satPosition(decRes[decInd].data, tTXcorr)

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
