projected_position_ECI2ECEF = (
    [4.0186749839886995e6, 427051.19422150403, 4.918252576909507e6, -2.340933478024591e7], 17.961339568218193)


@testset "User position ECI2ECEF" begin
    positions = map( i -> sat_position_ECI_2_ECEF(test_dcs[i], test_cops[i], test_caps[i]), 1 : length(test_dcs))
    pseudo_ranges = PVT.pseudo_ranges(test_dcs, test_cops, test_caps) #pseudo_range(dcs, cps)
    out = user_position(positions, pseudo_ranges)
    @test out == projected_position_ECI2ECEF


    out = false
    try 
        calc_PVT(test_dcs[1:3], test_cops[1:3], test_caps[1:3])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true


    out = false
    try 
        calc_PVT(test_dcs, zeros(5), test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true
end


projected_position_ECEF = (
    [4.0186749839887144e6, 427051.1942215096, 4.918252576909532e6, -2.3409334780245904e7], 17.961339568218765)

@testset "User position ECEF" begin
    positions = map( i -> sat_position_ECEF(test_dcs[i], test_cops[i], test_caps[i]), 1 : length(test_dcs))
    pseudo_ranges = PVT.pseudo_ranges(test_dcs, test_cops, test_caps) #pseudo_range(dcs, cps)
    out = user_position(positions, pseudo_ranges)
    @test out == projected_position_ECEF


    

    out = false
    try 
        calc_PVT(test_dcs[1:3], test_cops[1:3], test_caps[1:3])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true


    out = false
    try 
        calc_PVT(test_dcs, zeros(5), test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true
end

