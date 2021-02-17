projected_position_ECI2ECEF = PVTSolution(
    ECEF([
        4.0187039543983415e6, 
        427061.61867692455, 
        4.918310675640349e6
    ]),
    -2.1332062204026792e7, 
    2.346849505593897)


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
        calc_PVT(test_dcs, zeros(4), test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true
end


projected_position_ECEF = PVTSolution(
    ECEF([
        4.018703954398348e6, 
        427061.6186769259, 
        4.918310675640358e6
    ]), 
    -2.1332062204026792e7, 
    2.3468495055938994)
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
        calc_PVT(test_dcs, zeros(4), test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true
end

