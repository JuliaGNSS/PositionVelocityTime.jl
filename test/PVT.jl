projected_pos = PVTSolution(
    ECEF([
        4.018703954398348e6, 
        427061.6186769259, 
        4.918310675640358e6
    ]), 
    -2.1332062204026792e7, 
    2.3468495055938994)
    


@testset "Testing calc_PVT" begin
        user = calc_PVT(test_dcs, test_cops, test_caps)
        @test user == projected_pos



    out = false
    try 
        calc_PVT(test_dcs, test_cops[1:3], test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true



    out = false
    try 
        calc_PVT(test_dcs, test_cops[1:3], test_caps[1:3])
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true


    out = false
    try 
        calc_PVT(test_dcs[1:3], test_cops[1:3], test_caps)
    catch e
        if typeof(e) == PVT.IncompatibleData
            out = true
        end
    end
    @test out == true


    out = false
    try 
        calc_PVT(test_dcs[1:3], test_cops[1:3], test_caps[1:3])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true


    dcs = deepcopy(test_dcs)
    dcs[1].subframes_decoded = [0,0,0,0,0]
    dcs[2].subframes_decoded = [0,0,0,0,0]
    out = false
    try 
        calc_PVT(dcs, test_cops, test_caps)
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true
end