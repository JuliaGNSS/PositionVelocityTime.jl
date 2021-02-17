projected_pos = PVTSolution(ECEF([
        4.0186749839887214e6, 
        427051.19422151416, 
        4.918252576909554e6
    ]), 
    -2.3409334780245896e7, 
    17.961339568218655)
    


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