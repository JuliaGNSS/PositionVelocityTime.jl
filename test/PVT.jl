projected_pos = (
    [4.0186749839887144e6, 427051.1942215096, 4.918252576909532e6, -2.3409334780245904e7], 17.961339568218765)



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