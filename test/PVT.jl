projected_pos = PVTSolution(
    ECEF([
        4.018703954398348e6, 
        427061.6186769259, 
        4.918310675640358e6
    ]), 
    -2.1332062204026792e7, 
    2.3468495055938994)
    


@testset "Testing calc_PVT" begin
        user = calc_PVT(satellite_states)
        @test user == projected_pos


    out = false
    try 
        calc_PVT(satellite_states[1:3])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true


    svs_bad = deepcopy(satellite_states)
    svs_bad[1].decoder_state.subframes_decoded = [0,0,0,0,0]
    svs_bad[2].decoder_state.subframes_decoded = [0,0,0,0,0]

    out = false
    try 
        calc_PVT(svs_bad)
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true
end