projected_pos = PVTSolution(
    ECEF([4.0187039543983485e6, 427061.6186769258, 4.9183106756403595e6]), 
    -2.133206220402679e7, 
    2.346849505593898)
    

@testset "Testing calc_PVT" begin
    user = calc_PVT(satellite_states)
    @test user.pos ≈ projected_pos.pos
    @test user.receiver_time_correction ≈ projected_pos.receiver_time_correction
    @test user.GDOP ≈ projected_pos.GDOP


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