projected_pos = PVTSolution(
    ECEF([4.0186911757808267e6, 427026.7592380331, 4.91826139491441e6]), 
    -2.1332096090475585e7, 
    DOP(2.346836665386327, 0, 0, 0, 0), 1:length(satellite_states), [], [])

prev_PVT = PVTSolution()
min_elevation =  0

@testset "Testing calc_PVT" begin
    user = calc_PVT(satellite_states, prev_PVT, min_elevation)
    @test user.pos ≈ projected_pos.pos
    @test user.receiver_time_correction ≈ projected_pos.receiver_time_correction
    @test get_gdop(user) ≈ get_gdop(projected_pos)


    out = false
    try 
        calc_PVT(satellite_states[1:3], prev_PVT, min_elevation)
    catch e
        if typeof(e) == PositionVelocityTime.SmallData
            out = true
        end
    end
    @test out == true


    svs_bad = deepcopy(satellite_states)
    svs_bad[1].decoder_state.subframes_decoded = [0,0,0,0,0]
    svs_bad[2].decoder_state.subframes_decoded = [0,0,0,0,0]

    out = false
    try 
        calc_PVT(svs_bad, prev_PVT, min_elevation)
    catch e
        if typeof(e) == PositionVelocityTime.SmallData
            out = true
        end
    end
    @test out == true
end