projected_position_ECI2ECEF = PVTSolution(
    ECEF([4.018703954398344e6, 427061.6186769246, 4.918310675640356e6]), -2.1332062204026792e7, 2.346849505593898)


@testset "User position ECI2ECEF" begin
    positions = map(x -> sat_position_ECI_2_ECEF(x), satellite_states)
    pseudo_ranges = PVT.pseudo_ranges(satellite_states)
    out = user_position(positions, pseudo_ranges)
    @test out == projected_position_ECI2ECEF


    out = false
    try 
        calc_PVT([sv1_struct])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true
end


projected_position_ECEF = PVTSolution(
    ECEF([4.0187039543983485e6, 427061.6186769258, 4.9183106756403595e6]), 
    -2.133206220402679e7, 
    2.346849505593898)
@testset "User position ECEF" begin
    positions = map( x -> sat_position_ECEF(x), satellite_states)
    pseudo_ranges = PVT.pseudo_ranges(satellite_states)
    out = user_position(positions, pseudo_ranges)
    @test out == projected_position_ECEF


    out = false
    try 
        calc_PVT([sv1_struct])
    catch e
        if typeof(e) == PVT.SmallData
            out = true
        end
    end
    @test out == true
end

