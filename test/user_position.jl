projected_position_ECI2ECEF = PVTSolution(
    ECEF([4.0186869577910933e6, 427050.24116486, 4.91825870249732e6]), -2.133209701701266e7, 2.3468360693153523, 0, 0, 0, 0, length(satellite_states), nothing, nothing, nothing)
prev_PVT = null_PVT()
min_elevation =  0


@testset "User position ECI2ECEF" begin
    positions = map(x -> sat_position_ECI_2_ECEF(x), satellite_states)
    pseudo_ranges = PositionVelocityTime.pseudo_ranges(satellite_states)
    out = user_position(prev_PVT.pos, positions, pseudo_ranges)
    @test out.pos ≈ projected_position_ECI2ECEF.pos
    @test out.receiver_time_correction ≈ projected_position_ECI2ECEF.receiver_time_correction
    @test out.GDOP ≈ projected_position_ECI2ECEF.GDOP

    out = false
    try 
        calc_PVT(prev_PVT, [sv1_struct], min_elevation)
    catch e
        println(e)
        if typeof(e) == PositionVelocityTime.SmallData
            out = true
        end
    end
    @test out == true
end


projected_position_ECEF = PVTSolution(
    ECEF([4.0186869577910956e6, 427050.2411648569, 4.918258702497325e6]), 
    -2.133209701701266e7, 
    2.346836069315351, 0, 0, 0, 0, length(satellite_states), nothing, nothing, nothing)
@testset "User position ECEF" begin
    positions = map( x -> sat_position_ECEF(x), satellite_states)
    pseudo_ranges = PositionVelocityTime.pseudo_ranges(satellite_states)
    out = user_position(prev_PVT.pos, positions, pseudo_ranges)
    @test out.pos ≈ projected_position_ECEF.pos
    @test out.receiver_time_correction ≈ projected_position_ECEF.receiver_time_correction
    @test out.GDOP ≈ projected_position_ECEF.GDOP


    out = false
    try 
        calc_PVT(prev_PVT, [sv1_struct], min_elevation)
    catch e
        if typeof(e) == PositionVelocityTime.SmallData
            out = true
        end
    end
    @test out == true
end

