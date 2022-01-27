@testset "End of Week Crossover" begin
    t = 0
    @test t == (PositionVelocityTime.check_crossover(deepcopy(t)))

    t = 350000
    @test t - 604800 == PositionVelocityTime.check_crossover(deepcopy(t))

    t = -350000
    @test t + 604800 == PositionVelocityTime.check_crossover(deepcopy(t))

end

correct_raw_times =  [
    208863.71771887812,
    208863.72169284965,
    208863.72290867585,
    208863.7173908643,
    208863.7304948417]


@testset "Calc uncorrected Time" begin
    for i in 1:length(satellite_states)
        out = PositionVelocityTime.calc_uncorrected_time(satellite_states[i])
        @test out == correct_raw_times[i]
    end
end