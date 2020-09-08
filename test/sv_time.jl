@testset "End of Week Crossover" begin
    t = 0
    @test t == (PVT.check_crossover(deepcopy(t)))

    t = 350000
    @test t - 604800 == PVT.check_crossover(deepcopy(t))

    t = -350000
    @test t + 604800 == PVT.check_crossover(deepcopy(t))

end

correct_raw_times =  [
    208863.7577188781,
    208863.76169284966,
    208863.76290867585,
    208863.75739086428]


@testset "Calc uncorrected Time" begin
    for i in 1:length(test_dcs)
        out = PVT.calc_uncorrected_time(test_dcs[i], test_cops[i], test_caps[i])
        @test out == correct_raw_times[i]
    end

    for i in 1:length(test_dcs)
        out = PVT.calc_uncorrected_time(test_dcs[i], 0)
        @test out != correct_raw_times[i]
    end
end