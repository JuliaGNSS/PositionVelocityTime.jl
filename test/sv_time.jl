@testset "End of Week Crossover" begin
    t = 0
    @test t == (PVT.check_crossover(deepcopy(t)))

    t = 350000
    @test t == ((PVT.check_crossover(deepcopy(t))) - 604800)

    t = -350000
    @test t == ((PVT.check_crossover(deepcopy(t))) + 604800)

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
        out = PVT.calc_uncorrected_time(test_dcs[i], test_cops[i])
        @test out != correct_raw_times[i]
    end
end