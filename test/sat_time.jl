@testset "Week crossover" begin
    @test PositionVelocityTime.correct_week_crossovers(0) == 0

    @test PositionVelocityTime.correct_week_crossovers(350000) == 350000 - 604800

    @test PositionVelocityTime.correct_week_crossovers(-350000) == 604800 - 350000
end