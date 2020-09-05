projected_pos = (
    [4.020733079511654e6,
    427494.43969367474,
    4.922104601077217e6,
    -2.340813840566231e7],
    17.972637959899416
)



@testset "Testing calc_PVT" begin
        user = calc_PVT(test_dcs, test_cops, test_caps)
        @test user == projected_pos
end