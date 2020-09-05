projected_pos = (
    [4.0186749839887144e6, 427051.1942215096, 4.918252576909532e6, -2.3409334780245904e7], 17.961339568218765)



@testset "Testing calc_PVT" begin
        user = calc_PVT(test_dcs, test_cops, test_caps)
        @test user == projected_pos
end