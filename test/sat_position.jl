projected_positions_ECI2ECEF = [
    [-8.537268174201572e6, -1.2988350094779432e7, 2.1851582648954894e7],
    [7.172876999458652e6, -2.1714421682282157e7, 1.3296862881729363e7],
    [2.401937266433255e7, 1.1645127140291367e7, 214656.1633842468],
    [-5.021627715111267e6, 2.2734776424851105e7, 1.2106563114258626e7]
    ]
@testset "Sat Position ECI2ECEF" begin
    for i in 1:length(test_dcs)
        out = sat_position_ECI_2_ECEF(test_dcs[i], test_cops[i], test_caps[i])
        @test out == projected_positions_ECI2ECEF[i]
    end
end


projected_positions_ECEF = [
    [-8.537268174201585e6, -1.2988350094779423e7, 2.1851582648954894e7],
    [7.172876999458676e6, -2.171442168228215e7, 1.3296862881729363e7],
    [2.4019372664332554e7, 1.1645127140291365e7, 214656.1633842468],
    [-5.021627715111248e6, 2.273477642485111e7, 1.2106563114258626e7]
    ]
@testset "Sat Position ECEF" begin
    for i in 1:length(test_dcs)
        out = sat_position_ECEF(test_dcs[i], test_cops[i], test_caps[i])
        @test out == projected_positions_ECEF[i]
    end
end
