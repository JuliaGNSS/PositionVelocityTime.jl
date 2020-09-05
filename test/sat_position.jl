projected_positions_ECI2ECEF = [
    [-8.537268174165908e6, -1.2988350094802791e7, 2.185158264895441e7],
    [7.17287699944867e6, -2.171442168229333e7, 1.3296862881716367e7],
    [2.4019372664334804e7, 1.1645127140286792e7, 214656.1633768615],
    [-5.0216277151356945e6, 2.2734776424831044e7, 1.210656311428586e7]
    ]
@testset "Sat Position ECI2ECEF" begin
    for i in 1:length(test_dcs)
        out = sat_position_ECI_2_ECEF(test_dcs[i], test_cops[i], test_caps[i])
        @test out == projected_positions_ECI2ECEF[i]
    end
end


projected_positions_ECEF = [
    [-8.537591620812675e6, -1.2988831384034544e7, 2.1851170228041403e7],
    [7.172251721705176e6, -2.1714364550845448e7, 1.3297293456733963e7],
    [2.4019373920070715e7, 1.1645124554449795e7, 214654.89263821093],
    [-5.02018407963489e6, 2.2734511020904932e7, 1.210766018629351e7]
    ]
@testset "Sat Position ECEF" begin
    for i in 1:length(test_dcs)
        out = sat_position_ECEF(test_dcs[i], test_cops[i], test_caps[i])
        @test out == projected_positions_ECEF[i]
    end
end
