projected_position_ECI2ECEF = (
    [4.018674983951664e6, 
    427051.1942053369, 
    4.918252576825729e6, 
    -2.0419893533580706e7], 
    17.961339568104332
)

@testset "User position ECI2ECEF" begin
    positions = map( i -> sat_position_ECI_2_ECEF(test_dcs[i], test_cops[i], test_caps[i]), 1 : length(test_dcs))
    pseudo_ranges = pseudo_range(test_dcs, test_cops, test_caps) #pseudo_range(dcs, cps)
    out = user_pos(positions, pseudo_ranges)
    @test out == projected_position_ECI2ECEF

end


projected_position_ECEF = (
    [4.020733079511654e6, 
    427494.43969367474, 
    4.922104601077217e6, 
    -2.340813840566231e7], 
    17.972637959899416
    )

@testset "User position ECEF" begin
    positions = map( i -> sat_position_ECEF(test_dcs[i], test_cops[i], test_caps[i]), 1 : length(test_dcs))
    pseudo_ranges = pseudo_range(test_dcs, test_cops, test_caps) #pseudo_range(dcs, cps)
    out = user_pos(positions, pseudo_ranges)
    @test out == projected_position_ECI2ECEF
end

