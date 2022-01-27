projected_positions_ECI2ECEF = [
    [-8.53735973609123e6, -1.2988290225267656e7, 2.1851584427175216e7],
    [7.172861967216288e6, -2.171448849034389e7, 1.3296760809673272e7],
    [2.4019377846893914e7, 1.1645117858204663e7, 214527.98779279838],
    [-5.021601719182211e6, 2.2734837592605237e7, 1.2106460195444101e7],  
    [8.067370261443458e6, 1.280420712089501e7, 2.1814190881086327e7]
    ]
@testset "Sat Position ECI2ECEF" begin
    for i in 1:length(satellite_states)
        out = sat_position_ECI_2_ECEF(satellite_states[i])
        @test out ≈ projected_positions_ECI2ECEF[i]
    end

    sv_bad = deepcopy(sv1_struct)
    sv_bad.decoder_state.subframes_decoded = [0,1,1,1,1]

    out = false
    try 
        sat_position_ECEF(sv_bad)
    catch e
        if typeof(e) == PositionVelocityTime.BadData
            out = true
        end
    end
    @test out == true
end


projected_positions_ECEF = [
    [-8.537359736091228e6, -1.2988290225267658e7, 2.1851584427175216e7],
    [7.172861967216265e6, -2.1714488490343902e7, 1.3296760809673272e7],
    [2.40193778468939e7, 1.1645117858204681e7, 214527.98779279838],
    [-5.021601719182221e6, 2.2734837592605233e7, 1.2106460195444101e7],
    [8.067370261443453e6, 1.2804207120895013e7, 2.1814190881086327e7]
    ]
@testset "Sat Position ECEF" begin
    for i in 1:length(satellite_states)
        out = sat_position_ECEF(satellite_states[i])
        @test out == projected_positions_ECEF[i]
    end

    sv_bad = deepcopy(sv1_struct)
    sv_bad.decoder_state.subframes_decoded = [0,0,0,0,0]
    out = false
    try 
        sat_position_ECEF(sv_bad)
    catch e
        if typeof(e) == PositionVelocityTime.BadData
            out = true
        end
    end
    @test out == true
end
