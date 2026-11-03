@testset "observation preparation API" begin
    PVT = PositionVelocityTime
    gps = gps_l1_states(0.0Hz)
    # Supply broadcast coefficients to test model selection during preparation.
    gps = map(gps) do observation
        data = GNSSDecoder.GPSLNAVData(
            observation.decoder.data;
            α_0 = 1.118e-8,
            α_1 = 7.451e-9,
            α_2 = -5.961e-8,
            α_3 = -5.961e-8,
            β_0 = 90112.0,
            β_1 = 0.0,
            β_2 = -196608.0,
            β_3 = -65536.0,
        )
        SignalObservation(;
            decoder = GNSSDecoderState(observation.decoder; data),
            signal = observation.signal,
            code_phase = observation.code_phase,
            carrier_doppler = observation.carrier_doppler,
        )
    end
    observations = (gps = gps, galileo = galileo_e1b_states(0.0Hz))
    (; prepared, ionospheric_model, inter_system_time_offsets) =
        PVT.prepare_observations(observations; approximate_year = 2021)
    model = ionospheric_model
    @test model isa PVT.IonosphericModel
    @test model.parameters isa KlobucharParameters
    workspace = PVTWorkspace()
    output = solved_pvt!(PVTSolution(), workspace, observations, PVTSolution();
        approximate_year = 2021)
    @test prepared == workspace.prepared
    buffer = PVT.PreparedObservation[]
    buffered = PVT.prepare_observations!(buffer, observations; approximate_year = 2021)
    @test buffer == prepared
    @test buffered.ionospheric_model.parameters == model.parameters
    @test buffered.inter_system_time_offsets == inter_system_time_offsets

    @testset "failure statuses and output ownership" begin
        previous = copy(output)
        for scarce in (gps[1:0], gps[1:3])
            result = calc_pvt!(output, PVTWorkspace(), scarce, previous;
                approximate_year = 2021)
            @test result.status === PVT.insufficient_observations
            @test result.solution === output
            @test output.time == previous.time
            @test collect(pairs(output.observations)) ==
                  collect(pairs(previous.observations))
        end
        # Reject duplicate broadcast identities even when the epoch is too small to solve.
        @test_throws ArgumentError calc_pvt([gps[1], gps[1]]; approximate_year = 2021)
        @test_throws ArgumentError calc_pvt!(output, PVTWorkspace(),
            [gps[1], gps[1]], previous; approximate_year = 2021)
        # Preparation rejects them, so an external estimator never receives them either.
        @test_throws ArgumentError PVT.prepare_observations([gps[1], gps[1]];
            approximate_year = 2021)
        @test_throws ArgumentError PVT.prepare_observations!(PVT.PreparedObservation[],
            (collection_of(gps[1:2]), collection_of(gps[2:2])); approximate_year = 2021)
        # Prepared records belong to the measurement-model API, not the solve entry points.
        @test_throws ArgumentError calc_pvt(prepared; approximate_year = 2021)
        @test_throws ArgumentError calc_pvt!(output, PVTWorkspace(),
            prepared, previous; approximate_year = 2021)
        @test output.time == previous.time
        @test collect(pairs(output.observations)) == collect(pairs(previous.observations))
    end

    @testset "a missing decoded time anchor is not ready for positioning" begin
        first_observation = first(gps)
        unanchored = SignalObservation(;
            decoder = GNSSDecoderState(
                first_observation.decoder;
                num_bits_since_time_reference = nothing,
            ),
            signal = first_observation.signal,
            code_phase = first_observation.code_phase,
            carrier_doppler = first_observation.carrier_doppler,
        )
        @test !has_positioning_data(unanchored.decoder)
        @test isempty(PVT.prepare_observations([unanchored];
            approximate_year = 2021).prepared)
        @test calc_pvt([unanchored]; approximate_year = 2021).status ===
              PVT.insufficient_observations
    end
end
