@testset "prepared observation API" begin
    PVT = PositionVelocityTime
    gps = gps_l1_states(0.0Hz)
    # Supply broadcast coefficients so the public split tests model hand-over too.
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
            carrier_phase = observation.carrier_phase,
        )
    end
    groups = (gps = gps, galileo = galileo_e1b_states(0.0Hz))
    prepared, model = PVT.prepare_observations(groups; approximate_year = 2021)
    before = copy(prepared)
    @test model isa PVT.IonosphericModel
    @test model.parameters isa KlobucharParameters

    for ionosphere in (false, true), troposphere in (false, true)
        options = (;
            enable_ionospheric_correction = ionosphere,
            enable_tropospheric_correction = troposphere,
        )
        direct = solved_pvt(groups; approximate_year = 2021, options...)
        from_prepared = solved_pvt(prepared; ionospheric_model = model, options...)
        @test direct.position == from_prepared.position
        @test direct.velocity == from_prepared.velocity
        @test direct.time == from_prepared.time
        @test collect(pairs(direct.observations)) ==
              collect(pairs(from_prepared.observations))

        workspace = PVTWorkspace()
        output = solved_pvt!(
            PVTSolution(),
            workspace,
            prepared,
            PVTSolution();
            ionospheric_model = model,
            options...,
        )
        @test output.position == direct.position
        @test output.time == direct.time
        warm = solved_pvt!(
            output,
            workspace,
            prepared,
            output;
            ionospheric_model = model,
            options...,
        )
        # Both paths must recompute atmospheric delays from the same warm seed.
        # A cold fix evaluates them at its initial uncorrected position.
        warm_direct = solved_pvt(groups, direct; approximate_year = 2021, options...)
        @test norm(warm.position - warm_direct.position) < 1e-4
        allocations(out, ws, input, params) = @allocated calc_pvt!(
            out,
            ws,
            input,
            out;
            ionospheric_model = params,
            options...,
        )
        allocations(warm, workspace, prepared, model)
        @test allocations(warm, workspace, prepared, model) == 0
        # A direct prepared solve never fills or changes private preparation scratch.
        @test isempty(workspace.measurements)
        @test prepared == before
    end

    @testset "failure statuses and output ownership" begin
        output = solved_pvt(prepared; ionospheric_model = model)
        previous = copy(output)
        for scarce in (prepared[1:0], prepared[1:3])
            result = calc_pvt!(output, PVTWorkspace(), scarce, previous)
            @test result.status === PVT.insufficient_observations
            @test result.solution === output
            @test output.time == previous.time
            @test collect(pairs(output.observations)) ==
                  collect(pairs(previous.observations))
        end
        # Reject duplicate broadcast identities even when the epoch is too small to solve.
        @test_throws ArgumentError calc_pvt([prepared[1], prepared[1]])
        @test_throws ArgumentError calc_pvt!(
            output,
            PVTWorkspace(),
            [prepared[1], prepared[1]],
            previous,
        )
        @test output.time == previous.time
        @test collect(pairs(output.observations)) == collect(pairs(previous.observations))
    end

    @testset "a missing decoded time anchor is not a usable observation" begin
        first_observation = first(gps)
        unanchored = SignalObservation(;
            decoder = GNSSDecoderState(
                first_observation.decoder;
                num_bits_since_time_reference = nothing,
            ),
            signal = first_observation.signal,
            code_phase = first_observation.code_phase,
            carrier_doppler = first_observation.carrier_doppler,
            carrier_phase = first_observation.carrier_phase,
        )
        @test has_positioning_data(unanchored.decoder)
        kept, _ = PVT.prepare_observations([unanchored]; approximate_year = 2021)
        @test isempty(kept)
        @test calc_pvt([unanchored]; approximate_year = 2021).status ===
              PVT.insufficient_observations
    end
end
