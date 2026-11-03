# Concrete collections preserve specialization without a grouping wrapper.
@testset "observation collections" begin
    kw = (;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)

    @testset "single, tuple and named collections describe the same epoch" begin
        named = solved_pvt((gps = gps, galileo = gal); kw...)
        @test solved_pvt((gps, gal); kw...).position == named.position
        @test solved_pvt((other = gps, labels = gal); kw...).position == named.position
        @test solved_pvt(gps; kw...).position ==
              solved_pvt((gps,); kw...).position
        as_dictionary = Dictionary([state.decoder.prn for state in gps], gps)
        @test solved_pvt(as_dictionary; kw...).position == solved_pvt(gps; kw...).position
        @test solved_pvt((gps = as_dictionary, galileo = gal); kw...).position ==
              named.position
    end

    @testset "observation order fixes reference-system ties and result order" begin
        gps_first = solved_pvt((gps[1:3], gal[1:3]); kw...)
        gal_first = solved_pvt((gal[1:3], gps[1:3]); kw...)
        @test gps_first.reference_system === GPST()
        @test gal_first.reference_system === GST()
        @test collect(keys(gal_first.observations)) == [
            [(:GalileoE1B, state.decoder.prn) for state in gal[1:3]]
            [(:GPSL1CA, state.decoder.prn) for state in gps[1:3]]
        ]
    end

    @testset "abstract pooled vectors and malformed collections are refused" begin
        for input in
            ([gps; gal], SignalObservation[], (gps = gps, bad = 1), (gps, 1), 1, [1, 2])
            err = try
                solved_pvt(input; kw...)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        message = sprint(showerror, err)
        @test occursin("concrete", message)
        @test occursin("SignalObservation", message)
        @test occursin("calc_pvt", message)
        end
        @test_throws ArgumentError PositionVelocityTime.select_ionospheric_correction(
            (gps = SignalObservation[],))
    end

    @testset "collection is inferable and the solver compiles once" begin
        observations = (gps = gps, galileo = gal)
        prepared, _ =
            PositionVelocityTime.collect_measurements(observations; approximate_year = 2021)
        @test prepared isa Vector{PositionVelocityTime.PreparedObservation}
        @test isconcretetype(eltype(prepared))
        @test length(prepared) == length(gps) + length(gal)
        @test inferred_return_type(
            PositionVelocityTime.collect_measurements, Tuple{typeof(observations)}
        ).parameters[1] === Vector{PositionVelocityTime.PreparedObservation}
        @test inferred_return_type(calc_pvt, Tuple{typeof(observations)}) ===
              NamedTuple{(:status, :solution),Tuple{PVTStatus,PVTSolution}}
        for input in (
            gps, gal, (a = gps,), (gps, gal),
            (b = gal, a = gps))
            solved_pvt(input; approximate_year = 2021)
        end
        specializations =
            Base.specializations(only(methods(PositionVelocityTime._solve_pvt!)))
        @test count(!isnothing, collect(specializations)) == 1
    end

    @testset "empty collections retain their type" begin
        empty_gps = gps[1:0]
        for input in (empty_gps, (), (;), Dictionary(Int[], empty_gps))
        prepared, parameters, inter_system_time_offsets =
            PositionVelocityTime.collect_measurements(input; approximate_year = 2021)
        @test isempty(prepared)
        @test isnothing(parameters)
            @test isnothing(PositionVelocityTime.decide_bias_layout(prepared, inter_system_time_offsets))
        previous = solved_pvt(gps; kw...)
        result= calc_pvt(input, previous; kw...)
            @test result.status === PositionVelocityTime.insufficient_observations
            @test iszero(result.solution.position)
            @test isempty(result.solution.observations)
            @test solved_pvt((input = empty_gps, gps = gps); kw...).position ==
              previous.position
    end
    end

    @testset "a second solve does not allocate unboundedly" begin
        observations = (gps = gps, galileo = gal)
        previous = solved_pvt(observations; kw...)
        solved_pvt(observations, previous; kw...)
        @test (@allocated calc_pvt(observations, previous; kw...)) < 100 * 1024
    end
end
