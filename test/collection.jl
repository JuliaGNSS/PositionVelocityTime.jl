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
        named = calc_pvt((gps = gps, galileo = gal); kw...)
        @test calc_pvt((gps, gal); kw...).position == named.position
        @test calc_pvt((other = gps, labels = gal); kw...).position == named.position
        @test calc_pvt(gps; kw...).position == calc_pvt((gps,); kw...).position
        as_dictionary = Dictionary([state.decoder.prn for state in gps], gps)
        @test calc_pvt(as_dictionary; kw...).position == calc_pvt(gps; kw...).position
        @test calc_pvt((gps = as_dictionary, galileo = gal); kw...).position ==
              named.position
    end

    @testset "observation order fixes reference-system ties and result order" begin
        gps_first = calc_pvt((gps[1:3], gal[1:3]); kw...)
        gal_first = calc_pvt((gal[1:3], gps[1:3]); kw...)
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
                calc_pvt(input; kw...)
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
        @test_throws ArgumentError PositionVelocityTime.select_ionospheric_parameters((
            gps = SignalObservation[],
        ))
    end

    @testset "collection is inferable and the solver compiles once" begin
        groups = (gps = gps, galileo = gal)
        measurements, _ =
            PositionVelocityTime.collect_measurements(groups; approximate_year = 2021)
        @test measurements isa Vector{PositionVelocityTime.PreparedObservation}
        @test isconcretetype(eltype(measurements))
        @test length(measurements) == length(gps) + length(gal)
        @test inferred_return_type(
            PositionVelocityTime.collect_measurements,
            Tuple{typeof(groups)},
        ).parameters[1] === Vector{PositionVelocityTime.PreparedObservation}
        @test inferred_return_type(calc_pvt, Tuple{typeof(groups)}) === PVTSolution
        for input in (gps, gal, (a = gps,), (gps, gal), (b = gal, a = gps))
            calc_pvt(input; approximate_year = 2021)
        end
        specializations =
            Base.specializations(only(methods(PositionVelocityTime._solve_pvt!)))
        @test count(!isnothing, collect(specializations)) == 1
    end

    @testset "empty collections retain their type" begin
        empty_gps = gps[1:0]
        for input in (empty_gps, (), (;), Dictionary(Int[], empty_gps))
            measurements, parameters =
                PositionVelocityTime.collect_measurements(input; approximate_year = 2021)
            @test isempty(measurements)
            @test isnothing(parameters)
            @test isnothing(PositionVelocityTime.decide_bias_layout(measurements))
            previous = calc_pvt(gps; kw...)
            @test calc_pvt(input, previous; kw...) === previous
            @test calc_pvt((input = empty_gps, gps = gps); kw...).position ==
                  previous.position
        end
    end

    @testset "a second solve does not allocate unboundedly" begin
        groups = (gps = gps, galileo = gal)
        previous = calc_pvt(groups; kw...)
        calc_pvt(groups, previous; kw...)
        @test (@allocated calc_pvt(groups, previous; kw...)) < 100 * 1024
    end
end
