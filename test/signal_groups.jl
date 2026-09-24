# Signal groups: the container `calc_pvt` takes, the normalization that accepts the
# shorter spellings of it, the refusal of anything else, and the two properties the
# whole redesign exists for — one compiled solver for every constellation mix, and a
# collection pass that is statically dispatched when the groups are named.
#
# `SignalGroup` is written out fully qualified throughout, deliberately: `Tracking`
# exports a type of the same name and `test/tracking_ext.jl` brings it into this very
# session's `Main`.

@testset "signal groups" begin
    kw = (;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)
    gps_group = PositionVelocityTime.SignalGroup(GPSL1CA(), gps)
    gal_group = PositionVelocityTime.SignalGroup(GalileoE1B(), gal)

    @testset "normalization accepts every spelling of the same epoch" begin
        named = calc_pvt((gps = gps_group, galileo = gal_group); kw...)
        # A bare tuple of groups is numbered by position…
        @test calc_pvt((gps_group, gal_group); kw...).position == named.position
        # …and a lone group becomes the `:default` group, so a one-constellation solve
        # needs no NamedTuple around it.
        @test calc_pvt(gps_group; kw...).position ==
              calc_pvt((default = gps_group,); kw...).position
        # The numbering itself, which is what makes a bare tuple's order explicit.
        normalize = PositionVelocityTime._normalize_signal_groups
        @test keys(normalize((gps_group, gal_group))) == (:group1, :group2)
        @test keys(normalize(gps_group)) == (:default,)
        @test normalize((a = gps_group,)) === (a = gps_group,)
    end

    @testset "a group's satellites may be a Dictionary or a vector" begin
        as_dictionary = PositionVelocityTime.SignalGroup(
            GPSL1CA(), Dictionary([state.decoder.prn for state in gps], gps))
        @test calc_pvt(as_dictionary; kw...).position == calc_pvt(gps_group; kw...).position
        # A `Dictionary`-backed group rejects a duplicate PRN structurally, where the
        # old flat vector could only document the precondition.
        @test_throws Dictionaries.IndexError Dictionary(
            [state.decoder.prn for state in [gps[1], gps[1]]], [gps[1], gps[1]])
    end

    @testset "group order is significant, exactly as vector order was" begin
        # The primary system is the most-populated one, ties broken by first
        # appearance — so two groups of equal size swap which one dates the fix.
        # Three each: two systems need 3 + 2 measurements from as many distinct
        # satellites, so a pair each would not solve at all.
        three_gps = PositionVelocityTime.SignalGroup(GPSL1CA(), gps[1:3])
        three_gal = PositionVelocityTime.SignalGroup(GalileoE1B(), gal[1:3])
        gps_first = calc_pvt((a = three_gps, b = three_gal); kw...)
        gal_first = calc_pvt((a = three_gal, b = three_gps); kw...)
        @test gps_first.reference_system == GPST()
        @test gal_first.reference_system == GST()
        # And `sats` is built in flat order: group order × within-group order.
        @test collect(keys(gal_first.sats)) == [
            [(:GalileoE1B, state.decoder.prn) for state in gal[1:3]]
            [(:GPSL1CA, state.decoder.prn) for state in gps[1:3]]
        ]
    end

    @testset "a flat vector of satellite states is refused, with directions" begin
        # The 5.x input. It is not accepted and not silently converted: the signals a
        # pooled vector holds are a runtime property, so any grouping derived from it
        # would be inference-blind — which is the thing groups exist to remove. The
        # error says so and shows what to build instead.
        err = try
            calc_pvt([gps; gal]; kw...)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        message = sprint(showerror, err)
        @test occursin("takes signal groups", message)
        @test occursin("SignalGroup(GPSL1CA()", message)
        # And it is raised at the entry point, not somewhere inside the solver.
        @test occursin("calc_pvt", message)
        @test occursin("not a vector of `SatelliteState`s", message)
    end

    @testset "groups of bare state vectors are refused, with the same directions" begin
        # The first thing a 5.x call site is likely to be migrated to: the right keys,
        # but the vectors not yet wrapped. It must fail at the entry point with the
        # directions, not with a field error from deep inside the collection pass.
        for input in ((gps = gps, galileo = gal), (gps, gal), (gps = gps_group, galileo = gal))
            @test_throws ArgumentError calc_pvt(input; kw...)
            message = try
                calc_pvt(input; kw...)
            catch e
                sprint(showerror, e)
            end
            @test occursin("takes signal groups", message)
            @test occursin("SignalGroup(GPSL1CA()", message)
        end
        @test_throws ArgumentError PositionVelocityTime.select_ionospheric_correction(
            (gps = gps,))
        # An empty NamedTuple is an epoch with no groups, not a malformed one.
        @test PositionVelocityTime._normalize_signal_groups((;)) === (;)
    end

    @testset "collection is inferable and the solver compiles once" begin
        groups = (gps = gps_group, galileo = gal_group)
        rows, _ = PositionVelocityTime.collect_measurements(groups; approximate_year = 2021)
        # The rows are one concrete type with no type parameters, whatever the mix.
        @test rows isa Vector{PositionVelocityTime.SatelliteMeasurement}
        @test isconcretetype(eltype(rows))
        @test length(rows) == length(gps) + length(gal)
        # Named groups make the collection pass statically dispatched: the rows half of
        # its result is inferred exactly. (The ionospheric-correction half is a small
        # `Union` by design — `_solve_pvt` takes it `@nospecialize`d so that the model
        # choice cannot multiply the solver's compiled copies.)
        @test inferred_return_type(
            PositionVelocityTime.collect_measurements, Tuple{typeof(groups)}
        ).parameters[1] === Vector{PositionVelocityTime.SatelliteMeasurement}
        @test inferred_return_type(calc_pvt, Tuple{typeof(groups)}) === PVTSolution

        # The point of the flat row: one compiled body of the solver serves every
        # constellation mix. Solve six differently-shaped epochs and count the
        # specializations of `_solve_pvt` — there are two regardless, one for a live
        # ionospheric model and one for `nothing`, and neither depends on the mix.
        for input in (
            gps_group, gal_group, (a = gps_group,), (a = gps_group, b = gal_group),
            (b = gal_group, a = gps_group), (gps_group, gal_group),
        )
            calc_pvt(input; approximate_year = 2021)
        end
        specializations =
            Base.specializations(only(methods(PositionVelocityTime._solve_pvt)))
        @test count(!isnothing, collect(specializations)) <= 2
    end

    @testset "an epoch with no usable satellite is not a special case" begin
        empty_group = PositionVelocityTime.SignalGroup(GPSL1CA(), SatelliteState[])
        rows, correction =
            PositionVelocityTime.collect_measurements(empty_group; approximate_year = 2021)
        @test isempty(rows)
        @test isnothing(correction)
        # `decide_bias_layout` reports unsolvable, so `calc_pvt` hands `prev_pvt` back.
        @test isnothing(PositionVelocityTime.decide_bias_layout(rows))
        previous = calc_pvt(gps_group; kw...)
        @test calc_pvt(empty_group, previous; kw...) === previous
        @test calc_pvt((a = empty_group, b = gps_group); kw...).position ==
              previous.position
    end

    @testset "a second solve of the same epoch does not allocate unboundedly" begin
        # A regression ceiling, not a target: the flat row exists so a two-group solve
        # costs about what a one-group solve of the same satellites costs, instead of
        # the ~290 KB the pooled abstract vector used to spend on dispatch boxes. The
        # bound is generous enough to survive ordinary code motion and tight enough to
        # catch the abstraction coming back.
        groups = (gps = gps_group, galileo = gal_group)
        previous = calc_pvt(groups; kw...)
        calc_pvt(groups, previous; kw...)
        @test (@allocated calc_pvt(groups, previous; kw...)) < 100 * 1024
    end
end
