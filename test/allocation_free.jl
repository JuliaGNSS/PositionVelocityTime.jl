# The in-place solve: `calc_pvt!` writes the fix into a solution it is handed and works
# in a reusable `PVTWorkspace`. Three properties are pinned here — it is the same solve
# as `calc_pvt`, it overwrites exactly the solution it is given and nothing else, and once
# its buffers have grown to an epoch it allocates nothing at all.

using PositionVelocityTime:
    _precompile_states, _precompile_cnav, _precompile_cnav2,
    _precompile_fnav, _PRECOMPILE_GPS_L1CA_STATES, _PRECOMPILE_GALILEO_E1B_STATES

# The precompile workload's fixture satellites, re-encoded onto every navigation-data
# type the solver dispatches on — the same epochs the trimmed app solves.
precompile_gps(system, make_data = identity) = PositionVelocityTime.SignalGroup(
    system, _precompile_states(system, _PRECOMPILE_GPS_L1CA_STATES, make_data, GPSL1CA()))
precompile_galileo(system, make_data = identity) = PositionVelocityTime.SignalGroup(
    system,
    _precompile_states(system, _PRECOMPILE_GALILEO_E1B_STATES, make_data, GalileoE1B()))
first_satellites(group, n) =
    PositionVelocityTime.SignalGroup(group.signal, group.satellites[1:n])

# Every field of two solutions, containers by content and in order.
function same_solution(a::PVTSolution, b::PVTSolution)
    a.position == b.position && a.velocity == b.velocity &&
        a.course_over_ground == b.course_over_ground &&
        a.time_correction == b.time_correction && a.time == b.time &&
        a.relative_clock_drift == b.relative_clock_drift && a.dop == b.dop &&
        collect(pairs(a.sats)) == collect(pairs(b.sats)) &&
        a.reference_system === b.reference_system &&
        a.inter_system_biases == b.inter_system_biases &&
        a.inter_frequency_biases == b.inter_frequency_biases
end

snapshot(pvt::PVTSolution) = deepcopy(pvt)

# Measured inside a function, where the keyword call is compiled like any other call —
# `@allocated` at top level would count the boxing of its own global arguments.
allocated_solve(solution, workspace, groups, prev_pvt; kw...) =
    @allocated calc_pvt!(solution, workspace, groups, prev_pvt; kw...)

@testset "calc_pvt! (allocation-free solve)" begin
    kw = (; approximate_year = 2021)
    l1ca = precompile_gps(GPSL1CA())
    e1b = precompile_galileo(GalileoE1B())
    epochs = [
        "GPS L1 C/A" => l1ca,
        "GPS L5" => precompile_gps(GPSL5I(), _precompile_cnav),
        "GPS L2C" => precompile_gps(GPSL2CM(), _precompile_cnav),
        "GPS L1C" => precompile_gps(GPSL1C_D(), _precompile_cnav2),
        "Galileo E1B" => e1b,
        "Galileo E5a" => precompile_galileo(GalileoE5aI(), _precompile_fnav),
        "GPS + Galileo" => (gps = l1ca, galileo = e1b),
        # The same satellites on two bands: an inter-frequency-bias column.
        "GPS L1 + L5" => (l1 = l1ca, l5 = precompile_gps(GPSL5I(), _precompile_cnav)),
        # Too few satellites for two independent clocks.
        "scarce GPS + Galileo" =>
            (gps = first_satellites(l1ca, 3), galileo = first_satellites(e1b, 2)),
        # The test suite's own fixtures, Dictionary-backed.
        "fixture GPS + Galileo" =>
            (gps = signal_group(gps_l1_states(0.0Hz)),
                galileo = signal_group(galileo_e1b_states(0.0Hz))),
    ]
    corrections = [
        (enable_ionospheric_correction = true, enable_tropospheric_correction = true),
        (enable_ionospheric_correction = false, enable_tropospheric_correction = false),
    ]

    @testset "the same solve as calc_pvt: $name" for (name, groups) in epochs
        for correction in corrections
            cold = calc_pvt(groups; kw..., correction...)
            warm = calc_pvt(groups, cold; kw..., correction...)
            workspace = PVTWorkspace()
            solution = PVTSolution()
            calc_pvt!(solution, workspace, groups, PVTSolution(); kw..., correction...)
            @test same_solution(solution, cold)
            # Warm, with the solution as its own seed — the receiver loop.
            calc_pvt!(solution, workspace, groups, solution; kw..., correction...)
            @test same_solution(solution, warm)
        end
    end

    @testset "overwrites the solution it is given, and only that" begin
        groups = (gps = l1ca, galileo = e1b)
        prev_pvt = calc_pvt(groups; kw...)
        before = snapshot(prev_pvt)
        solution = calc_pvt(l1ca; kw...)          # something else, to be overwritten
        returned = calc_pvt!(solution, PVTWorkspace(), groups, prev_pvt; kw...)
        @test returned === solution
        @test same_solution(solution, calc_pvt(groups, prev_pvt; kw...))
        @test same_solution(prev_pvt, before)      # `prev_pvt` is only read
        # The containers are refilled, not replaced.
        sats = solution.sats
        calc_pvt!(solution, PVTWorkspace(), l1ca, prev_pvt; kw...)
        @test solution.sats === sats
        @test length(solution.sats) == length(l1ca.satellites)
        # `calc_pvt` still never touches its `prev_pvt`.
        calc_pvt(l1ca, prev_pvt; kw...)
        @test same_solution(prev_pvt, before)
    end

    @testset "an unsolvable epoch leaves a copy of prev_pvt" begin
        prev_pvt = calc_pvt(l1ca; kw...)
        unsolvable = first_satellites(l1ca, 3)
        @test calc_pvt(unsolvable, prev_pvt; kw...) === prev_pvt
        solution = calc_pvt(e1b; kw...)
        calc_pvt!(solution, PVTWorkspace(), unsolvable, prev_pvt; kw...)
        @test same_solution(solution, prev_pvt)
        @test solution.sats !== prev_pvt.sats
        # Aliased, it is left as it was.
        before = snapshot(prev_pvt)
        calc_pvt!(prev_pvt, PVTWorkspace(), unsolvable, prev_pvt; kw...)
        @test same_solution(prev_pvt, before)
    end

    @testset "a (signal, PRN) pair twice in one epoch is refused" begin
        # Two groups sharing a ranging signal, as `Dictionary(keys, values)` refused it
        # before the solution was filled in place.
        twice = (a = l1ca, b = l1ca)
        @test_throws Dictionaries.IndexError calc_pvt(twice; kw...)
        @test_throws Dictionaries.IndexError calc_pvt!(
            PVTSolution(), PVTWorkspace(), twice, PVTSolution(); kw...)
    end

    @testset "one workspace serves epochs of any size" begin
        workspace = PVTWorkspace()
        solution = PVTSolution()
        # Growing, shrinking and growing again, across layouts with and without an IFB
        # column: every solve matches a fresh `calc_pvt`.
        for (_, groups) in epochs[[1, 9, 7, 8, 5, 10, 1]]
            calc_pvt!(solution, workspace, groups, PVTSolution(); kw...)
            @test same_solution(solution, calc_pvt(groups; kw...))
        end
    end

    @testset "allocates nothing once warm: $name" for (name, groups) in epochs
        for correction in corrections
            workspace = PVTWorkspace()
            solution = PVTSolution()
            cold = PVTSolution()
            calc_pvt!(solution, workspace, groups, cold; kw..., correction...)
            calc_pvt!(solution, workspace, groups, solution; kw..., correction...)
            @test allocated_solve(solution, workspace, groups, cold; kw..., correction...) == 0
            @test allocated_solve(solution, workspace, groups, solution; kw..., correction...) ==
                  0
        end
    end

    @testset "an unsolvable epoch allocates nothing either" begin
        prev_pvt = calc_pvt(l1ca; kw...)
        workspace = PVTWorkspace()
        solution = PVTSolution()
        calc_pvt!(solution, workspace, l1ca, prev_pvt; kw...)
        unsolvable = first_satellites(l1ca, 3)
        calc_pvt!(solution, workspace, unsolvable, prev_pvt; kw...)
        @test allocated_solve(solution, workspace, unsolvable, prev_pvt; kw...) == 0
    end

    @testset "every ionospheric model is predicted without allocating" begin
        # None of the fixtures above carries broadcast coefficients, so each model is
        # driven directly: the in-place prediction must match the allocating one, and
        # a solve through each must allocate nothing.
        models = (
            PositionVelocityTime.KlobucharParams(
                1.118e-8, 7.451e-9, -5.961e-8, -5.961e-8, 90112.0, 0.0, -196608.0,
                -65536.0),
            PositionVelocityTime.BeiDouKlobucharParams(
                1.118e-8, 7.451e-9, -5.961e-8, -5.961e-8, 90112.0, 0.0, -196608.0,
                -65536.0),
            PositionVelocityTime.NTCMGParams(236.831641, -0.39362878, 0.00402826613, 1130),
            PositionVelocityTime.BDGIMParams(
                5.25, 3.0, -0.75, 1.5, -2.0, 0.5, 0.25, -0.125, 0.5, 800),
            nothing,
        )
        rows = measurement_rows((gps = l1ca, galileo = e1b))
        ξ = [4.0186e6, 427035.0, 4.918e6, 0.0, 0.0]
        predict! = PositionVelocityTime.predict_atmospheric_delays!
        allocated_predict(delays, model) =
            @allocated predict!(delays, ξ, rows, model, 259200.0, 151, true)
        workspace = PVTWorkspace()
        solution = PVTSolution()
        calc_pvt!(solution, workspace, (gps = l1ca, galileo = e1b), PVTSolution(); kw...)
        cold = PVTSolution()
        solve!(model) = PositionVelocityTime._solve_pvt!(
            solution, workspace, PositionVelocityTime.IonosphericModel(model), cold, true)
        allocated_solve!(model) = @allocated solve!(model)
        for correction in models
            model = PositionVelocityTime.IonosphericModel(correction)
            delays = Float64[]
            predict!(delays, ξ, rows, model, 259200.0, 151, true)
            @test delays == PositionVelocityTime.predict_atmospheric_delays(
                ξ, rows, correction, 259200.0, 151, true)
            @test all(isfinite, delays)
            @test allocated_predict(delays, model) == 0
            @test solve!(correction)
            @test allocated_solve!(correction) == 0
        end
    end

    @testset "a Dictionary is emptied in place, keeping its capacity" begin
        # `empty_keeping_capacity!` reaches into `Dictionaries.Indices`; this pins that
        # the dictionary it leaves behind is a working, empty one.
        empty_in_place! = PositionVelocityTime.empty_keeping_capacity!
        dict = Dictionary{Tuple{Symbol,Int},Int}()
        fill_dict!(dict, n) = (for i in 1:n
            insert!(dict, (:GPSL1CA, i), 10i)
        end; dict)
        fill_dict!(dict, 40)
        empty_in_place!(dict)
        @test isempty(dict)
        @test !haskey(dict, (:GPSL1CA, 1))
        fill_dict!(dict, 25)
        @test collect(keys(dict)) == [(:GPSL1CA, i) for i in 1:25]
        @test all(dict[(:GPSL1CA, i)] == 10i for i in 1:25)
        @test !haskey(dict, (:GPSL1CA, 26))
        refill!(dict) = @allocated (empty_in_place!(dict); fill_dict!(dict, 30))
        refill!(dict)
        @test refill!(dict) == 0
        # Deletions leave holes; an in-place empty clears those too.
        delete!(dict, (:GPSL1CA, 3))
        empty_in_place!(dict)
        fill_dict!(dict, 5)
        @test collect(keys(dict)) == [(:GPSL1CA, i) for i in 1:5]
    end
end
