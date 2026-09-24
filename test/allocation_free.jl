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
function allocated_loop(pvt, workspace, groups; kw...)
    @allocated for _ in 1:10
        pvt = calc_pvt!(pvt, workspace, groups, pvt; kw...)
    end
end

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
            pvt = calc_pvt!(PVTSolution(), workspace, groups, PVTSolution(); kw...,
                correction...)
            @test same_solution(pvt, cold)
            # Warm, with the solution as its own seed — the receiver loop.
            pvt = calc_pvt!(pvt, workspace, groups, pvt; kw..., correction...)
            @test same_solution(pvt, warm)
        end
    end

    @testset "PVTSolution stays immutable" begin
        @test !ismutabletype(PVTSolution)
        @test_throws ErrorException calc_pvt(l1ca; kw...).position = ECEF(0.0, 0.0, 0.0)
    end

    @testset "reuses the solution it is given, and only that" begin
        groups = (gps = l1ca, galileo = e1b)
        prev_pvt = calc_pvt(groups; kw...)
        before = snapshot(prev_pvt)
        solution = calc_pvt(l1ca; kw...)          # something else, to be reused
        returned = calc_pvt!(solution, PVTWorkspace(), groups, prev_pvt; kw...)
        @test same_solution(returned, calc_pvt(groups, prev_pvt; kw...))
        @test same_solution(prev_pvt, before)      # `prev_pvt` is only read
        # The returned solution took over `solution`'s containers, refilled.
        @test returned.sats === solution.sats
        @test returned.inter_system_biases === solution.inter_system_biases
        @test returned.inter_frequency_biases === solution.inter_frequency_biases
        again = calc_pvt!(returned, PVTWorkspace(), l1ca, prev_pvt; kw...)
        @test again.sats === solution.sats
        @test length(again.sats) == length(l1ca.satellites)
        # `calc_pvt` still never touches its `prev_pvt`.
        calc_pvt(l1ca, prev_pvt; kw...)
        @test same_solution(prev_pvt, before)
    end

    @testset "an unsolvable epoch returns prev_pvt" begin
        unsolvable = first_satellites(l1ca, 3)
        # With inter-system and inter-frequency biases, so every container is copied.
        for prev_pvt in (calc_pvt(l1ca; kw...), calc_pvt(epochs[8][2]; kw...),
            calc_pvt((gps = l1ca, galileo = e1b); kw...))
            @test calc_pvt(unsolvable, prev_pvt; kw...) === prev_pvt
            before = snapshot(prev_pvt)
            # Into another solution's containers: a copy, and `prev_pvt` untouched.
            solution = calc_pvt(e1b; kw...)
            copied = calc_pvt!(solution, PVTWorkspace(), unsolvable, prev_pvt; kw...)
            @test same_solution(copied, prev_pvt)
            @test copied.sats === solution.sats
            @test same_solution(prev_pvt, before)
            # As its own output: `prev_pvt` itself, as it was.
            @test calc_pvt!(prev_pvt, PVTWorkspace(), unsolvable, prev_pvt; kw...) ===
                  prev_pvt
            @test same_solution(prev_pvt, before)
        end
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
        pvt = PVTSolution()
        # Growing, shrinking and growing again, across layouts with and without an IFB
        # column: every solve matches a fresh `calc_pvt`.
        for (_, groups) in epochs[[1, 9, 7, 8, 5, 10, 1]]
            pvt = calc_pvt!(pvt, workspace, groups, PVTSolution(); kw...)
            @test same_solution(pvt, calc_pvt(groups; kw...))
        end
    end

    @testset "allocates nothing once warm: $name" for (name, groups) in epochs
        for correction in corrections
            workspace = PVTWorkspace()
            cold = PVTSolution()
            pvt = calc_pvt!(PVTSolution(), workspace, groups, cold; kw..., correction...)
            pvt = calc_pvt!(pvt, workspace, groups, pvt; kw..., correction...)
            @test allocated_solve(pvt, workspace, groups, cold; kw..., correction...) == 0
            @test allocated_solve(pvt, workspace, groups, pvt; kw..., correction...) == 0
            # The receiver loop itself, many epochs through one binding.
            @test allocated_loop(pvt, workspace, groups; kw..., correction...) == 0
        end
    end

    @testset "an unsolvable epoch allocates nothing either" begin
        prev_pvt = calc_pvt(epochs[8][2]; kw...)
        workspace = PVTWorkspace()
        solution = calc_pvt!(PVTSolution(), workspace, epochs[8][2], prev_pvt; kw...)
        unsolvable = first_satellites(l1ca, 3)
        calc_pvt!(solution, workspace, unsolvable, prev_pvt; kw...)
        @test allocated_solve(solution, workspace, unsolvable, prev_pvt; kw...) == 0
        @test allocated_solve(prev_pvt, workspace, unsolvable, prev_pvt; kw...) == 0
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
            @test first(solve!(correction))
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

    @testset "the allocating forms delegate to the in-place ones" begin
        # `curve_fit` is `curve_fit!` on a fresh workspace: a small straight-line fit.
        LM = PositionVelocityTime.LevenbergMarquardt
        xs = collect(0.0:4.0)
        ys = 2.0 .* xs .+ 1.0
        line!(out, x, p) = (out .= p[1] .* x .+ p[2]; out)
        line_jacobian!(J, x, p) = (J[:, 1] .= x; J[:, 2] .= 1.0; J)
        fit = LM.curve_fit(line!, line_jacobian!, xs, ys, [0.0, 0.0])
        @test fit.converged
        @test fit.param ≈ [2.0, 1.0]
        workspace = LM.LMWorkspace()
        in_place = LM.curve_fit!(workspace, line!, line_jacobian!, xs, ys, [0.0, 0.0])
        @test in_place.param == fit.param
        @test in_place.param === workspace.x
        # A normal matrix Cholesky cannot factor (here: NaN) is a rejected step — the
        # damping grows and the parameters stay put — rather than a thrown error.
        nan_jacobian!(J, x, p) = fill!(J, NaN)
        stuck = LM.curve_fit(line!, nan_jacobian!, xs, ys, [0.5, 0.5]; maxIter = 3)
        @test !stuck.converged
        @test stuck.param == [0.5, 0.5]

        # The time-system and ENU helpers.
        @test PositionVelocityTime.unique_time_systems((GPST(), GST(), GPST(), BDT())) ==
              [GPST(), GST(), BDT()]
        user = ECEF(4.0186e6, 427035.0, 4.918e6)
        sat = ECEF(1.5e7, 1.0e7, 1.8e7)
        @test get_sat_enu(user, sat) == get_sat_enu(ENUfromECEF(user, wgs84), sat)

        # Anything but one of the four coefficient sets (or `nothing`) is refused.
        rows = measurement_rows(l1ca)
        @test_throws ArgumentError PositionVelocityTime.predict_atmospheric_delays(
            [4.0186e6, 427035.0, 4.918e6, 0.0], rows, 1.0, 259200.0, 151, true)
    end
end
