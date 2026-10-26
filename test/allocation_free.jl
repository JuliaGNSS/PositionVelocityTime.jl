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
        a.clock_bias == b.clock_bias && a.time == b.time &&
        a.clock_drift == b.clock_drift && a.dop == b.dop &&
        collect(pairs(a.satellites)) == collect(pairs(b.satellites)) &&
        a.reference_system === b.reference_system &&
        a.inter_system_biases == b.inter_system_biases &&
        a.inter_frequency_biases == b.inter_frequency_biases &&
        a.applied_broadcast_inter_system_biases == b.applied_broadcast_inter_system_biases
end

# `copy` copies a solution's containers and shares everything else, which keeps the
# copy independent only while every other field is immutable, and every container
# holds immutable entries. `Symbol` is a mutable type to Julia, but interned and never
# changed.
immutable_value(T::Union) = immutable_value(T.a) && immutable_value(T.b)
immutable_value(T::DataType) =
    T === Symbol || !ismutabletype(T) && all(immutable_value, fieldtypes(T))
container_of_immutables(T) =
    T <: Union{Dictionary,Dict} && immutable_value(keytype(T)) &&
    immutable_value(valtype(T)) || T <: Vector && immutable_value(eltype(T))

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
        # Constrained: one Galileo satellite carrying a GGTO, which ties its clock to
        # the GPS one through an extra observation row.
        "constrained GPS + Galileo" => (gps = first_satellites(l1ca, 3),
            galileo = first_satellites(precompile_galileo(GalileoE1B(), data ->
                    GNSSDecoder.GalileoINAVData(data; A_0G = 5.0e-9, A_1G = 0.0,
                        t_0G = 0, WN_0G = 0)), 1)),
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
        @test returned.satellites === solution.satellites
        @test returned.inter_system_biases === solution.inter_system_biases
        @test returned.inter_frequency_biases === solution.inter_frequency_biases
        @test returned.applied_broadcast_inter_system_biases ===
              solution.applied_broadcast_inter_system_biases
        again = calc_pvt!(returned, PVTWorkspace(), l1ca, prev_pvt; kw...)
        @test again.satellites === solution.satellites
        @test length(again.satellites) == length(l1ca.satellites)
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
            @test copied.satellites === solution.satellites
            @test same_solution(prev_pvt, before)
            # As its own output: `prev_pvt` itself, as it was.
            @test calc_pvt!(prev_pvt, PVTWorkspace(), unsolvable, prev_pvt; kw...) ===
                  prev_pvt
            @test same_solution(prev_pvt, before)
        end
    end

    @testset "copy keeps a solution the loop reuses" begin
        # With inter-system and inter-frequency biases, so every container is filled.
        for pvt in (calc_pvt((gps = l1ca, galileo = e1b); kw...),
            calc_pvt(epochs[8][2]; kw...))
            kept = copy(pvt)
            @test same_solution(kept, pvt)
            # No container shared — including any added to `PVTSolution` later.
            for name in fieldnames(PVTSolution)
                value = getfield(pvt, name)
                value isa Union{Dictionary,Dict,Vector} &&
                    @test getfield(kept, name) !== value
            end
            before = snapshot(pvt)
            calc_pvt!(pvt, PVTWorkspace(), e1b, pvt; kw...)   # refills `pvt`'s containers
            @test same_solution(kept, before)
        end
        @test same_solution(copy(PVTSolution()), PVTSolution())
    end

    @testset "every PVTSolution field is copied or immutable: $name" for (name, T) in
        zip(fieldnames(PVTSolution), fieldtypes(PVTSolution))
        @test container_of_immutables(T) || immutable_value(T)
    end

    @testset "a (signal, PRN) pair twice in one epoch is refused" begin
        # Two groups sharing a ranging signal, as `Dictionary(keys, values)` refused it
        # before the solution was filled in place.
        twice = (a = l1ca, b = l1ca)
        @test_throws Dictionaries.IndexError calc_pvt(twice; kw...)
        @test_throws Dictionaries.IndexError calc_pvt!(
            PVTSolution(), PVTWorkspace(), twice, PVTSolution(); kw...)
        # As its own output, the refused epoch leaves `prev_pvt` as it was.
        prev_pvt = calc_pvt((gps = l1ca, galileo = e1b); kw...)
        before = snapshot(prev_pvt)
        @test_throws Dictionaries.IndexError calc_pvt!(
            prev_pvt, PVTWorkspace(), twice, prev_pvt; kw...)
        @test same_solution(prev_pvt, before)
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
        parameter_sets = (
            GNSSDecoder.KlobucharParameters(
                1.118e-8, 7.451e-9, -5.961e-8, -5.961e-8, 90112.0, 0.0, -196608.0,
                -65536.0),
            GNSSDecoder.BeiDouKlobucharParameters(
                1.118e-8, 7.451e-9, -5.961e-8, -5.961e-8, 90112.0, 0.0, -196608.0,
                -65536.0),
            GNSSDecoder.NeQuickGParameters(236.831641, -0.39362878, 0.00402826613),
            GNSSDecoder.BDGIMParameters(
                5.25, 3.0, -0.75, 1.5, -2.0, 0.5, 0.25, -0.125, 0.5),
            nothing,
        )
        measurements = measurements_of((gps = l1ca, galileo = e1b))
        ξ = [4.0186e6, 427035.0, 4.918e6, 0.0, 0.0]
        predict! = PositionVelocityTime.predict_atmospheric_delays!
        epoch = get_week_start_epoch(GPST(), 2160) + 259200.0
        allocated_predict(delays, model) =
            @allocated predict!(delays, ξ, measurements, model, epoch, true)
        workspace = PVTWorkspace()
        solution = PVTSolution()
        calc_pvt!(solution, workspace, (gps = l1ca, galileo = e1b), PVTSolution(); kw...)
        cold = PVTSolution()
        solve!(model) = PositionVelocityTime._solve_pvt!(
            solution, workspace, PositionVelocityTime.IonosphericModel(model), cold, true)
        allocated_solve!(model) = @allocated solve!(model)
        for parameters in parameter_sets
            model = PositionVelocityTime.IonosphericModel(parameters)
            delays = Float64[]
            predict!(delays, ξ, measurements, model, epoch, true)
            @test delays == PositionVelocityTime.predict_atmospheric_delays(
                ξ, measurements, model, epoch, true)
            @test all(isfinite, delays)
            @test allocated_predict(delays, model) == 0
            @test first(solve!(parameters))
            @test allocated_solve!(parameters) == 0
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

    @testset "a solution's containers are emptied in place" begin
        solution = calc_pvt!(PVTSolution(), PVTWorkspace(),
            (l1 = l1ca, l5 = precompile_gps(GPSL5I(), _precompile_cnav),
                galileo = e1b), PVTSolution(); kw...)
        @test !isempty(solution.satellites)
        @test !isempty(solution.inter_system_biases)
        @test !isempty(solution.inter_frequency_biases)
        sats = solution.satellites
        @test PositionVelocityTime.empty_keeping_capacity!(solution) === solution
        @test solution.satellites === sats
        @test isempty(solution.satellites)
        @test isempty(solution.inter_system_biases)
        @test isempty(solution.inter_frequency_biases)
        @test isempty(solution.applied_broadcast_inter_system_biases)
    end

    @testset "the measurements are collected without allocating" begin
        groups = (gps = l1ca, galileo = e1b)
        measurements = PositionVelocityTime.SatelliteMeasurement[]
        parameters = PositionVelocityTime.collect_measurements!(
            measurements, groups; approximate_year = 2021)
        expected_measurements, expected_parameters =
            PositionVelocityTime.collect_measurements(groups; approximate_year = 2021)
        @test measurements == expected_measurements
        @test parameters == expected_parameters
        collect_into!(measurements, groups) = @allocated(
            PositionVelocityTime.collect_measurements!(measurements, groups;
                approximate_year = 2021))
        collect_into!(measurements, groups)
        @test collect_into!(measurements, groups) == 0
    end

    @testset "the allocating forms delegate to the in-place ones" begin
        # `curve_fit!` on a small straight-line fit; its result is the workspace's own
        # buffers.
        LM = PositionVelocityTime.LevenbergMarquardt
        xs = collect(0.0:4.0)
        ys = 2.0 .* xs .+ 1.0
        line!(out, x, p) = (out .= p[1] .* x .+ p[2]; out)
        line_jacobian!(J, x, p) = (J[:, 1] .= x; J[:, 2] .= 1.0; J)
        workspace = LM.LMWorkspace()
        fit = LM.curve_fit!(workspace, line!, line_jacobian!, xs, ys, [0.0, 0.0])
        @test fit.converged
        @test fit.param ≈ [2.0, 1.0]
        @test fit.param === workspace.x
        # A normal matrix Cholesky cannot factor (here: NaN) is a rejected step — the
        # damping grows and the parameters stay put — rather than a thrown error.
        nan_jacobian!(J, x, p) = fill!(J, NaN)
        stuck = LM.curve_fit!(
            LM.LMWorkspace(), line!, nan_jacobian!, xs, ys, [0.5, 0.5]; maxIter = 3)
        @test !stuck.converged
        @test stuck.param == [0.5, 0.5]

        # The time-system and ENU helpers.
        @test PositionVelocityTime.unique_time_systems!(
            PositionVelocityTime.SupportedTimeSystem[GST()],
            (GPST(), GST(), GPST(), BDT())) == [GPST(), GST(), BDT()]
        mixed_measurements = measurements_of((gps = l1ca, galileo = e1b))
        mixed_layout = PositionVelocityTime.decide_bias_layout(mixed_measurements)
        select_reference(layout) =
            @allocated PositionVelocityTime.select_reference_system(layout)
        select_reference(mixed_layout)
        @test select_reference(mixed_layout) == 0
        # Columns without constraints share one empty vector rather than allocate one.
        # Stored, so the compiler cannot drop the unused columns and their vector.
        indices = [1, 1]
        sink = Ref{PositionVelocityTime.BiasColumns}()
        columns(sink, indices) =
            @allocated sink[] = PositionVelocityTime.BiasColumns(indices, 1, indices, 0)
        columns(sink, indices)
        @test columns(sink, indices) == 0
        @test isempty(PositionVelocityTime.BiasColumns(indices, 1, indices, 0).constraints)
        user = ECEF(4.0186e6, 427035.0, 4.918e6)
        sat = ECEF(1.5e7, 1.0e7, 1.8e7)
        @test calc_elevation_and_azimuth(user, sat) ==
              calc_elevation_and_azimuth(ENUfromECEF(user, wgs84), sat)

        # Anything but one of the four coefficient sets (or `nothing`) is refused.
        @test_throws MethodError PositionVelocityTime.IonosphericModel(1.0)
    end
end
