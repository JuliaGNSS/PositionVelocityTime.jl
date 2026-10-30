# `calc_pvt` and `calc_pvt!`: the collection pass, then the solve behind it.

# PDOP above which a previous solution is not trusted as a warm-start seed. A PDOP of 20 is
# already an unusable fix, while the spurious far-away roots this catches show hundreds
# (the observed incident: 587); a genuine high-DOP fix discarded only costs a cold start.
const MAX_TRUSTED_WARM_START_PDOP = 50.0

"""
    PVTWorkspace()

The reusable buffers of [`calc_pvt!`](@ref), in the order the solve uses them. Each grows
to the largest epoch seen; their contents are not part of any interface.

The least-squares state, whose length depends on the epoch's systems and bands, is a plain
vector rather than a static one, so the solver compiles once rather than once per layout;
the matrices sized by it are kept at their largest and used through views.
"""
mutable struct PVTWorkspace
    measurements::Vector{PreparedObservation}
    layout::BiasLayout
    pseudoranges::Vector{Float64}
    prev_ξ::Vector{Float64}
    sat_positions::Vector{SVector{3,Float64}}
    atmospheric_delays::Vector{Float64}
    corrected_ranges::Vector{Float64}
    least_squares::LMWorkspace
    residuals::Vector{Float64}
    design_matrix::Matrix{Float64}
    normal_matrix::Matrix{Float64}
    velocity_design::Vector{SVector{4,Float64}}
    rate_residuals::Vector{Float64}
end

PVTWorkspace() = PVTWorkspace(
    PreparedObservation[],
    BiasLayout(),
    Float64[],
    Float64[],
    SVector{3,Float64}[],
    Float64[],
    Float64[],
    LMWorkspace(),
    Float64[],
    Matrix{Float64}(undef, 0, 0),
    Matrix{Float64}(undef, 0, 0),
    SVector{4,Float64}[],
    Float64[],
)

"""
    calc_pvt(groups, prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate position, velocity and time from GNSS satellite measurements.

Satellites of several GNSS and frequency bands may be combined. One receiver clock bias is
estimated per time system and one inter-frequency bias per band beyond a reference band,
so the state is `[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]`; where the satellites cannot
separate two clocks, a broadcast time offset ties them (see
[`decide_bias_layout`](@ref), which also gives the satellites a fix needs). Position and
time are found by least squares, velocity and clock drift from carrier Doppler.

The ionospheric delay is corrected from the navigation messages' own coefficients, with one
model for the whole solve (see [`select_ionospheric_parameters`](@ref)), and the
tropospheric delay with a blind model (see [`tropospheric_delay`](@ref)).

# Arguments
- `groups`: this epoch's measurements, as signal groups: a `NamedTuple` of
  [`SignalGroup`](@ref)s, one per ranging signal. A bare `Tuple` of groups (numbered
  `group1`, `group2`, …) or a single `SignalGroup` (the `:default` group) also does:

  ```julia
  using PositionVelocityTime: SignalGroup
  calc_pvt(SignalGroup(GPSL1CA(), gps_states))
  calc_pvt((gps = SignalGroup(GPSL1CA(), gps_states),
            galileo = SignalGroup(GalileoE1B(), galileo_states)))
  ```

  Grouping keeps the solve type-stable (see [`SignalGroup`](@ref)); a pooled
  `Vector{SignalObservation}` is refused, and [`signal_groups`](@ref) builds the groups from a
  Tracking.jl `TrackState`.

  Each `(signal, PRN)` pair must appear at most once, or it enters the fit twice; two
  groups sharing a ranging signal are not checked for collisions.

  **Order is significant**: the measurements are in group order × within-group order,
  which breaks the ties of the reference system ("most satellites, first appearance") and
  of [`ifb_layout`](@ref)'s reference bands, and orders
  `PVTSolution.observations`.
- `prev_pvt`: the previous solution, the warm start (default: the origin). One whose DOP
  is implausible (`NaN`, or `PDOP` above `MAX_TRUSTED_WARM_START_PDOP` = 50) is not
  used, so a spurious fix cannot re-seed itself.

# Keyword Arguments
- `approximate_year`: the calendar year of the observation, within ±9 years, which
  resolves the GPS L1 C/A 1024-week rollover (see GNSSDecoder's `get_broadcast_week`). The default,
  the current year, suits live signals; pass the rough year of an archived recording.
- `enable_ionospheric_correction`: `false` skips the ionospheric correction.
- `enable_tropospheric_correction`: `false` skips the tropospheric correction.

# Returns
A new [`PVTSolution`](@ref); `prev_pvt` is never modified. An epoch that cannot be solved
— too few healthy satellites (see [`decide_bias_layout`](@ref)) or a rank-deficient
geometry — returns `prev_pvt` as passed, without throwing, so a receiver can hand over
whatever it tracks and carry the last fix forward.
"""
function calc_pvt(
    groups,
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    # A fresh solution lends the solve its containers.
    _, solution = _calc_pvt!(PVTSolution(), PVTWorkspace(), groups, prev_pvt,
        approximate_year, enable_ionospheric_correction, enable_tropospheric_correction)
    solution
end

"""
    calc_pvt!(solution::PVTSolution, workspace::PVTWorkspace, groups, prev_pvt::PVTSolution;
              approximate_year::Integer = year(now(UTC)),
              enable_ionospheric_correction::Bool = true,
              enable_tropospheric_correction::Bool = true) -> PVTSolution

[`calc_pvt`](@ref) — the same solve and result — reusing the containers of `solution`
and the buffers of `workspace`, so a stream of epochs is solved without allocating. The
fix is a new, immutable [`PVTSolution`](@ref); keep the return value:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for groups in epochs
    pvt = calc_pvt!(pvt, workspace, groups, pvt)
end
```

**`solution` is the output, and the only thing overwritten**: the returned solution takes
over its four containers, emptied and refilled. Afterwards `solution`'s scalar fields still
describe its own epoch but its containers the new one, so keep a
[`copy`](@ref Base.copy(::PVTSolution)) of it beforehand if you need it (a `Dictionary` made by `map` over
`solution.observations` shares its keys and sees the change too). `prev_pvt` is only read,
and may be `solution` itself, as in the loop above.

An epoch that cannot be solved returns `prev_pvt`: itself where it shares `solution`'s
containers, otherwise a copy of it in them.

It allocates nothing once `workspace` and `solution` have held an epoch of as many
satellites and biases; that covers everything but the `groups`, which the caller builds
(and can reuse: a `SignalGroup`'s satellites may be any vector). The `workspace` carries no
state between epochs, but serves one solve at a time; give each task its own.
"""
function calc_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    groups,
    prev_pvt::PVTSolution;
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    solved, result = _calc_pvt!(solution, workspace, groups, prev_pvt, approximate_year,
        enable_ionospheric_correction, enable_tropospheric_correction)
    solved ? result : copy_into(solution, prev_pvt)
end

# The collection pass and the solve: `(true, fix)` in `solution`'s containers, or
# `(false, prev_pvt)` with `solution` untouched. Specialises on the group shape.
function _calc_pvt!(
    solution,
    workspace,
    groups,
    prev_pvt,
    approximate_year,
    enable_ionospheric_correction,
    enable_tropospheric_correction,
)
    # The function barrier (see collection.jl): the measurements in, one compiled solver
    # behind.
    ionospheric_parameters =
        _collect_measurements!(workspace.measurements, groups, approximate_year)
    _solve_pvt!(
        solution,
        workspace,
        IonosphericModel(enable_ionospheric_correction ? ionospheric_parameters : nothing),
        prev_pvt,
        enable_tropospheric_correction,
    )
end

# The solver, compiled once: every argument is concretely typed, the ionospheric model
# wrapped in an `IonosphericModel`. Returns `(true, fix)` in `solution`'s containers or
# `(false, prev_pvt)`; the containers are written only once the epoch is solvable and
# everything is read from `prev_pvt`, which may be the same solution. A tuple, as a
# `Union{Nothing,PVTSolution}` would be returned boxed.
function _solve_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    ionospheric_model,
    prev_pvt::PVTSolution,
    enable_tropospheric_correction::Bool,
)
    measurements = workspace.measurements
    num_sats = length(measurements)

    # Solvability by the counts; a degenerate geometry is caught after the solve by the DOP.
    layout = workspace.layout
    decide_bias_layout!(layout, measurements) || return (false, prev_pvt)
    (; clock_systems, extra_bands, reference_bands, constraints) = layout
    bias_columns = fit_bias_columns(layout)
    (; clock_columns, num_clock_biases) = bias_columns
    num_params = num_parameters(bias_columns)

    # The reference system dates the reported time (see `select_reference_system`).
    reference_clock_column = select_reference_clock_column(layout)
    reference_system = clock_systems[reference_clock_column]
    # `something` keeps the index an `Int`: every clock column has a satellite.
    reference_index = something(findfirst(==(reference_clock_column), clock_columns))

    # The pseudoranges: how far each transmit time trails the reference measurement's, on
    # the reference system's count (`calc_time_scale_offset`; raw, a BeiDou
    # measurement would carry 4.2e9 m) and folded across a week crossover
    # (`calc_pseudorange`). The common reference cancels into the receiver clock biases,
    # which then date `reference_time`, the reference measurement's own transmit time, in
    # its own week.
    reference_time = measurements[reference_index].transmit_time
    pseudoranges = resize!(workspace.pseudoranges, num_sats + length(constraints))
    for (j, measurement) in enumerate(measurements)
        pseudoranges[j] = calc_pseudorange(reference_time, measurement.transmit_time,
            calc_time_scale_offset(measurement, reference_system))
    end
    # The clock constraints' known values are observations too, after the satellites'.
    for (i, constraint) in enumerate(constraints)
        pseudoranges[num_sats+i] = constraint.value
    end

    # The reference week's start dates the epoch: the atmospheric delays' here, before the
    # clock correction, and the reported time at the end.
    reference_measurement = measurements[reference_index]
    week_start =
        get_week_start_epoch(reference_measurement.time_system, reference_measurement.week)
    epoch = week_start + reference_time

    # A spurious far-away root would re-seed itself epoch after epoch; it shows an
    # impossible DOP, and is not used as a seed.
    distrusted =
        !isnothing(prev_pvt.dop) &&
        (isnan(prev_pvt.dop.GDOP) || prev_pvt.dop.PDOP > MAX_TRUSTED_WARM_START_PDOP)

    # Seed each clock bias from the previous solution, as its reference bias plus its
    # inter-system bias; a system it lacks gets the reference bias, metres off rather than
    # the whole clock bias. A distrusted solution seeds nothing: the cold start, all zeros.
    prev_abs_bias(sys) =
        sys === prev_pvt.reference_system ? prev_pvt.clock_bias :
        prev_pvt.clock_bias + get(prev_pvt.inter_system_biases, sys, 0.0m)
    prev_ξ = fill!(resize!(workspace.prev_ξ, num_params), 0.0)
    if !distrusted
        prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
        for (k, sys) in enumerate(clock_systems)
            prev_ξ[3+k] = ustrip(m, prev_abs_bias(sys))
        end
        # An IFB is seeded only against the same reference band, as another reference makes
        # it a different quantity. `haskey` and an index, not `get(…, nothing)`: a `Union`
        # of an `InterFrequencyBias` (holding a `Symbol`) and `nothing` is returned boxed.
        prev_ifbs = prev_pvt.inter_frequency_biases
        for (i, band) in enumerate(extra_bands)
            prev_ξ[3+num_clock_biases+i] =
                haskey(prev_ifbs, band) && prev_ifbs[band].reference_band == reference_bands[i] ?
                ustrip(m, prev_ifbs[band].value) : 0.0
        end
    end

    sat_positions = resize!(workspace.sat_positions, num_sats)
    for j in 1:num_sats
        sat_positions[j] = measurements[j].position
    end

    # The atmospheric delays come off the pseudoranges, unless neither model applies.
    correct_atmosphere =
        !isnothing(ionospheric_model.parameters) || enable_tropospheric_correction
    function corrected_ranges(ξ_prediction)
        delays = predict_atmospheric_delays!(
            workspace.atmospheric_delays, ξ_prediction, measurements, ionospheric_model,
            epoch, enable_tropospheric_correction)
        ranges = resize!(workspace.corrected_ranges, length(pseudoranges)) .= pseudoranges
        view(ranges, 1:num_sats) .-= delays
        ranges
    end
    least_squares = workspace.least_squares

    # `ξ` is the least-squares workspace's own vector, valid until its next solve.
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: the delays need a position (Klobuchar is undefined near the
        # geocentre), so an uncorrected fix first, then one corrected solve.
        ξ_uncorrected, residuals_uncorrected = estimate_position!(least_squares,
            workspace.residuals, sat_positions, pseudoranges, bias_columns, prev_ξ)
        if correct_atmosphere
            estimate_position!(least_squares, workspace.residuals, sat_positions,
                corrected_ranges(ξ_uncorrected), bias_columns, ξ_uncorrected)
        else
            (ξ_uncorrected, residuals_uncorrected)
        end
    else
        # Warm start: the delays from the previous, metre-accurate position.
        estimate_position!(least_squares, workspace.residuals, sat_positions,
            correct_atmosphere ? corrected_ranges(prev_ξ) : pseudoranges, bias_columns,
            prev_ξ)
    end
    position = ECEF(ξ[1], ξ[2], ξ[3])

    # The geometry at the converged position: a rank deficiency (a `NaN` DOP) makes `ξ`
    # meaningless, and an epoch unsolvable. The velocity solve backs this up with its own
    # rank test (see `estimate_velocity`). The constraints count as exact
    # (`dop_bias_columns`).
    dop_columns = dop_bias_columns(layout)
    num_dop_params = num_parameters(dop_columns)
    workspace.design_matrix = grown(workspace.design_matrix, num_sats, num_dop_params)
    H = calc_design_matrix!(view(workspace.design_matrix, 1:num_sats, 1:num_dop_params),
        sat_positions, ξ, dop_columns)
    workspace.normal_matrix = grown(workspace.normal_matrix, num_dop_params, num_dop_params)
    dop = calc_dop!(view(workspace.normal_matrix, 1:num_dop_params, 1:num_dop_params), H,
        position, layout.clock_group_columns[reference_index])

    isnan(dop.GDOP) && return (false, prev_pvt)

    # The velocity solve starts from the previous solution where the position solve trusts
    # it, from rest otherwise (see `estimate_velocity` for why one step suffices).
    previous_velocity_and_drift = SVector(prev_pvt.velocity[1], prev_pvt.velocity[2],
        prev_pvt.velocity[3], ustrip(m / s, prev_pvt.clock_drift))
    velocity_start =
        !distrusted && all(isfinite, previous_velocity_and_drift) ?
        previous_velocity_and_drift : zero(SVector{4,Float64})
    user_velocity_and_clock_drift, rate_residuals =
        estimate_velocity!(workspace.rate_residuals,
            workspace.velocity_design, measurements, position, velocity_start)
    all(isfinite, user_velocity_and_clock_drift) || return (false, prev_pvt)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    clock_bias = ξ[3+reference_clock_column]
    # A receiver clock running ahead reads the epoch late by its bias.
    corrected_reference_time = reference_time - clock_bias / SPEED_OF_LIGHT

    # Everything is read from `prev_pvt` by now; `solution`'s containers are written from
    # here on. A duplicate `(signal, PRN)` is refused before any is, so a throw leaves
    # `solution` whole; not by `insert!`, whose message formats the key, a dynamic `show`
    # a `juliac --trim` build rejects.
    for j in 2:num_sats, k in 1:j-1
        measurements[j].signal_id === measurements[k].signal_id &&
            measurements[j].prn == measurements[k].prn && throw(IndexError(
                "a (signal, PRN) pair appears twice in one epoch; each must appear once"))
    end
    observations = empty_keeping_capacity!(solution.observations)
    for (j, measurement) in enumerate(measurements)
        key = (measurement.signal_id, measurement.prn)
        set!(observations, key, ObservationResult(measurement.position,
            measurement.transmit_time * s, residuals[j] * m, rate_residuals[j] * (m / s)))
    end

    # The reference system's own bias is `clock_bias`.
    inter_system_biases = empty_keeping_capacity!(solution.inter_system_biases)
    for (k, sys) in enumerate(clock_systems)
        sys === reference_system && continue
        set!(inter_system_biases, sys, (ξ[3+k] - clock_bias) * m)
    end

    # The reference band's delay is folded into the clock biases.
    inter_frequency_biases = empty_keeping_capacity!(solution.inter_frequency_biases)
    for (i, band) in enumerate(extra_bands)
        set!(inter_frequency_biases, band,
            InterFrequencyBias(ξ[3+num_clock_biases+i] * m, reference_bands[i]))
    end

    # The receiver-clock differences constrained by broadcast offsets, by time system.
    applied_biases = empty!(solution.broadcast_clock_constraints)
    for constraint in constraints
        push!(applied_biases, BroadcastClockConstraint(
            clock_systems[constraint.column], clock_systems[constraint.target_column],
            constraint.value * m))
    end

    return true,
    PVTSolution(
        position,
        velocity,
        calc_course_over_ground(position, velocity),
        week_start + corrected_reference_time,
        clock_bias * m,
        user_velocity_and_clock_drift[4] * (m / s),
        reference_system,
        inter_system_biases,
        applied_biases,
        inter_frequency_biases,
        dop,
        observations,
    )
end
