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
    prepared::Vector{PreparedObservation}
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
    calc_pvt(observations, prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> (; status, solution)

Calculate position, velocity and time from GNSS signal observations.

Satellites of several GNSS and frequency bands may be combined. One receiver clock bias is
estimated per time system and one inter-frequency bias per band beyond a reference band,
so the state is `[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]`; where the satellites cannot
separate two clocks, a broadcast time offset ties them (see
[`decide_bias_layout`](@ref), which also gives the satellites a fix needs). Position and
time are found by least squares, velocity and clock drift from carrier Doppler.

The ionospheric delay is corrected from the navigation messages' own coefficients, with one
model for the whole solve (see [`select_ionospheric_correction`](@ref)), and the
tropospheric delay with a blind model (see [`tropospheric_delay`](@ref)).

# Arguments
- `observations`: a concrete-typed vector or `Dictionary` of [`SignalObservation`](@ref)s,
  or a tuple / NamedTuple of those collections for several decoder or signal types:

  ```julia
  calc_pvt(gps_observations)
  calc_pvt((gps = gps_observations, galileo = galileo_observations))
  ```

  Collections keep the preparation pass statically dispatched. A pooled vector with
  an abstract element type is refused. NamedTuple keys are only caller labels; the
  solver uses each observation's signal and decoder PRN. Empty collections retain
  their concrete element type (for example `gps_observations[1:0]`).

  Each `(signal, PRN)` pair must appear at most once, or it enters the fit twice; two
  collections sharing a ranging signal are not checked for collisions.

  **Order is significant**: the observations are prepared in collection order ×
  within-collection order, which breaks the ties of the reference system ("most
  satellites, first appearance") and of [`ifb_layout`](@ref)'s reference bands, and orders
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
`(; status, solution)`: a [`PVTStatus`](@ref) and, when it is `PositionVelocityTime.solved`,
the new [`PVTSolution`](@ref). On failure `solution` is an empty `PVTSolution()`, not a
fix, and nothing throws. `prev_pvt` only seeds the solve and is never modified; whether to
carry an earlier fix forward is the caller's choice.
"""
function calc_pvt(
    observations,
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    # A fresh solution lends the solve its containers.
    _calc_pvt!(PVTSolution(), PVTWorkspace(), observations, prev_pvt,
        approximate_year, enable_ionospheric_correction, enable_tropospheric_correction)
end

"""
    calc_pvt!(solution::PVTSolution, workspace::PVTWorkspace, observations,
              prev_pvt::PVTSolution;
              approximate_year::Integer = year(now(UTC)),
              enable_ionospheric_correction::Bool = true,
              enable_tropospheric_correction::Bool = true) -> (; status, solution)

[`calc_pvt`](@ref) — the same solve and result — reusing the containers of `solution`
and the buffers of `workspace`, so a stream of epochs is solved without allocating. On
failure the returned `solution` is the `solution` argument, untouched. A receiver loop:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for observations in epochs
    result = calc_pvt!(pvt, workspace, observations, pvt)
    if result.status == PositionVelocityTime.solved
        pvt = result.solution
    end
end
```

**`solution` is the output, and the only thing overwritten**: the returned solution takes
over its four containers, emptied and refilled. Afterwards `solution`'s scalar fields still
describe its own epoch but its containers the new one, so keep a
[`copy`](@ref Base.copy(::PVTSolution)) of it beforehand if you need it (a `Dictionary` made by `map` over
`solution.observations` shares its keys and sees the change too). `prev_pvt` is only read,
and may be `solution` itself, as in the loop above.

It allocates nothing once `workspace` and `solution` have held an epoch of as many
satellites and biases; that covers everything but the `observations`, which the caller
builds (which can be reused vectors or dictionaries). The `workspace` carries no
state between epochs, but serves one solve at a time; give each task its own.
"""
function calc_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    observations,
    prev_pvt::PVTSolution;
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    _calc_pvt!(solution, workspace, observations, prev_pvt, approximate_year,
        enable_ionospheric_correction, enable_tropospheric_correction,
    )
end

# The collection pass and the solve, specialized on the observation collection shape.
function _calc_pvt!(
    solution,
    workspace,
    observations,
    prev_pvt,
    approximate_year,
    enable_ionospheric_correction,
    enable_tropospheric_correction,
)
    # The function barrier (see collection.jl): the prepared observations in, one compiled
    # solver behind.
    collections = observation_collections(observations)
    (; ionospheric_model, time_offsets) =
        _collect_measurements!(workspace.prepared, collections, approximate_year)
    _solve_pvt!(
        solution,
        workspace,
        enable_ionospheric_correction ? ionospheric_model : IonosphericModel(nothing),
        time_offsets,
        prev_pvt,
        enable_tropospheric_correction,
    )
end

# The solver, compiled once behind the typed preparation pass. Containers are written
# only after both geometry checks pass and everything has been read from `prev_pvt`,
# which may be the output solution itself.
function _solve_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    ionospheric_model,
    time_offsets::BroadcastTimeOffsets,
    prev_pvt::PVTSolution,
    enable_tropospheric_correction::Bool,
)
    prepared = workspace.prepared
    num_observations = length(prepared)

    # Solvability by the counts; a degenerate geometry is caught after the solve by the DOP.
    layout = workspace.layout
    decide_bias_layout!(layout, prepared, time_offsets) ||
        return (; status = insufficient_observations, solution)
    (; clock_systems, extra_bands, reference_bands, constraints) = layout
    bias_columns = fit_bias_columns(layout)
    (; clock_columns, num_clock_biases) = bias_columns
    num_params = num_parameters(bias_columns)

    # The reference system dates the reported time (see `select_reference_system`).
    reference_clock_column = select_reference_clock_column(layout)
    reference_system = clock_systems[reference_clock_column]
    # `something` keeps the index an `Int`: every clock column has a satellite.
    reference_index = something(findfirst(==(reference_clock_column), clock_columns))

    # The pseudoranges: how far each transmit time trails the reference observation's, on
    # the reference system's count (GNSSSignals' `get_time_system_offset`; raw, a BeiDou
    # observation would carry 4.2e9 m) and folded across a week crossover. A crossover can
    # fall between the times: within the transmit-time spread, or, in a mixed solve,
    # between BeiDou times on the GPS count (`SOW + 14`, past 604800 for the first 14 s of
    # a GPS week) and GPS times already wrapped. The common reference cancels into the
    # receiver clock biases, which then date `reference_time`, the reference observation's
    # own transmit time, in its own week. Only the differenced times move: the satellite is
    # propagated on its own system's scale, where the 14 s would put 55 km of along-track
    # error into the ephemeris.
    reference_time = prepared[reference_index].transmit_time
    pseudoranges = resize!(workspace.pseudoranges, num_observations)
    for (j, observation) in enumerate(prepared)
        correction =
            ustrip(s, get_time_system_offset(observation.time_system, reference_system))
        pseudoranges[j] =
            fold_week_crossover(reference_time - (correction + observation.transmit_time)) *
            SPEED_OF_LIGHT
    end
    # The clock constraints' known values are observations too, after the satellites'.
    resize!(pseudoranges, num_observations + length(constraints))
    for (i, constraint) in enumerate(constraints)
        pseudoranges[num_observations+i] = constraint.value
    end

    # The reference week's start dates the epoch: the atmospheric delays' here, before the
    # clock correction, and the reported time at the end.
    reference_observation = prepared[reference_index]
    week_start =
        get_week_start_epoch(reference_observation.time_system, reference_observation.week)
    epoch = week_start + reference_time

    # A spurious far-away root would re-seed itself epoch after epoch; it shows an
    # impossible DOP, and is not used as a seed.
    distrusted =
        !isnothing(prev_pvt.dop) &&
        (isnan(prev_pvt.dop.GDOP) || prev_pvt.dop.PDOP > MAX_TRUSTED_WARM_START_PDOP)

    # Seed the reference clock bias from the previous position: the reference observation's
    # pseudorange is zero, so its bias is minus its range, up to the atmosphere. Each other
    # clock adds its previous inter-system bias, taken against the new reference; a system
    # the previous solution lacks gets the reference bias, metres off rather than the whole
    # clock bias. A distrusted solution, or none (the origin), seeds nothing: the cold
    # start, all zeros.
    prev_inter_system_bias(sys) = ustrip(m, get(prev_pvt.inter_system_biases, sys, 0.0m))
    prev_ξ = fill!(resize!(workspace.prev_ξ, num_params), 0.0)
    if !distrusted
        prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
        if !iszero(prev_pvt.position)
            reference_bias = -norm(prev_pvt.position - reference_observation.position) -
                             prev_inter_system_bias(reference_system)
            for (k, sys) in enumerate(clock_systems)
                prev_ξ[3+k] = reference_bias + prev_inter_system_bias(sys)
            end
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

    sat_positions = resize!(workspace.sat_positions, num_observations)
    for j in 1:num_observations
        sat_positions[j] = prepared[j].position
    end

    # The atmospheric delays come off the pseudoranges, unless neither model applies.
    correct_atmosphere =
        !isnothing(ionospheric_model.parameters) || enable_tropospheric_correction
    function corrected_ranges(ξ_prediction)
        delays = predict_atmospheric_delays!(
            workspace.atmospheric_delays, ξ_prediction, prepared, ionospheric_model,
            epoch, enable_tropospheric_correction)
        ranges = resize!(workspace.corrected_ranges, length(pseudoranges)) .= pseudoranges
        view(ranges, 1:num_observations) .-= delays
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
    workspace.design_matrix =
        grown(workspace.design_matrix, num_observations, num_dop_params)
    H = calc_design_matrix!(
        view(workspace.design_matrix, 1:num_observations, 1:num_dop_params),
        sat_positions, ξ, dop_columns)
    workspace.normal_matrix = grown(workspace.normal_matrix, num_dop_params, num_dop_params)
    dop = calc_dop!(view(workspace.normal_matrix, 1:num_dop_params, 1:num_dop_params), H,
        position, layout.clock_group_columns[reference_index])

    isnan(dop.GDOP) && return (; status = degenerate_geometry, solution)

    # The velocity solve starts from the previous solution where the position solve trusts
    # it, from rest otherwise (see `estimate_velocity` for why one step suffices).
    previous_velocity_and_drift = SVector(prev_pvt.velocity[1], prev_pvt.velocity[2],
        prev_pvt.velocity[3], ustrip(m / s, prev_pvt.clock_drift))
    velocity_start =
        !distrusted && all(isfinite, previous_velocity_and_drift) ?
        previous_velocity_and_drift : zero(SVector{4,Float64})
    user_velocity_and_clock_drift, rate_residuals =
        estimate_velocity!(workspace.rate_residuals,
            workspace.velocity_design, prepared, position, velocity_start)
    all(isfinite, user_velocity_and_clock_drift) ||
        return (; status = degenerate_geometry, solution)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    # The reference clock bias counts from the reference observation's transmit time, so
    # it is minus the signal's travel time, not a receiver clock offset: it dates the
    # epoch and is not reported.
    clock_bias = ξ[3+reference_clock_column]
    corrected_reference_time = reference_time - clock_bias / SPEED_OF_LIGHT

    # Everything is read from `prev_pvt` by now; `solution`'s containers are written from
    # here on. A duplicate `(signal, PRN)` is refused before any is, so a throw leaves
    # `solution` whole; not by `insert!`, whose message formats the key, a dynamic `show`
    # a `juliac --trim` build rejects.
    for j in 2:num_observations, k in 1:j-1
        prepared[j].signal_id === prepared[k].signal_id &&
            prepared[j].prn == prepared[k].prn && throw(IndexError(
                "a (signal, PRN) pair appears twice in one epoch; each must appear once"))
    end
    observations = empty_keeping_capacity!(solution.observations)
    for (j, observation) in enumerate(prepared)
        key = (observation.signal_id, observation.prn)
        set!(observations, key, ObservationResiduals(observation.position,
            observation.transmit_time * s, residuals[j] * m, rate_residuals[j] * (m / s)))
    end

    # Each other clock against the reference one.
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
    applied_biases = empty!(solution.applied_broadcast_inter_system_biases)
    for constraint in constraints
        push!(applied_biases, BroadcastInterSystemBias(
            clock_systems[constraint.column], clock_systems[constraint.target_column],
            constraint.value * m))
    end

    fix = PVTSolution(
        position,
        velocity,
        calc_course_over_ground(position, velocity),
        week_start + corrected_reference_time,
        user_velocity_and_clock_drift[4] * (m / s),
        reference_system,
        inter_system_biases,
        applied_biases,
        inter_frequency_biases,
        dop,
        observations,
    )
    (; status = solved, solution = fix)
end
