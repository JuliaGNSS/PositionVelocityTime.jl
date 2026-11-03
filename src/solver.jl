# `calc_pvt` and `calc_pvt!`: the collection pass, then the solve behind it.

# PDOP above which a previous solution is distrusted as a warm-start seed and
# discarded (see the gate at the top of `calc_pvt`). Genuine geometries this
# side of useless stay well under it — a PDOP of 20 is already an unusable fix —
# while the spurious far-away roots the gate exists to catch show hundreds (the
# observed incident: 587). Discarding a fix that honestly earned a high DOP is
# harmless — the cold solve lands in the same place — so the exact value only
# sets how often that happens.
const MAX_TRUSTED_WARM_START_PDOP = 50.0

"""
    PVTWorkspace()

The reusable scratch storage of [`calc_pvt!`](@ref): the observations, the bias layout,
the pseudoranges and their corrections, the least-squares buffers, the design and
normal-equations matrices, the velocity solve's design rows and the residuals. Every
buffer grows to the largest epoch it has seen and is reused from then on; its contents
after a solve are not part of any interface.

Satellite positions are kept as a `Vector{SVector{3,Float64}}` — a fixed-size quantity
per satellite, in a vector that is resized in place. The least-squares state, whose
length `3 + M + B` depends on the epoch's time systems and bands, is a plain vector
instead of a static one, so the solver compiles once rather than once per layout; the
matrices sized by it are kept at their largest and used through views.
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
    calc_pvt(groups, prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS signal observations.

Satellites from different GNSS and frequency bands may be combined. Because each
constellation references its broadcasts to its own system time, one receiver clock
bias is estimated per GNSS time system; and because the receiver's RF chain delays
each band differently, one receiver inter-frequency bias is estimated per frequency
band beyond a reference band (shared across constellations on that band; see
`get_band_id`). The state vector is therefore
`[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]` for `M` distinct time systems and `B`
extra bands. Position and time are found by least squares; velocity and clock drift
are solved from carrier Doppler.

A solution requires `n ≥ 3 + M + B` healthy satellite observations (each system needs at
least one satellite, and a system contributing a single satellite spends it entirely
on that system's clock bias), of which `3 + M` must come from *distinct* satellites: a
satellite tracked on several bands supplies one line of sight, and its extra
observations constrain the inter-frequency biases rather than the geometry. When either
condition fails, constellations whose messages carry a broadcast offset to another
tracked system — Galileo's GGTO (Galileo–GPS Time Offset), BeiDou's BGTO (BDT–GNSS
Time Offset, toward GPS or Galileo), or the GGTO GPS itself broadcasts on CNAV/CNAV-2
— have their clock tied to that system's by the broadcast offset, a constraint the fit
meets exactly, which makes a 4-satellite mixed fix possible (see
[`decide_bias_layout`](@ref)). Estimating an independent bias is preferred whenever the
geometry allows it, since it avoids the broadcast offset's own error.

Unless disabled via `enable_ionospheric_correction`, the ionospheric delay is
corrected automatically using only the coefficients decoded from the navigation
messages. A single model is chosen for the whole solve and applied to every
satellite, preferring the more accurate global TEC models: NTCM-G if Galileo
Effective Ionisation Level coefficients have been decoded, else BDGIM if a BDS-3
B-CNAV (B1C/B2a/B2b) coefficient set has, else the GPS Klobuchar model, else
BeiDou's own Klobuchar variant from its legacy B1I/B3I message, else no
correction. See
[`select_ionospheric_correction`](@ref) and [`ionospheric_delay`](@ref).

Unless disabled via `enable_tropospheric_correction`, the tropospheric delay is
corrected with a blind model (no broadcast coefficients needed): Saastamoinen
zenith delays mapped to the line of sight by the Niell mapping functions.
See [`tropospheric_delay`](@ref).

# Arguments
- `groups`: this epoch's observations, as [`SignalGroups`](@ref) — a `NamedTuple` of
  [`SignalGroup`](@ref)s, one per ranging signal, each holding that signal's
  [`SignalObservation`](@ref)s. A bare `Tuple` of groups is accepted (and numbered
  `group1`, `group2`, …), as is a single `SignalGroup` (which becomes the `:default`
  group), so a one-constellation solve stays a one-liner:

  ```julia
  using PositionVelocityTime: SignalGroup
  calc_pvt(SignalGroup(GPSL1CA(), gps_states))
  calc_pvt((gps = SignalGroup(GPSL1CA(), gps_states),
            galileo = SignalGroup(GalileoE1B(), galileo_states)))
  ```

  Grouping is what makes the solve type-stable: within a group every satellite shares
  one concrete state type, so [`collect_measurements`](@ref) dispatches statically, and
  the solver behind it compiles once for every constellation mix. Build the groups where
  the satellites are tracked — with Tracking.jl 9.1 or later loaded,
  [`signal_groups`](@ref)`(track_state, decoders)` builds them from a whole `TrackState`.
  A pooled `Vector{SignalObservation}`, which is what this function took before 6.0, is
  refused with an error saying what to build instead.

  Each `(signal, PRN)` pair must appear at most once — a receiver produces one
  observation per signal per satellite, and a duplicate would enter the least-squares
  solve twice. A `Dictionary`-backed group prevents this within itself; two groups that
  share a ranging signal can still collide, and are not checked for.

  **Order is significant**, as vector order was before it: the flat observation order
  is group order × within-group order, and it fixes the reference system's tie-break
  ("most satellites, ties by first appearance"), [`ifb_layout`](@ref)'s
  reference-band tie-break, and the insertion order of `PVTSolution.observations`.
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin). A
  previous solution whose own DOP is implausible (`NaN`, or `PDOP` above
  `MAX_TRUSTED_WARM_START_PDOP` = 50) is not used as a seed and the epoch is
  solved from cold, so one spurious fix cannot re-seed itself epoch after epoch.

# Keyword Arguments
- `approximate_year`: Calendar year of the observation, used to resolve the
  GPS L1 C/A 1024-week rollover ambiguity (legacy LNAV broadcasts only a
  10-bit week number, so the receiver needs external information to
  determine which 1024-week cycle the recording is in). Anything within
  ±9 years of the actual observation date works. Defaults to the current
  UTC year, which is correct for live signals; for processing archived
  recordings, pass the rough year of the recording.
- `enable_ionospheric_correction`: when `true` (default), apply the automatic
  ionospheric correction described above. Set to `false` to skip it entirely
  and solve from the raw pseudoranges (e.g. for diagnostics or when an external
  correction is applied elsewhere).
- `enable_tropospheric_correction`: when `true` (default), apply the Saastamoinen
  tropospheric correction. Set to `false` to skip it.

# Returns
A new [`PVTSolution`](@ref) containing position, velocity, time, DOP values, and
satellite information; `prev_pvt` is never modified. [`calc_pvt!`](@ref) is the same
solve reusing the containers of a solution passed to it, without allocating. Returns `prev_pvt` if the epoch cannot be solved: too few healthy
satellites to solve the constellation (including the GGTO fallback and the
distinct-satellite condition — a satellite tracked on several bands supplies one line of
sight, so observations alone are not enough), a geometry whose solved design matrix is
rank deficient (reported as a `NaN` DOP). None of these throw, so a receiver can pass
whatever it currently tracks each epoch and carry `prev_pvt` forward. Distrust of
`prev_pvt` (see above) affects only the seed: an unsolvable epoch still returns
`prev_pvt` exactly as passed.
"""
function calc_pvt(
    groups,
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    # A fresh solution lends the solve its (empty) containers; on an unsolvable epoch the
    # solver hands `prev_pvt` itself back.
    _, solution = _calc_pvt!(PVTSolution(), PVTWorkspace(), groups, prev_pvt,
        approximate_year, enable_ionospheric_correction, enable_tropospheric_correction)
    solution
end

"""
    calc_pvt!(solution::PVTSolution, workspace::PVTWorkspace, groups, prev_pvt::PVTSolution;
              approximate_year::Integer = year(now(UTC)),
              enable_ionospheric_correction::Bool = true,
              enable_tropospheric_correction::Bool = true) -> PVTSolution

[`calc_pvt`](@ref), reusing the containers of `solution` and the buffers of `workspace`,
so that a steady stream of epochs is solved without allocating. The solve is the same —
same arguments, same keywords, same numbers — and so is the result it returns.

[`PVTSolution`](@ref) is immutable, so the fix is a new solution, returned; keep the
return value:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for groups in epochs
    pvt = calc_pvt!(pvt, workspace, groups, pvt)
end
```

**`solution` is the output argument, and the only thing overwritten**: the returned
solution takes over its `satellites`, `inter_system_biases` and `inter_frequency_biases`
containers, which are emptied and refilled with this epoch's values. `solution` itself
is therefore spent — its own fields still describe the epoch it was computed for, but
its containers now hold the new one — so use the returned solution from here on, and
keep a [`copy`](@ref Base.copy(::PVTSolution)) of `solution` if you need it after the
call (a `Dictionary` made by `map` over `solution.observations` shares its keys, and sees the
change too). `prev_pvt` is only read. Passing the same solution as both is allowed and
is the usual receiver loop above: every value the seed needs is read before anything is
written.

Where [`calc_pvt`](@ref) returns `prev_pvt` — an epoch that cannot be solved — so does
this: `prev_pvt` itself when it is `solution` (or shares its containers), and otherwise a
copy of it in `solution`'s containers, so that `prev_pvt` stays untouched.

It allocates nothing once `workspace` and `solution` have held an epoch with at least as
many satellites and as many estimated biases; a larger epoch grows them, once. That
covers everything from the collection pass to the returned solution, but not the
`groups` themselves, which the caller builds (and can reuse too — a `SignalGroup`'s
satellites may be any vector).

The `workspace` holds no state between epochs that affects the result — only buffers —
but it must not be used by two solves at once; give each task its own.
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

# The collection pass and the solve: `(true, fix)` with the fix in `solution`'s
# containers, or `(false, prev_pvt)` for an epoch that cannot be solved, with `solution`
# untouched. Specialises on the group shape, like `collect_measurements`, and is small;
# the solver behind it compiles once.
function _calc_pvt!(
    solution,
    workspace,
    groups,
    prev_pvt,
    approximate_year,
    enable_ionospheric_correction,
    enable_tropospheric_correction,
)
    # The function barrier. `_collect_measurements!` specialises on the group shape and
    # reduces it to observations; `_solve_pvt!` below is one compiled body for every
    # constellation mix there is.
    (; ionospheric_model, inter_system_time_offsets) =
        _collect_measurements!(workspace.prepared, groups, approximate_year)
    _solve_pvt!(
        solution,
        workspace,
        enable_ionospheric_correction ? ionospheric_model : IonosphericModel(nothing),
        inter_system_time_offsets,
        prev_pvt,
        enable_tropospheric_correction,
    )
end

# The solver: ~200 lines that must compile exactly once, so every argument is of a
# concrete type — the ionospheric model included, which is wrapped in an
# `IonosphericModel` rather than passed as the `Union` it is (see there). The inner
# barrier that does specialise on the model, once per model, where the work actually is,
# sits behind `predict_atmospheric_delays!`.
#
# Returns `(true, fix)`, the fix a new solution holding `solution`'s containers, or
# `(false, prev_pvt)` for an epoch that cannot be solved. The containers are written only
# once the epoch is known to be solvable, and only after everything is read from
# `prev_pvt`, which may be the same solution. A tuple rather than a
# `Union{Nothing,PVTSolution}`: the solution holds references, so such a `Union` would
# be returned boxed — an allocation per epoch.
function _solve_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    ionospheric_model,
    inter_system_time_offsets::BroadcastInterSystemTimeOffsets,
    prev_pvt::PVTSolution,
    enable_tropospheric_correction::Bool,
)
    prepared = workspace.prepared
    num_observations = length(prepared)

    # Solvability by the counts; a degenerate geometry is caught after the solve by the DOP.
    layout = workspace.layout
    decide_bias_layout!(layout, prepared, inter_system_time_offsets) || return (false, prev_pvt)
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
            get_time_system_offset(observation.time_system, reference_system)
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

    isnan(dop.GDOP) && return (false, prev_pvt)

    # The velocity solve starts from the previous solution where the position solve trusts
    # it, from rest otherwise (see `estimate_velocity` for why one step suffices).
    prev_vel_and_drift = SVector(prev_pvt.velocity[1], prev_pvt.velocity[2],
        prev_pvt.velocity[3], ustrip(m / s, prev_pvt.clock_drift))
    velocity_start =
        !distrusted && all(isfinite, prev_vel_and_drift) ?
        prev_vel_and_drift : zero(SVector{4,Float64})
    user_velocity_and_clock_drift, rate_residuals =
        estimate_velocity!(workspace.rate_residuals,
            workspace.velocity_design, prepared, position, velocity_start)
    all(isfinite, user_velocity_and_clock_drift) || return (false, prev_pvt)
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

    return true,
    PVTSolution(
        position,
        velocity,
        calc_course_over_ground(position, velocity),
        week_start + corrected_reference_time,
        user_velocity_and_clock_drift[4] * (m / s),
        reference_system,
        inter_system_biases,
        inter_frequency_biases,
        dop,
        observations,
    )
end
