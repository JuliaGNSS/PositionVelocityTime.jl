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

The reusable scratch storage of [`calc_pvt!`](@ref): the flat measurement rows, the bias
layout, the pseudoranges and their corrections, the least-squares buffers, the design
and normal-equations matrices, the velocity solve's design rows and the residuals. Every buffer grows to the largest epoch
it has seen and is reused from then on; its contents after a solve are not part of any
interface.

Satellite positions are kept as a `Vector{SVector{3,Float64}}` — a fixed-size quantity
per satellite, in a vector that is resized in place. The least-squares state, whose
length `3 + M + B` depends on the epoch's time systems and bands, is a plain vector
instead of a static one, so the solver compiles once rather than once per layout; the
matrices sized by it are kept at their largest and used through views.
"""
mutable struct PVTWorkspace
    measurements::Vector{SatelliteMeasurement}
    layout::BiasLayoutWorkspace
    unique_systems::Vector{SupportedTimeSystem}
    sat_positions::Vector{SVector{3,Float64}}
    pseudoranges::Vector{Float64}
    hub_offsets::Vector{Float64}
    atmospheric_delays::Vector{Float64}
    corrected_ranges::Vector{Float64}
    prev_ξ::Vector{Float64}
    residuals::Vector{Float64}
    rate_residuals::Vector{Float64}
    velocity_design::Vector{SVector{4,Float64}}
    design_matrix::Matrix{Float64}
    normal_matrix::Matrix{Float64}
    least_squares::LMWorkspace
end

PVTWorkspace() = PVTWorkspace(
    SatelliteMeasurement[],
    BiasLayoutWorkspace(),
    SupportedTimeSystem[],
    SVector{3,Float64}[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    SVector{4,Float64}[],
    Matrix{Float64}(undef, 0, 0),
    Matrix{Float64}(undef, 0, 0),
    LMWorkspace(),
)

"""
    calc_pvt(groups, prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS satellite measurements.

Satellites from different GNSS and frequency bands may be combined. Because each
constellation references its broadcasts to its own system time, one receiver clock
bias is estimated per GNSS time system; and because the receiver's RF chain delays
each band differently, one receiver inter-frequency bias is estimated per frequency
band beyond a reference band (shared across constellations on that band; see
`get_band_id`). The state vector is therefore
`[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]` for `M` distinct time systems and `B`
extra bands. Position and time are found by least squares; velocity and clock drift
are solved from carrier Doppler.

A solution requires `n ≥ 3 + M + B` healthy satellite measurements (each system needs at
least one satellite, and a system contributing a single satellite spends it entirely
on that system's clock bias), of which `3 + M` must come from *distinct* satellites: a
satellite tracked on several bands supplies one line of sight, and its extra
measurements constrain the inter-frequency biases rather than the geometry. When either
condition fails, constellations whose messages carry a broadcast offset to another
tracked system — Galileo's GGTO (Galileo–GPS Time Offset), BeiDou's BGTO (BDT–GNSS
Time Offset, toward GPS or Galileo), or the GGTO GPS itself broadcasts on CNAV/CNAV-2
— have their clock bias collapsed onto that hub system using the broadcast offset,
which makes a 4-satellite mixed fix possible (see [`decide_bias_layout`](@ref) for the
hub order). Estimating an independent bias is preferred whenever the geometry allows
it, since it avoids the broadcast offset's own error.

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
- `groups`: this epoch's measurements, as [`SignalGroups`](@ref) — a `NamedTuple` of
  [`SignalGroup`](@ref)s, one per ranging signal, each holding that signal's
  [`SatelliteState`](@ref)s. A bare `Tuple` of groups is accepted (and numbered
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
  A pooled `Vector{SatelliteState}`, which is what this function took before 6.0, is
  refused with an error saying what to build instead.

  Each `(signal, PRN)` pair must appear at most once — a receiver produces one
  measurement per signal per satellite, and a duplicate would enter the least-squares
  solve twice. A `Dictionary`-backed group prevents this within itself; two groups that
  share a ranging signal can still collide, and are not checked for.

  **Order is significant**, as vector order was before it: the flat measurement order
  is group order × within-group order, and it fixes the primary system's tie-break
  ("most satellites, ties by first appearance"), [`band_ifb_layout`](@ref)'s
  reference-band tie-break, and the insertion order of `PVTSolution.sats`.
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin). A
  previous solution whose own DOP is implausible (`GDOP < 0`, or `PDOP` above
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
sight, so measurements alone are not enough), a geometry whose solved design matrix is
rank deficient (reported as a negative GDOP). None of these throw, so a receiver can pass
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
solution takes over its `sats`, `inter_system_biases` and `inter_frequency_biases`
containers, which are emptied and refilled with this epoch's values. `solution` itself
is therefore spent — its own fields still describe the epoch it was computed for, but
its containers now hold the new one — so use the returned solution from here on, and
keep a [`copy`](@ref Base.copy(::PVTSolution)) of `solution` if you need it after the
call (a `Dictionary` made by `map` over `solution.sats` shares its keys, and sees the
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
    # reduces it to flat rows; `_solve_pvt!` below is one compiled body for every
    # constellation mix there is.
    ionospheric_correction =
        _collect_measurements!(workspace.measurements, groups, approximate_year)
    _solve_pvt!(
        solution,
        workspace,
        IonosphericModel(enable_ionospheric_correction ? ionospheric_correction : nothing),
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
    prev_pvt::PVTSolution,
    enable_tropospheric_correction::Bool,
)
    measurements = workspace.measurements

    # Gauss-Newton converges within its seed's basin, so a spurious far-away root
    # would re-seed itself epoch after epoch. Such roots flag themselves with an
    # impossible geometry (negative GDOP, exploded PDOP): solve from a cold seed
    # instead. Only the seed is affected — an unsolvable epoch still returns
    # `prev_pvt` unchanged.
    distrusted =
        !isnothing(prev_pvt.dop) &&
        (prev_pvt.dop.GDOP < 0 || prev_pvt.dop.PDOP > MAX_TRUSTED_WARM_START_PDOP)

    num_sats = length(measurements)

    # Every satellite that was not decoded-complete-and-healthy has already been
    # dropped by the collection pass, and each surviving row carries the GNSSSignals
    # keys that drive the solution: `time_system` (`GPST()`/`GST()`/`BDT()`) groups the
    # receiver clock bias — one per time system, ordered by first appearance —
    # and `band_id` (`:L1`, `:L5`, …) groups the receiver inter-frequency bias — one
    # per band beyond a per-coverage-component reference. (`signal_id`, `:GPSL1CA` …, is
    # the per-signal `sats` identity used below, not a grouping key.)
    #
    # Solvability is decided here: `decide_bias_layout!` keeps only observable IFBs and
    # falls back to a hub collapse when the geometry is disconnected or
    # under-determined. A degenerate geometry, which no count can see, is caught after
    # the solve by the DOP.
    decide_bias_layout!(workspace.layout, measurements) || return (false, prev_pvt)
    (; bias_columns, clock_systems, extra_bands, reference_bands, hub_system, hub_rows) =
        bias_layout(workspace.layout)
    (; clock_bias_indices, num_clock_biases) = bias_columns
    num_params = num_lsq_params(bias_columns)

    # The propagated ephemerides, already evaluated per satellite by the collection pass.
    sat_positions = resize!(workspace.sat_positions, num_sats)
    for j in 1:num_sats
        sat_positions[j] = measurements[j].position
    end

    # Primary system — its clock bias, reference time, week and start epoch
    # define the reported time (see `select_primary_system`).
    unique_systems = unique_time_systems!(workspace.unique_systems,
        measurement.time_system for measurement in measurements)
    primary_system = select_primary_system(measurements, hub_system)
    # `something`, not a bare `findfirst`: the primary system is either one of the
    # systems present or the hub, which `decide_bias_layout!` only picks when it is
    # present — so this never fails, and saying so keeps the index an `Int` rather than
    # a `Union{Int,Nothing}` that would widen the clock column and the DOP call below.
    primary_index = something(
        findfirst(measurement -> measurement.time_system === primary_system, measurements))
    primary_clock_index = clock_bias_indices[primary_index]

    # The common reference cancels out of the reported time (the primary clock
    # bias absorbs it), so any latest-transmit-time reference works — but the times
    # must first be put on one count (`calc_time_scale_offsets!`).
    # A BeiDou second-of-week reads 14 s below a GPS time of week for the same instant,
    # so differencing them raw hands every BeiDou measurement 4.2e9 m of structural
    # offset: not merely a biased BeiDou clock column, but a parameter nine orders of
    # magnitude above the others in the normal equations. Zero for GPS and Galileo, and
    # for any single-constellation epoch.
    # The offsets are written straight into the pseudorange buffer, and the transmit
    # times added onto them, so no separate offset buffer is needed.
    pseudoranges =
        calc_time_scale_offsets!(workspace.pseudoranges, measurements, primary_system)
    for (j, measurement) in enumerate(measurements)
        pseudoranges[j] += measurement.time
    end
    _, reference_time = calc_pseudoranges!(pseudoranges)
    # The known per-satellite broadcast steering offset (zero unless that system was
    # collapsed onto the hub), as a measurement correction like the atmospheric delays
    # below. Kept for the inter-system-bias readout at the end, which reports this
    # term alone: the time-count difference already removed above is a convention, not
    # a bias, and belongs in neither the readout nor the solve.
    hub_offsets =
        calc_hub_range_offsets!(workspace.hub_offsets, measurements, hub_rows, hub_system)
    pseudoranges .-= hub_offsets

    # The primary system's week and start epoch date the epoch absolutely: the day of
    # year feeds the tropospheric mapping's seasonal term here, and week/start epoch
    # date the reported time after the solve.
    primary_measurement = measurements[primary_index]
    week = primary_measurement.week
    start_time = primary_measurement.system_start_epoch
    doy = day_of_year(primary_measurement.system_start_time, week, reference_time)

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias. A
    # distrusted previous solution seeds nothing: the state stays all zeros, the cold
    # start.
    #
    # Per column, from the column's own system, not per satellite: a collapsed system's
    # satellites share the hub's column, and their own previous bias (the hub's plus the
    # collapsed inter-system bias) would overwrite the hub's seed.
    #
    # A system the previous solution has no bias for (it was not tracked then) is seeded
    # with the previous reference clock. That is off by the inter-system bias, typically
    # metres; the alternative, zero, is off by the whole receiver clock bias, often
    # orders of magnitude more.
    prev_abs_bias(sys) =
        sys === prev_pvt.reference_system ? prev_pvt.time_correction :
        prev_pvt.time_correction + get(prev_pvt.inter_system_biases, sys, 0.0m)
    prev_ξ = fill!(resize!(workspace.prev_ξ, num_params), 0.0)
    if !distrusted
        prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
        for (k, sys) in enumerate(clock_systems)
            prev_ξ[3+k] = ustrip(m, prev_abs_bias(sys))
        end
        # Only warm-start an IFB column from the previous solution when that band's bias
        # was measured against the same reference band; a reference change (the anchor
        # differs across epochs) makes the stored value refer to a different quantity,
        # so seed at 0 instead. The IFB enters the design linearly, so a stale seed only
        # costs iterations, but this keeps the starting point meaningful.
        #
        # `haskey` and an index, not `get(…, nothing)`: an `InterFrequencyBias` holds a
        # `Symbol`, so a `Union` of it with `nothing` is returned boxed.
        prev_ifbs = prev_pvt.inter_frequency_biases
        for (i, band) in enumerate(extra_bands)
            prev_ξ[3+num_clock_biases+i] =
                haskey(prev_ifbs, band) && prev_ifbs[band].reference == reference_bands[i] ?
                ustrip(m, prev_ifbs[band].value) : 0.0
        end
    end

    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model was chosen for the whole solve by the
    # collection pass. The prediction is `predict_atmospheric_delays!`, a function
    # barrier that specialises on the model this body deliberately does not.
    #
    # Is any atmospheric correction active at all? When neither the ionosphere (no
    # model selected) nor the troposphere contributes, skip the prediction entirely and
    # the solve runs on the raw pseudoranges.
    correct_atmosphere =
        !isnothing(ionospheric_model.correction) || enable_tropospheric_correction
    function corrected_ranges(ξ_prediction)
        delays = predict_atmospheric_delays!(
            workspace.atmospheric_delays, ξ_prediction, measurements, ionospheric_model,
            reference_time, doy, enable_tropospheric_correction)
        resize!(workspace.corrected_ranges, num_sats) .= pseudoranges .- delays
    end
    least_squares = workspace.least_squares

    # `ξ` is the least-squares workspace's own parameter vector, so it stays valid until
    # the next solve on that workspace — which is after everything below.
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, residuals_uncorrected = user_position!(least_squares,
            workspace.residuals, sat_positions, pseudoranges, bias_columns, prev_ξ)
        if correct_atmosphere
            user_position!(least_squares, workspace.residuals, sat_positions,
                corrected_ranges(ξ_uncorrected), bias_columns, ξ_uncorrected)
        else
            (ξ_uncorrected, residuals_uncorrected)
        end
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction. With no
        # active correction, solve directly on the raw pseudoranges.
        user_position!(least_squares, workspace.residuals, sat_positions,
            correct_atmosphere ? corrected_ranges(prev_ξ) : pseudoranges, bias_columns,
            prev_ξ)
    end
    workspace.design_matrix = grown(workspace.design_matrix, num_sats, num_params)
    H = calc_H!(view(workspace.design_matrix, 1:num_sats, 1:num_params), sat_positions, ξ,
        bias_columns)
    position = ECEF(ξ[1], ξ[2], ξ[3])

    # Check the geometry at the converged position — the DOP is reported to the caller, and
    # a rank deficiency here (negative GDOP) means `ξ` is meaningless, so nothing should be
    # derived from it. The satellite conditions in `decide_bias_layout!` already reject the
    # count-shaped degeneracies before the solve; this catches what a count cannot see.
    #
    # It also stands ahead of the velocity solve, whose design is this one's but for a few
    # 1e-6 relative; that solve tests its own design against the same tolerance as a
    # backstop and returns `NaN` should it fail, an unsolvable epoch as well (see
    # `calc_user_velocity_and_clock_drift`).
    workspace.normal_matrix = grown(workspace.normal_matrix, num_params, num_params)
    dop = calc_DOP!(view(workspace.normal_matrix, 1:num_params, 1:num_params), H, position,
        primary_clock_index)

    dop.GDOP < 0 && return (false, prev_pvt)

    # The velocity solve takes one Gauss–Newton step, which is the least-squares fit to
    # rounding from the previous epoch's velocity and clock drift, and to ~1e-8 m/s from
    # rest (see `calc_user_velocity_and_clock_drift`). So it starts from the previous
    # solution whenever the position solve trusts that as its seed, and from rest
    # otherwise: on a cold start (the default `prev_pvt` is at rest anyway), after a
    # distrusted solution, or should the previous one carry a non-finite value.
    previous_velocity_and_drift = SVector(prev_pvt.velocity[1], prev_pvt.velocity[2],
        prev_pvt.velocity[3], prev_pvt.relative_clock_drift * SPEED_OF_LIGHT)
    velocity_start =
        !distrusted && all(isfinite, previous_velocity_and_drift) ?
        previous_velocity_and_drift : zero(SVector{4,Float64})
    user_velocity_and_clock_drift, rate_residuals =
        calc_user_velocity_and_clock_drift!(workspace.rate_residuals,
            workspace.velocity_design, measurements, position, velocity_start)
    all(isfinite, user_velocity_and_clock_drift) || return (false, prev_pvt)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    time_correction = ξ[3+primary_clock_index]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEED_OF_LIGHT

    # Everything is read from `prev_pvt` by now, so the containers of `solution` —
    # possibly `prev_pvt`'s own — are written from here on.
    #
    # Per-satellite `sats` key: (signal id, PRN) — signal-level (not time system), so a
    # satellite tracked on two signals of one constellation stays distinct; the
    # receiver-clock grouping is separate, by time system.
    #
    # A duplicate key is refused, as `Dictionary(keys, values)` refused it before. It is
    # checked before any container is written, so the throw leaves `solution` — possibly
    # `prev_pvt` itself — as it was, not half-refilled. Not with `insert!`, which formats
    # the key into its message — a dynamic `show` a `juliac --trim` build rejects — so
    # the message here is fixed.
    for j in 2:num_sats, k in 1:j-1
        measurements[j].signal_id === measurements[k].signal_id &&
            measurements[j].prn == measurements[k].prn && throw(IndexError(
                "a (signal, PRN) pair appears twice in one epoch; each must appear once"))
    end
    sats = empty_keeping_capacity!(solution.sats)
    for (j, measurement) in enumerate(measurements)
        key = (measurement.signal_id, measurement.prn)
        set!(sats, key, SatInfo(measurement.position, measurement.time, residuals[j] * m,
            rate_residuals[j] * (m / s)))
    end

    # Inter-system biases relative to the reference (primary) system's clock, in
    # meters. The reference is omitted (its bias is `time_correction`); for a
    # collapsed system this is the broadcast offset −c·Δt_systems, `hub_offsets` as
    # `inter_system_bias` computed it, read at the system's first satellite (the
    # per-satellite offsets differ only by the offset polynomial's drift term over the
    # transmit-time spread — sub-millimetre).
    inter_system_biases = empty!(solution.inter_system_biases)
    for sys in unique_systems
        sys === primary_system && continue
        j = something(findfirst(measurement -> measurement.time_system === sys, measurements))
        inter_system_biases[sys] =
            (ξ[3+clock_bias_indices[j]] + hub_offsets[j] - time_correction) * m
    end

    # Receiver inter-frequency biases relative to the reference band, in meters
    # (the reference band is omitted; its bias is folded into the clock biases).
    inter_frequency_biases = empty!(solution.inter_frequency_biases)
    for (i, band) in enumerate(extra_bands)
        inter_frequency_biases[band] =
            InterFrequencyBias(ξ[3+num_clock_biases+i] * m, reference_bands[i])
    end

    return true,
    PVTSolution(
        position,
        velocity,
        calc_course_over_ground(position, velocity),
        time_correction * m,
        # Assumes `start_time.fraction == 0` (true for GPS/Galileo: integer-second
        # origins).
        TAITime(
            week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) +
            start_time.second,
            corrected_reference_time - floor(Int, corrected_reference_time),
        ),
        user_velocity_and_clock_drift[4] / SPEED_OF_LIGHT,
        dop,
        sats,
        primary_system,
        inter_system_biases,
        inter_frequency_biases,
    )
end
