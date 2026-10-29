@testset "calc_dop" begin
    # ECEF line-of-sight unit vector for a satellite at azimuth `az` (rad, from North) and
    # elevation `el` (rad), seen from `user`. Built by rotating the local ENU direction back
    # to ECEF via Geodesy, so the geometry matches what the solver feeds calc_dop (rows of
    # negated ECEF lines of sight, as calc_design_matrix! builds).
    enu_los(az, el) = [cos(el) * sin(az), cos(el) * cos(az), sin(el)]
    function ecef_los(user, az, el)
        e = enu_los(az, el)
        Vector(ECEFfromENU(user, wgs84)(ENU(e[1], e[2], e[3]))) .- Vector(user)
    end

    # Non-degenerate location (ECEF ≠ ENU axes, so an un-rotated ECEF DOP would
    # give visibly wrong HDOP/VDOP).
    user = ECEFfromLLA(wgs84)(LLA(50.1, 8.7, 120.0))
    azels = [
        (0.0, deg2rad(80)),
        (deg2rad(90), deg2rad(30)),
        (deg2rad(180), deg2rad(45)),
        (deg2rad(270), deg2rad(20)),
        (deg2rad(45), deg2rad(60)),
    ]

    # Single-system geometry: 3 position columns + 1 clock column.
    H_ecef = reduce(vcat, [[ecef_los(user, az, el)' 1.0] for (az, el) in azels])
    dop = PositionVelocityTime.calc_dop(H_ecef, user)

    # Reference DOP computed directly in the local ENU frame (RTKLIB-style: rows
    # are the ENU line-of-sight components). calc_dop must reproduce these from
    # the ECEF H — this is what pins the ECEF→ENU rotation (issue: HDOP/VDOP were
    # previously taken in ECEF and would fail here).
    H_enu = reduce(vcat, [[enu_los(az, el)' 1.0] for (az, el) in azels])
    Q = inv(H_enu' * H_enu)
    @test dop.HDOP ≈ sqrt(Q[1, 1] + Q[2, 2])
    @test dop.VDOP ≈ sqrt(Q[3, 3])
    @test dop.PDOP ≈ sqrt(Q[1, 1] + Q[2, 2] + Q[3, 3])
    @test dop.TDOP ≈ sqrt(Q[4, 4])
    @test dop.GDOP ≈ sqrt(tr(Q))

    # Internal consistency / ordering.
    @test dop.PDOP ≈ hypot(dop.HDOP, dop.VDOP)
    @test dop.HDOP ≤ dop.PDOP ≤ dop.GDOP
    @test dop.VDOP ≤ dop.PDOP
    @test all(>(0), (dop.GDOP, dop.PDOP, dop.HDOP, dop.VDOP, dop.TDOP))

    # Rank-deficient geometry (identical rows) → singular HᵀH → sentinel
    # `NaN` DOPs returned gracefully, not an exception.
    H_singular = repeat([1.0 0.0 0.0 1.0], 4)
    dop_singular = PositionVelocityTime.calc_dop(H_singular, user)
    @test all(isnan, (dop_singular.GDOP, dop_singular.PDOP, dop_singular.HDOP,
        dop_singular.VDOP, dop_singular.TDOP))

    # Multi-GNSS: 3 position + 2 clock columns (3 sats per system). GDOP/HDOP are
    # independent of which clock is the reference; TDOP reports the selected system's
    # clock variance, so it changes with reference_clock_column.
    azels6 = [azels; (deg2rad(300), deg2rad(35))]
    sys = [1, 1, 1, 2, 2, 2]
    H_multi = reduce(
        vcat,
        [
            [ecef_los(user, az, el)' (sys[i]==1 ? 1.0 : 0.0) (sys[i]==2 ? 1.0 : 0.0)] for
            (i, (az, el)) in enumerate(azels6)
        ],
    )
    d1 = PositionVelocityTime.calc_dop(H_multi, user, 1)
    d2 = PositionVelocityTime.calc_dop(H_multi, user, 2)
    @test d1.GDOP ≈ d2.GDOP
    @test d1.HDOP ≈ d2.HDOP
    @test d1.PDOP ≈ d2.PDOP
    @test !(d1.TDOP ≈ d2.TDOP)
    @test d1.GDOP > 0
end

@testset "calc_fit_covariance!" begin
    H = [1.0 0.0 1.0; 0.0 2.0 1.0; 1.0 1.0 1.0; -1.0 0.5 1.0]
    normal_matrix = fill(NaN, 3, 3)
    @test PositionVelocityTime.calc_fit_covariance!(normal_matrix, H)
    # The upper triangle holds the inverse; the lower is left as the factorisation had it.
    @test Symmetric(normal_matrix, :U) ≈ inv(H' * H)
    @test !PositionVelocityTime.calc_fit_covariance!(normal_matrix,
        repeat([1.0 0.0 1.0], 4))
end

# Regression: a degenerate geometry must be reported, not thrown out of a least-squares
# solve as a SingularException. Observed while processing a dual-band recording, where an
# epoch's geometry left the velocity normal equations singular and the 4×4 solve threw
# before the epoch could be rejected.
@testset "degenerate geometry is rejected without throwing" begin
    # The rank test `calc_dop` decides the geometry by: the normal-equations matrix `HᵀH`
    # of the design, Cholesky-factorised with a relative tolerance on its pivots.
    full_rank(H) = PositionVelocityTime.factorize_normal_matrix!(
        Matrix{Float64}(undef, size(H, 2), size(H, 2)), H)

    @testset "the normal-matrix factorisation rejects what is not positive definite" begin
        @test full_rank([sqrt(2.0) 0.0; 0.0 sqrt(3.0)])          # HᵀH positive definite
        @test !full_rank([1.0 1.0])                              # semidefinite: zero pivot
        # Ill-conditioned to the point of being rank deficient in Float64 (cond(HᵀH) =
        # 1e14): the pivot is positive, so a successful factorisation alone would accept
        # it and return rounding noise.
        @test !full_rank([1.0 0.0; 0.0 1e-7])
        # A poor but genuinely solvable geometry is still accepted (cond(HᵀH) = 1e6).
        @test full_rank([1.0 0.0; 0.0 1e-3])
    end

    @testset "the geometry checks follow the design matrix rank" begin
        # `calc_dop` and the velocity solve both decide rank through the normal-equations
        # matrix, so that is what these cases exercise.
        @test full_rank([1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 1.0 0.0])
        @test !full_rank(repeat([1.0 0.0 1.0], 3))               # identical rows
        # A single-band, single-system position design over three satellites: four
        # columns constrained by three lines of sight is always rank deficient — the shape
        # the layout's distinct-satellite condition rejects before the solve.
        @test !full_rank([1.0 0.0 0.0 1.0; 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0])

        # The case no satellite count can see: four *distinct* satellites — plenty for the
        # 3 + 1 unknowns — whose lines of sight lie on a common cone (equal projection onto
        # z), which puts the clock column in their span.
        cosθ = 0.5
        sinθ = sqrt(1 - cosθ^2)
        cone = [sinθ 0.0 cosθ; 0.0 sinθ cosθ; -sinθ 0.0 cosθ; 0.0 -sinθ cosθ]
        @test !full_rank([cone ones(4)])
        # Tilting one satellite off the cone restores full rank.
        tilted = [cone[1:3, :]; normalize([0.0, -sinθ, 2cosθ])']
        @test full_rank([tilted ones(4)])
    end

    # The velocity solve assumes a full-rank position design (its caller establishes that
    # with `calc_dop` first), so what is pinned here is the implication it relies on: a
    # degenerate line-of-sight set is degenerate for the velocity design too, and would
    # reach a singular solve if the DOP check were ever moved after it.
    @testset "a degenerate line-of-sight set is degenerate for the velocity design" begin
        states = gps_l1_states(0.0Hz)[1:4]
        measurements, _ = PositionVelocityTime.collect_measurements(
            signal_group(states); approximate_year = 2021)
        # The solve takes its lines of sight from the satellite positions, so the rows'
        # satellites are moved onto synthetic directions from a nearby user position.
        user = SVector(3.9074e6, 3.0684e5, 5.0150e6)
        with_position(m, position) = PositionVelocityTime.SatelliteMeasurement(
            (name === :position ? position : getfield(m, name)
             for name in fieldnames(PositionVelocityTime.SatelliteMeasurement))...)
        placed(dirs) = [with_position(m, user + 2.0e7 * normalize(SVector{3}(d)))
                        for (m, d) in zip(measurements, dirs)]
        # The position design's columns — the negated lines of sight, as
        # `calc_design_matrix!` builds them — and one column of ones for the per-system
        # clock columns, which the velocity solve collapses into one drift.
        design(dirs) = reduce(vcat, [[-normalize(d)' 1.0] for d in dirs])

        # Four well-spread directions ⇒ solvable: a finite [vx, vy, vz, ċ], plus one
        # finite post-fit range-rate residual per satellite.
        spread = [[1.0, 0.2, 0.9], [-0.5, 1.0, 0.7], [0.3, -1.0, 0.5], [0.0, 0.1, 1.0]]
        measurements = placed(spread)
        solution, rate_residuals =
            PositionVelocityTime.estimate_velocity(measurements, user)
        @test length(solution) == 4
        @test all(isfinite, solution)
        @test length(rate_residuals) == length(measurements)
        @test all(isfinite, rate_residuals)
        # A post-fit least-squares residual is orthogonal to every column of its design
        # — the gradient of the modelled rate at the solved state. That is the position
        # design's columns and a clock-drift column of ones but for the Earth's rotation
        # of the line of sight and the light-time factor's `s/c`, ~1e-5 here, where the
        # satellites' real velocities meet synthetic directions. That is what makes it a
        # residual rather than an arbitrary measured-minus-modeled difference, so it is
        # pinned directly instead of only by magnitude.
        let user_vel = SVector(solution[1], solution[2], solution[3])
            A = reduce(vcat, [last(PositionVelocityTime.pseudorange_rate_and_gradient(user,
                user_vel, solution[4], m.position, m.velocity, m.satellite_clock_drift))'
                              for m in measurements])
            @test norm(A - design(spread)) < 1e-4
            @test norm(A' * rate_residuals) < 1e-6
        end

        # The degenerate cases are stated on the design matrix, which `calc_dop` tests
        # before the velocity solve. Two distinct directions (a satellite tracked on a
        # second band repeats its line of sight) leave the position columns rank deficient.
        solvable(dirs) = full_rank(design(dirs))
        repeated = [[1.0, 0.2, 0.9], [1.0, 0.2, 0.9], [0.3, -1.0, 0.5], [0.3, -1.0, 0.5]]
        @test !solvable(repeated)

        # Three independent directions on a common cone (equal projection onto z) put the
        # drift column in their span, leaving the fourth pivot exactly zero — the shape
        # that produced SingularException(4) from the 4×4 solve.
        cosθ = 0.5
        sinθ = sqrt(1 - cosθ^2)
        cone = [
            [sinθ, 0.0, cosθ], [0.0, sinθ, cosθ], [-sinθ, 0.0, cosθ], [0.0, -sinθ, cosθ]]
        @test !solvable(cone)

        # Should such a geometry reach the velocity solve all the same, it tests its own
        # design — which differs from the position design by a few 1e-6 relative, too
        # little to restore its rank — and reports `NaN` instead of throwing.
        for dirs in (repeated, cone)
            solution, rate_residuals =
                PositionVelocityTime.estimate_velocity(placed(dirs), user)
            @test all(isnan, solution)
            @test all(isnan, rate_residuals)
        end
    end
end

@testset "residuals are observed minus computed" begin
    # The reported orientation is a convention rather than a by-product: `LsqFit`
    # residuates `model - data`, and the velocity solve's own measurement `yⱼ` carries
    # `-doppler * λ`, so either domain would report the opposite sign if left to itself.
    # Both are pinned here in the sense RTKLIB's `rescode` / `resdop` define — a
    # measurement exceeding the model reads positive — because a silent flip is
    # invisible to every magnitude-based assertion elsewhere in this suite.

    # ── Pseudorange ──────────────────────────────────────────────────────────────
    user = [3.9074e6, 3.0684e5, 5.0150e6]
    dirs = [[0.3, 0.1, 0.95], [-0.6, 0.3, 0.75], [0.2, -0.8, 0.55],
        [0.7, 0.5, 0.5], [-0.3, -0.5, 0.8], [0.0, 0.9, 0.45]]
    sat_positions = reduce(hcat, [user .+ 2.0e7 .* normalize(d) for d in dirs])
    bias_columns = PositionVelocityTime.BiasColumns(ones(Int, 6), 1, zeros(Int, 6), 0)
    ρ = PositionVelocityTime.predict_pseudoranges!(
        zeros(6), sat_positions, [user; 30.0], bias_columns)
    _, resid = PositionVelocityTime.estimate_position(sat_positions, ρ, bias_columns)
    @test maximum(abs, resid) < 1e-6   # self-consistent data ⇒ no residual

    # ±100 m on satellite 3's measurement alone. Least squares spreads part of it over
    # the four unknowns, so the satellite keeps `1 - P₃₃` of it — a fraction, but with
    # the sign of the perturbation: measuring long reads positive, short negative. (How
    # large the fraction is depends on the satellite's leverage, so only the sign and
    # the bound are pinned.)
    perturbed(δ) = last(PositionVelocityTime.estimate_position(
        sat_positions, (ρ_δ = copy(ρ); ρ_δ[3] += δ; ρ_δ), bias_columns))[3]
    @test 0.0 < perturbed(100.0) < 100.0
    @test -100.0 < perturbed(-100.0) < 0.0
    # Antisymmetric to sub-mm: the model is nonlinear in position, so the two do not
    # mirror each other exactly.
    @test perturbed(100.0) ≈ -perturbed(-100.0) atol = 1e-3

    # ── Range rate ───────────────────────────────────────────────────────────────
    # Reported in `yⱼ`'s geometric sense (positive while the satellite recedes), so a
    # *higher* measured Doppler — closing faster — reads NEGATIVE. This is the pairing
    # a receiver forming the same residual from its loops' `λ · carrier_doppler` has to
    # respect; getting it wrong inverts the field without changing any magnitude.
    # The fixture's own satellites, seen from `user` (near the fixture's fix).
    states = gps_l1_states(0.0Hz)[1:5]
    with_doppler(state, doppler) = PositionVelocityTime.SatelliteState(;
        decoder = state.decoder,
        signal = state.signal,
        code_phase = state.code_phase,
        carrier_phase = state.carrier_phase,
        carrier_doppler = doppler,
    )

    measurements(v) = measurements_of(signal_group(v))
    _, base = PositionVelocityTime.estimate_velocity(measurements(states), user)
    Δdoppler = 50.0Hz
    faster = copy(states)
    faster[3] = with_doppler(states[3], states[3].carrier_doppler + Δdoppler)
    _, closing = PositionVelocityTime.estimate_velocity(measurements(faster), user)
    λ = PositionVelocityTime.SPEED_OF_LIGHT / ustrip(Hz, get_center_frequency(states[3].signal))
    @test -λ * ustrip(Hz, Δdoppler) < closing[3] - base[3] < 0.0
end

# The velocity solve fits `calc_pseudorange_rate` itself: at the solved state, the
# measured rate minus `calc_pseudorange_rate` is the reported rate residual, to the
# precision of its one step (~1e-12 m/s from rest for this resting receiver). Real
# observables, so the residuals are not zero.
@testset "the velocity solve fits calc_pseudorange_rate's model" begin
    c = PositionVelocityTime.SPEED_OF_LIGHT
    groups = signal_group(gps_l1_states(0.0Hz))
    pvt = calc_pvt(groups; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    user_pos = SVector{3,Float64}(pvt.position)
    user_vel = SVector{3,Float64}(pvt.velocity)
    for measurement in measurements_of(groups)
        modelled = PositionVelocityTime.calc_pseudorange_rate(user_pos, user_vel,
            ustrip(m / s, pvt.clock_drift), measurement.position, measurement.velocity,
            measurement.satellite_clock_drift)
        measured = -measurement.carrier_doppler * c / measurement.center_frequency
        key = (measurement.signal_id, measurement.prn)
        residual = ustrip(m / s, pvt.satellites[key].rate_residual)
        @test abs(residual) > 1e-3
        @test measured - modelled ≈ residual atol = 1e-9
    end
end

# `measurement` with its carrier Doppler set so that its measured range rate is `rate`.
with_pseudorange_rate(measurement, rate) = PositionVelocityTime.SatelliteMeasurement(
    (name === :carrier_doppler ? -rate * measurement.center_frequency /
                                 PositionVelocityTime.SPEED_OF_LIGHT :
     getfield(measurement, name)
     for name in fieldnames(PositionVelocityTime.SatelliteMeasurement))...)

# The solve takes one Gauss–Newton step and relies on the model being nearly linear in
# `[v; ċ]`. Measured rates generated by `calc_pseudorange_rate` itself at a known state,
# for a receiver at an aircraft's speed with a large clock drift, pin how far one step
# from rest and from a nearby start lands from that state.
@testset "one Gauss–Newton step solves the velocity: $name" for (name, states) in (
    ("GPS L1 C/A", gps_l1_states(0.0Hz)),
    ("Galileo E1B", galileo_e1b_states(0.0Hz)),
)
    groups = signal_group(states)
    pvt = calc_pvt(groups; approximate_year = 2021)
    user_pos = SVector{3,Float64}(pvt.position)
    measurements = measurements_of(groups)
    for (user_vel, clock_drift) in ((SVector(150.0, -200.0, 50.0), 3000.0),
        (SVector(-600.0, 400.0, -300.0), -30000.0))
        truth = SVector(user_vel..., clock_drift)
        exact = [with_pseudorange_rate(measurement,
            PositionVelocityTime.calc_pseudorange_rate(user_pos, user_vel, clock_drift,
                measurement.position, measurement.velocity,
                measurement.satellite_clock_drift)) for measurement in measurements]
        # From rest: hundreds of m/s and kilometres per second of drift in one step.
        solved, residuals =
            PositionVelocityTime.estimate_velocity(exact, user_pos)
        @test norm(solved - truth) < 1e-6
        @test maximum(abs, residuals) < 1e-6
        # From a start a few m/s off, the previous epoch's solution, to rounding.
        start = truth + SVector(2.0, -1.5, 0.5, 0.3)
        solved, residuals =
            PositionVelocityTime.estimate_velocity(exact, user_pos, start)
        @test norm(solved - truth) < 1e-10
        @test maximum(abs, residuals) < 1e-10
    end
end
