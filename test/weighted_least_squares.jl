# Weighting the least-squares solve by measurement uncertainty (C/N₀ + elevation)
# instead of treating every satellite as equally precise. Two things are pinned here:
# the variance model itself, and that the weights actually buy accuracy — a fix whose
# satellites carry the noise their C/N₀ implies is measurably better when weighted.

using Random: MersenneTwister

const C_LIGHT = PositionVelocityTime.SPEEDOFLIGHT

# A satellite's code-chip length in metres: one chip of code phase is one chip length
# of pseudorange (the pseudorange is `(t_ref − t_transmit)·c` and the code phase
# advances the transmit time, so *more* code phase is *less* range).
chip_length(system) = C_LIGHT / ustrip(Hz, get_code_frequency(system))

"Copy of `state` with the given fields replaced."
function with_state(
    state;
    code_phase = state.code_phase,
    carrier_doppler = state.carrier_doppler,
    cn0 = state.cn0,
    pseudorange_variance = state.pseudorange_variance,
)
    SatelliteState(;
        decoder = state.decoder,
        system = state.system,
        code_phase,
        carrier_doppler,
        carrier_phase = state.carrier_phase,
        cn0,
        pseudorange_variance,
    )
end

"Copy of `state` whose pseudorange carries an additional error of `error_metres`."
function with_range_error(state, error_metres)
    chips = error_metres / chip_length(state.system)
    with_state(state; code_phase = state.code_phase - chips)
end

"The same states with their C/N₀ removed, i.e. the epoch as unweighted least squares."
strip_cn0(states) = [with_state(s; cn0 = nothing) for s in states]

rms(x) = sqrt(sum(abs2, x) / length(x))

"Copy of `state` whose range rate carries an additional error of `error_m_per_s`."
function with_range_rate_error(state, error_m_per_s)
    λ = C_LIGHT / ustrip(Hz, get_center_frequency(state.system))
    with_state(state; carrier_doppler = state.carrier_doppler - error_m_per_s / λ * Hz)
end

@testset "pseudorange variance model" begin
    pseudorange_variance = PositionVelocityTime.pseudorange_variance
    gpsl1 = GPSL1CA()
    zenith = deg2rad(90)

    # Noise falls with C/N₀ — the whole point of weighting. The C/N₀ term is the
    # coherent-DLL thermal jitter `(c/f_chip)·√(B_n·d/2)·10^(−CN0/20)`, ≈ 1.2 m at
    # 45 dBHz for GPS L1 C/A, so a 30 dBHz satellite is several times noisier.
    σ(cn0, el = zenith, system = gpsl1) = sqrt(pseudorange_variance(system, cn0, el))
    @test σ(30.0dBHz) > σ(35.0dBHz) > σ(45.0dBHz) > σ(55.0dBHz)
    thermal(cn0) = sqrt(σ(cn0)^2 - σ(nothing)^2)
    @test thermal(45.0dBHz) ≈ chip_length(gpsl1) * sqrt(0.5) * 10^(-45 / 20) rtol = 1e-12
    # 15 dB of C/N₀ is a factor 10^(15/20) ≈ 5.6 in the thermal term.
    @test thermal(30.0dBHz) / thermal(45.0dBHz) ≈ 10^(15 / 20)

    # A wider-band signal has a shorter chip and thus a sharper correlation peak: the
    # same C/N₀ on GPS L5 (10.23 Mcps) is ten times less code noise than on L1 C/A.
    thermal_l5 = sqrt(pseudorange_variance(GPSL5I(), 45.0dBHz, zenith) -
                      pseudorange_variance(GPSL5I(), nothing, zenith))
    @test thermal(45.0dBHz) / thermal_l5 ≈ 10 rtol = 1e-6

    # Elevation mapping: noise grows towards the horizon, saturates below the
    # tropospheric model's low-elevation bound rather than diverging (including for a
    # satellite that the current position estimate puts below the horizon).
    @test σ(45.0dBHz, deg2rad(90)) < σ(45.0dBHz, deg2rad(30)) < σ(45.0dBHz, deg2rad(5))
    @test σ(45.0dBHz, deg2rad(1)) == σ(45.0dBHz, deg2rad(-30)) ==
          σ(45.0dBHz, PositionVelocityTime._LOW_ELEVATION_THRESHOLD)

    # Unknown inputs drop their term instead of guessing: no elevation (cold start, no
    # position yet) and no C/N₀ (a caller that does not report it) are both fine, and
    # leave the floor that keeps a strong satellite from taking over the fix.
    @test pseudorange_variance(gpsl1, nothing, nothing) ==
          PositionVelocityTime._RESIDUAL_UERE_SIGMA^2
    @test σ(nothing, deg2rad(30)) < σ(30.0dBHz, deg2rad(30))
    @test σ(45.0dBHz, nothing) < σ(45.0dBHz, deg2rad(30))

    # Robustness ceiling: a satellite whose C/N₀ reads absurdly low is heavily
    # de-weighted but never silently dropped, and the floor keeps every σ meaningful.
    # `0.0dBHz` — what Tracking's estimator reports before it has integrated a prompt —
    # counts as "not reported", not as an unusable signal.
    lo, hi = PositionVelocityTime._PSEUDORANGE_SIGMA_BOUNDS
    @test σ(5.0dBHz) == hi
    @test lo < σ(90.0dBHz, deg2rad(90)) < σ(-10.0dBHz, deg2rad(0)) ≤ hi
    @test σ(0.0dBHz) == σ(nothing)
    @test σ(-5.0dBHz) == σ(nothing)
    @test !PositionVelocityTime._is_usable_cn0(nothing)
    @test !PositionVelocityTime._is_usable_cn0(0.0dBHz)
    @test PositionVelocityTime._is_usable_cn0(20.0dBHz)

    # C/N₀ enters as a linear ratio: 45 dBHz is 10^4.5 Hz.
    @test PositionVelocityTime._linear_cn0(45.0dBHz) ≈ 10^4.5

    # The per-satellite form takes the state's C/N₀, and an explicit variance replaces
    # the model outright (elevation and C/N₀ then do not matter).
    state = with_state(gps_l1_states(0.0Hz)[1]; cn0 = 30.0dBHz)
    @test pseudorange_variance(state, zenith) ==
          pseudorange_variance(gpsl1, 30.0dBHz, zenith)
    override = with_state(state; pseudorange_variance = (7.0m)^2)
    @test pseudorange_variance(override, zenith) == 49.0
    @test pseudorange_variance(override, deg2rad(5)) == 49.0
    # A supplied variance is bounded like a modeled one, so no caller can hand a
    # satellite an infinite (or `NaN`) weight and let it dictate the fix.
    bounded(variance) = pseudorange_variance(
        with_state(state; pseudorange_variance = variance), zenith)
    @test bounded(0.0m^2) == lo^2
    @test bounded((1e6m)^2) == hi^2
    @test bounded(NaN * m^2) == hi^2
    @test PositionVelocityTime.has_measurement_uncertainty(override)
    @test PositionVelocityTime.has_measurement_uncertainty(state)
    unreported(cn0) = !PositionVelocityTime.has_measurement_uncertainty(
        with_state(state; cn0))
    @test unreported(nothing)
    @test unreported(0.0dBHz)
end

@testset "range rate variance model" begin
    range_rate_variance = PositionVelocityTime.range_rate_variance
    gpsl1 = GPSL1CA()
    σv(cn0, system = gpsl1) = sqrt(range_rate_variance(system, cn0))

    # Doppler noise follows the carrier loop, so it has its own C/N₀ dependence — and
    # its own floor, which is all that is left when no C/N₀ is reported.
    @test σv(30.0dBHz) > σv(40.0dBHz) > σv(50.0dBHz) > σv(nothing)
    @test σv(nothing) == PositionVelocityTime._RANGE_RATE_SIGMA_FLOOR
    @test σv(0.0dBHz) == σv(nothing)
    thermal(cn0, system = gpsl1) = sqrt(σv(cn0, system)^2 - σv(nothing)^2)
    @test thermal(30.0dBHz) / thermal(45.0dBHz) ≈ 10^(15 / 20)
    # A metre-per-second scale, not a metre one: the FLL sees the carrier, not the code.
    @test 0.05 < σv(30.0dBHz) < 0.2
    # The clamp is a guard, so every reading stays inside it.
    lo, hi = PositionVelocityTime._RANGE_RATE_SIGMA_BOUNDS
    @test all(cn0 -> lo ≤ σv(cn0) ≤ hi, (0.1dBHz, 5.0dBHz, 45.0dBHz, 90.0dBHz))
    # Longer wavelength, more metres per Hz of Doppler error: L5 (1176 MHz) over
    # L1 (1575 MHz).
    @test thermal(45.0dBHz, GPSL5I()) > thermal(45.0dBHz)
    @test range_rate_variance(with_state(gps_l1_states(0.0Hz)[1]; cn0 = 30.0dBHz)) ==
          range_rate_variance(gpsl1, 30.0dBHz)
end

@testset "an epoch without reported uncertainty is solved exactly as before" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    states = gps_l1_states(0.0Hz)
    pvt = calc_pvt(states; kwargs...)

    # No satellite reports a C/N₀ or a variance, so nothing is weighted: every
    # satellite is credited the same nominal UERE …
    nominal = PositionVelocityTime._NOMINAL_PSEUDORANGE_SIGMA * m
    @test all(info -> info.pseudorange_sigma == nominal, values(pvt.sats))
    # … and the formal accuracy is then exactly the familiar DOP × UERE product, which
    # is what makes it comparable to the geometric DOP it is reported alongside.
    @test pvt.accuracy.horizontal ≈ pvt.dop.HDOP * nominal
    @test pvt.accuracy.vertical ≈ pvt.dop.VDOP * nominal
    @test pvt.accuracy.position ≈ pvt.dop.PDOP * nominal
    @test pvt.accuracy.time ≈ pvt.dop.TDOP * nominal

    # Reporting a C/N₀ leaves the DOP untouched — it is a geometric quantity and must
    # not silently become a weighted one — while the accuracy does change. (The DOPs are
    # equal to a few parts in 1e8, not bit-identical: they are read off the design matrix
    # at each solve's own converged position, and those sit a few centimetres apart.)
    weighted = calc_pvt([with_state(s; cn0 = 45.0dBHz) for s in states]; kwargs...)
    @test weighted.dop.GDOP ≈ pvt.dop.GDOP rtol = 1e-6
    @test weighted.dop.HDOP ≈ pvt.dop.HDOP rtol = 1e-6
    @test weighted.dop.PDOP ≈ pvt.dop.PDOP rtol = 1e-6
    @test weighted.accuracy.position < pvt.accuracy.position
    @test all(info -> info.pseudorange_sigma < nominal, values(weighted.sats))
end

@testset "formal accuracy" begin
    calc_formal_accuracy = PositionVelocityTime.calc_formal_accuracy
    user = ECEF(ECEFfromLLA(wgs84)(LLA(50.1, 8.7, 120.0)))
    function ecef_los(az, el)
        enu = ENU(cos(el) * sin(az), cos(el) * cos(az), sin(el))
        Vector(ECEFfromENU(user, wgs84)(enu)) .- Vector(user)
    end
    azels = [(0.0, deg2rad(80)), (deg2rad(90), deg2rad(30)), (deg2rad(180), deg2rad(45)),
        (deg2rad(270), deg2rad(20)), (deg2rad(45), deg2rad(60))]
    H = reduce(vcat, [[ecef_los(az, el)' 1.0] for (az, el) in azels])

    # Uniform variances reproduce DOP × σ, the classic accuracy estimate.
    dop = PositionVelocityTime.calc_DOP(H, user)
    accuracy = calc_formal_accuracy(H, fill(4.0, 5), user)
    @test accuracy.horizontal ≈ dop.HDOP * 2.0m
    @test accuracy.vertical ≈ dop.VDOP * 2.0m
    @test accuracy.position ≈ dop.PDOP * 2.0m
    @test accuracy.time ≈ dop.TDOP * 2.0m
    # Horizontal/vertical split of the same 3D quantity, taken in the ENU tangent plane.
    @test accuracy.position ≈ hypot(accuracy.horizontal, accuracy.vertical)

    # It is a covariance, so it scales with the assumed measurement σ …
    @test calc_formal_accuracy(H, fill(16.0, 5), user).position ≈ 2 * accuracy.position
    # … and a satellite trusted less can only make the fix less certain.
    degraded = calc_formal_accuracy(H, [100.0, 4.0, 4.0, 4.0, 4.0], user)
    @test degraded.position > accuracy.position
    # A rank-deficient geometry is reported, not thrown — as in `calc_DOP`.
    @test isnothing(calc_formal_accuracy(repeat([1.0 0.0 0.0 1.0], 4), fill(4.0, 4), user))
end

@testset "a satellite trusted less pulls the fix less" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    states = gps_l1_states(0.0Hz)

    # Five satellites, one of which is handed a 60 m pseudorange error. Weighted with an
    # (over-)confident σ on the four good ones and the model ceiling on the broken one,
    # the fix should end up close to what the four good satellites alone give — the
    # broken satellite contributes its geometry, not its error.
    good = [with_state(s; pseudorange_variance = (0.5m)^2) for s in states[1:4]]
    broken = with_state(with_range_error(states[5], 60.0); pseudorange_variance = (50.0m)^2)
    without_broken = calc_pvt(states[1:4]; kwargs...)
    unweighted = calc_pvt([states[1:4]; [with_range_error(states[5], 60.0)]]; kwargs...)
    weighted = calc_pvt([good; [broken]]; kwargs...)
    @test norm(weighted.position - without_broken.position) <
          0.05 * norm(unweighted.position - without_broken.position)

    # The residual reported per satellite stays the raw metre residual, weighted or not,
    # so it remains comparable across epochs; the σ it was weighted by is reported
    # alongside it, which is what turns it into a normalised residual for fault
    # detection. The broken satellite is the one that stands out.
    broken_info = weighted.sats[(:GPSL1CA, broken.decoder.prn)]
    @test broken_info.pseudorange_sigma == 50.0m
    @test abs(broken_info.residual) > 40m
    @test all(
        info -> abs(info.residual) < 5m,
        [weighted.sats[(:GPSL1CA, s.decoder.prn)] for s in good],
    )
end

@testset "the weighted velocity solve keeps reporting raw range-rate residuals" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    # A weak satellite whose Doppler is off by 1 m/s. Its C/N₀ earns it a small weight, so
    # it barely moves the velocity — but the residual reported for it is the raw m/s
    # disagreement, not a weight-normalised one, so it still stands out as the outlier.
    states = gps_l1_states(0.0Hz)[1:6]
    reported = [with_state(s; cn0 = i == 1 ? 30.0dBHz : 48.0dBHz)
                for (i, s) in enumerate(states)]
    clean = calc_pvt(reported; kwargs...)
    broken_prn = reported[1].decoder.prn
    noisy = [i == 1 ? with_range_rate_error(s, 1.0) : s for (i, s) in enumerate(reported)]
    weighted = calc_pvt(noisy; kwargs...)

    rate_residual(pvt, prn) = pvt.sats[(:GPSL1CA, prn)].rate_residual
    @test all(info -> isfinite(ustrip(info.rate_residual)), values(weighted.sats))
    # The 1 m/s error shows up as ~1 m/s of residual on that satellite …
    @test abs(rate_residual(weighted, broken_prn) - rate_residual(clean, broken_prn)) >
          0.8m / s
    # … and the velocity barely moves, because that satellite is the one trusted least.
    unweighted = calc_pvt(strip_cn0(noisy); kwargs...)
    @test norm(weighted.velocity - clean.velocity) <
          0.5 * norm(unweighted.velocity - calc_pvt(strip_cn0(reported); kwargs...).velocity)
end

# The point of the exercise: satellites whose C/N₀ says they are noisy, carrying exactly
# that noise, degrade the fix less when the solve knows about it. Everything below is
# driven by one seeded RNG, so it is a fixed computation rather than a flaky one.
#
# Three marginal satellites among nine — the case the issue describes: a receiver admits
# a weak-but-tracked satellite (four are needed for a fix at all) and then has to live
# with its noise. Each estimator's error is measured against its *own* noiseless fix,
# because that isolates what is being compared: how much of the injected measurement
# noise reaches the solution. (The fixture measurements carry real-world errors of their
# own, so the two noiseless fixes are not identical — they differ by ~0.9 m — and neither
# is a ground truth.)
const WEAK_CN0 = 30.0dBHz
const STRONG_CN0 = 48.0dBHz
graded_cn0_states() = [
    with_state(s; cn0 = i ≤ 3 ? WEAK_CN0 : STRONG_CN0)
    for (i, s) in enumerate(gps_l1_states(0.0Hz))
]

@testset "a degraded-C/N₀ subset degrades a weighted fix less" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    reported = graded_cn0_states()

    # Noiseless reference fixes, and the σ the model assigns each satellite (elevation
    # included). The noise below is drawn with exactly the σ the solve is told to
    # expect, so what is measured is the weighting, not a model mismatch.
    truth = calc_pvt(reported; kwargs...)
    ols_truth = calc_pvt(strip_cn0(reported); kwargs...)
    σ = [ustrip(m, truth.sats[(:GPSL1CA, s.decoder.prn)].pseudorange_sigma)
         for s in reported]
    # The marginal satellites really are the noisy ones — every one of them noisier than
    # every strong one, ~7 m against ~2 m (the 30 dBHz / 48 dBHz code-noise gap).
    @test minimum(σ[1:3]) > maximum(σ[4:end])

    # One realisation, spelled out: +2σ on the three weak satellites and −2σ on the
    # strong ones, a jolt no estimator can average away — the weighted fix stays closer.
    errors = [i ≤ 3 ? 2σ[i] : -2σ[i] for i in eachindex(σ)]
    noisy = [with_range_error(s, e) for (s, e) in zip(reported, errors)]
    @test norm(calc_pvt(noisy; kwargs...).position - truth.position) <
          norm(calc_pvt(strip_cn0(noisy); kwargs...).position - ols_truth.position)

    # And as an RMS error over 150 realisations — the claim that matters, since a single
    # draw can favour either estimator. Weighted least squares is the minimum-variance
    # estimator when the weights are the true inverse variances (Gauss-Markov), so it
    # must win on average; here it does by ~21 % of the RMS error (6.5 m against 8.2 m),
    # and the velocity below by ~60 %.
    rng = MersenneTwister(20260805)
    ols_errors = Float64[]
    wls_errors = Float64[]
    for _ in 1:150
        draw = [with_range_error(s, σⱼ * randn(rng)) for (s, σⱼ) in zip(reported, σ)]
        push!(ols_errors,
            norm(calc_pvt(strip_cn0(draw); kwargs...).position - ols_truth.position))
        push!(wls_errors, norm(calc_pvt(draw; kwargs...).position - truth.position))
    end
    @test rms(wls_errors) < 0.85 * rms(ols_errors)
    # Not just an average: the weighted fix is the closer of the two in most individual
    # realisations as well. (A per-realisation *guarantee* is not on offer — a single
    # noise draw can favour either estimator, which is why the RMS above is the claim.)
    @test count(wls_errors .< ols_errors) > 0.6 * length(wls_errors)

    # The reported formal accuracy is a usable prediction of that spread, not merely an
    # ordering: the empirical 3D RMS error lands within 50 % of it.
    predicted = ustrip(m, truth.accuracy.position)
    @test 0.5 * predicted < rms(wls_errors) < 1.5 * predicted
    # The unweighted solve's own spread is *not* predicted by its accuracy figure, which
    # assumes a uniform nominal UERE — one more reason to report the weighted one.
    @test rms(ols_errors) > ustrip(m, ols_truth.accuracy.position) * 0.5
end

@testset "a degraded-C/N₀ subset degrades a weighted velocity less" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    reported = graded_cn0_states()
    truth = calc_pvt(reported; kwargs...)
    ols_truth = calc_pvt(strip_cn0(reported); kwargs...)
    # Doppler noise comes from the carrier loop, so the velocity solve is weighted by
    # `range_rate_variance` — its own model — rather than by the pseudorange σ.
    σv = [sqrt(PositionVelocityTime.range_rate_variance(s)) for s in reported]
    @test minimum(σv[1:3]) > maximum(σv[4:end])

    rng = MersenneTwister(20260806)
    ols_errors = Float64[]
    wls_errors = Float64[]
    for _ in 1:150
        draw = [with_range_rate_error(s, σⱼ * randn(rng)) for (s, σⱼ) in zip(reported, σv)]
        push!(ols_errors,
            norm(calc_pvt(strip_cn0(draw); kwargs...).velocity - ols_truth.velocity))
        push!(wls_errors, norm(calc_pvt(draw; kwargs...).velocity - truth.velocity))
    end
    @test rms(wls_errors) < 0.6 * rms(ols_errors)
    @test count(wls_errors .< ols_errors) > 0.75 * length(wls_errors)
end

# The weighted solve is only worth anything if it actually converges to the weighted
# optimum. It does from anywhere — but only because of the solver settings in
# `user_position`: with LsqFit's defaults the iteration stops at the first step shorter
# than `x_tol·‖ξ‖`, and `‖ξ‖` is dominated by the ~2e7 m clock bias, so a heavily damped
# first step would end the solve up to a metre short on an ordinary warm start.
@testset "the solve converges to the same fix from any starting point" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    reported = [with_state(s; cn0 = 40.0dBHz) for s in gps_l1_states(0.0Hz)]
    cold = calc_pvt(reported; kwargs...)

    # Warm-starting from the answer, and from previous fixes displaced by 1 m or 100 m,
    # all land on the same solution to well under a millimetre. A 100 km-stale start is
    # allowed a few centimetres: the a-priori weights are evaluated about the *reference*
    # position (as the atmospheric delays are), and 100 km moves the elevations enough to
    # move the σ they map to — so that is a slightly differently weighted epoch, not an
    # unconverged solve.
    for (offset, tolerance) in (
        [0.0, 0.0, 0.0] => 1e-4,
        [1.0, -1.0, 1.0] => 1e-4,
        [100.0, 0.0, -60.0] => 1e-4,
        [1e5, 1e5, -1e5] => 0.05,
    )
        stale = PVTSolution(;
            position = ECEF((Vector(cold.position) .+ offset)...),
            time_correction = cold.time_correction,
            reference_system = cold.reference_system,
        )
        warm = calc_pvt(reported, stale; kwargs...)
        @test norm(warm.position - cold.position) < tolerance
    end
end
