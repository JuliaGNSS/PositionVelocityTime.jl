# The modelled range rate must be the exact time derivative of the pseudorange the
# position solve models — Earth rotation over the travel time included — or the velocity
# solve inherits millimetres per second of bias from the model alone. The references
# below are derived from that pseudorange directly, by finite difference, so they share
# none of the model's algebra. `calc_pseudorange_rate` is that model for an estimator of its
# own (vector tracking), and the velocity solve fits the same model, so the tests also
# tie the two together.

const ONE_CLOCK_COLUMN = PositionVelocityTime.BiasColumns([1], 1, [0], 0)

# The pseudorange at receiver time `t`, as `calc_ρ_hat!` models it: the satellite taken
# at its transmit time `t_t = t − ρ/c` (found by fixed-point iteration), the receiver at
# `user_at(t)` with clock bias `bias_at(t)` (m), and the satellite clock offset
# `sat_clock_at(t_t)` (m) subtracted.
function modelled_pseudorange(sat_at, sat_clock_at, user_at, bias_at, t)
    c = PositionVelocityTime.SPEED_OF_LIGHT
    ρ = 0.0
    t_t = t
    for _ in 1:10
        t_t = t - ρ / c
        ρ = PositionVelocityTime.calc_ρ_hat!(
            zeros(1), [sat_at(t_t)], [user_at(t); bias_at(t)], ONE_CLOCK_COLUMN)[1] -
            sat_clock_at(t_t)
    end
    return ρ, t_t
end

central_difference(f, t, h) = (f(t + h) - f(t - h)) / 2h

# Measured inside a function — `@allocated` at top level would count the boxing of its own
# global arguments (Julia 1.10 does).
allocated_rate(args) = @allocated PositionVelocityTime.calc_pseudorange_rate(args...)

@testset "calc_pseudorange_rate is the derivative of the modelled pseudorange" begin
    c = PositionVelocityTime.SPEED_OF_LIGHT
    # A satellite in the equatorial plane, where the Earth's rotation acts in full, at
    # 3.2 km/s and receding from an equatorial user at a few hundred m/s.
    θ = deg2rad(20)
    sat_pos₀ = 2.656e7 * SVector(cos(θ), sin(θ), 0.0)
    sat_vel = 3200.0 * SVector(-sin(θ), cos(θ), 0.0)
    sat_clock_drift = 2e-11
    user_pos₀ = SVector(6.378e6, 0.0, 0.0)
    clock_bias, clock_drift = 1.0e3, 30.0

    for user_vel in (SVector(0.0, 0.0, 0.0), SVector(15.0, -40.0, 8.0))
        sat_at(t_t) = sat_pos₀ + sat_vel * t_t
        user_at(t) = user_pos₀ + user_vel * t
        pseudorange(t) = first(modelled_pseudorange(sat_at,
            t_t -> c * sat_clock_drift * t_t, user_at, t -> clock_bias + clock_drift * t, t))
        reference = central_difference(pseudorange, 0.0, 0.05)

        sat_pos = sat_at(last(modelled_pseudorange(sat_at,
            t_t -> c * sat_clock_drift * t_t, user_at, t -> clock_bias, 0.0)))
        rate = PositionVelocityTime.calc_pseudorange_rate(
            user_pos₀, user_vel, clock_drift, sat_pos, sat_vel, sat_clock_drift)
        @test 300 < dot(normalize(sat_pos - user_pos₀), sat_vel - user_vel) < 400
        @test rate ≈ reference atol = 1e-5

        # What the released model computed — the rotated line of sight against the
        # unrotated satellite velocity — is off by ~15 mm/s here, and rotating the
        # velocity alone still leaves ~0.6 mm/s; the test tells both apart.
        e = PositionVelocityTime.calc_line_of_sight(sat_pos, [user_pos₀; 0.0])
        τ = norm(sat_pos - user_pos₀) / c
        released = dot(e, sat_vel - user_vel) - c * sat_clock_drift + clock_drift
        rotated_sat_vel = PositionVelocityTime.earth_rotation(τ) * sat_vel
        rotated = dot(e, rotated_sat_vel - user_vel) - c * sat_clock_drift + clock_drift
        @test abs(released - reference) > 5e-3
        @test abs(rotated - reference) > 3e-4
    end

    # On the rotation axis the Earth's turn moves nothing, which leaves the light-time
    # factor alone: a satellite closing head-on at 1000 m/s reads 1000·(1 + 1000/c),
    # ~3.3 mm/s more than its speed.
    closing = PositionVelocityTime.calc_pseudorange_rate(SVector(0.0, 0.0, 6.4e6),
        SVector(0.0, 0.0, 0.0), 0.0, SVector(0.0, 0.0, 2.6e7), SVector(0.0, 0.0, -1000.0),
        0.0)
    @test closing ≈ -1000 * (1 + 1000 / c) rtol = 1e-14

    # Usable inside a consumer's own allocation-free loop.
    args = (user_pos₀, SVector(15.0, -40.0, 8.0), clock_drift, sat_pos₀, sat_vel,
        sat_clock_drift)
    allocated_rate(args)
    @test allocated_rate(args) == 0
end

# `state` re-timed so that its corrected transmit time is `transmit_time`, carrying
# `doppler` (Hz). Only the code phase and bit count move, as in `_precompile_states`;
# the carrier phase, which also counts into the transmit time, is kept.
function retimed_state(state, transmit_time, doppler)
    uncorrected = PositionVelocityTime.calc_uncorrected_time(state)
    # `correct_clock` is a slowly varying offset, so a few fixed-point steps invert it.
    for _ in 1:5
        uncorrected += transmit_time -
                       PositionVelocityTime.correct_clock(state.decoder, state.system, uncorrected)
    end
    rate = ustrip(Hz, get_data_frequency(state.decoder))
    chips = ustrip(Hz, get_code_frequency(state.system))
    offset = state.decoder.num_bits_after_valid_syncro_sequence / rate +
             mod(state.code_phase, chips / rate) / chips +
             uncorrected - PositionVelocityTime.calc_uncorrected_time(state)
    num_bits = floor(Int, offset * rate)
    SatelliteState(;
        decoder = GNSSDecoder.GNSSDecoderState(state.decoder;
            num_bits_after_valid_syncro_sequence = num_bits),
        system = state.system,
        code_phase = (offset - num_bits / rate) * chips,
        carrier_doppler = doppler * Hz,
        carrier_phase = state.carrier_phase,
    )
end

@testset "calc_pvt recovers the true velocity from exact observables: $name" for (name, states) in (
    ("GPS L1 C/A", gps_l1_states(0.0Hz)),
    ("Galileo E1B", galileo_e1b_states(0.0Hz)),
)
    c = PositionVelocityTime.SPEED_OF_LIGHT
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    # The fixture's own fix sets a realistic geometry, and its latest transmit time plus
    # a typical travel time the receiver clock reading at reception.
    fixture = calc_pvt(signal_group(states); kwargs...)
    user_pos = SVector{3,Float64}(fixture.position)
    # An aircraft's speed: the travel-time rate's user-velocity share of the Earth-rotation
    # term is ~2e-6·|v| per satellite, and leaving it out of the design costs ~0.1 mm/s of
    # velocity here — ten times the bounds below, which the solve meets to ~1e-6 m/s.
    user_vel = SVector(150.0, -200.0, 50.0)
    clock_bias, clock_drift = 0.0, 30.0
    reception = maximum(PositionVelocityTime.calc_corrected_time, states) + 0.075

    exact_states = map(states) do state
        orbit_at(t_t) = calc_satellite_position_and_velocity(state.decoder, t_t).position
        # The transmit time that puts this satellite on the true position…
        _, transmit_time = modelled_pseudorange(orbit_at, t_t -> 0.0, t -> user_pos,
            t -> clock_bias, reception)
        # …and the exact Doppler: the negated rate of the pseudorange there, with the
        # receiver moving and drifting and the satellite clock running at its drift.
        sat_clock_drift = PositionVelocityTime.calc_satellite_clock_drift(state.decoder,
            transmit_time)
        pseudorange(t) = first(modelled_pseudorange(orbit_at,
            t_t -> c * sat_clock_drift * (t_t - transmit_time),
            t -> user_pos + user_vel * (t - reception),
            t -> clock_bias + clock_drift * (t - reception), t))
        λ = c / ustrip(Hz, get_center_frequency(state.system))
        doppler = -central_difference(pseudorange, reception, 0.1) / λ
        # The public model at the true state is the exact rate, satellite by satellite.
        sat = calc_satellite_position_and_velocity(state.decoder, transmit_time)
        @test PositionVelocityTime.calc_pseudorange_rate(user_pos, user_vel, clock_drift,
            sat.position, sat.velocity, sat_clock_drift) ≈ -doppler * λ atol = 1e-5
        retimed_state(state, transmit_time, doppler)
    end

    pvt = calc_pvt(signal_group(exact_states); kwargs...)
    # The transmit times are seconds-of-week counts, resolved in Float64 to ~6e-11 s —
    # ~2 cm of range — so the position is exact only to millimetres.
    @test norm(pvt.position - user_pos) < 1e-2
    @test norm(pvt.velocity - user_vel) < 1e-5
    @test pvt.relative_clock_drift * c ≈ clock_drift atol = 1e-5
    @test maximum(info -> abs(ustrip(m / s, info.rate_residual)), pvt.sats) < 1e-5
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
    for row in measurement_rows(groups)
        modelled = PositionVelocityTime.calc_pseudorange_rate(user_pos, user_vel,
            pvt.relative_clock_drift * c, row.position, row.velocity, row.clock_drift)
        measured = -row.carrier_doppler * c / row.center_frequency
        residual = ustrip(m / s, pvt.sats[(row.signal_id, row.prn)].rate_residual)
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
    rows = measurement_rows(groups)
    for (user_vel, clock_drift) in ((SVector(150.0, -200.0, 50.0), 3000.0),
        (SVector(-600.0, 400.0, -300.0), -30000.0))
        truth = SVector(user_vel..., clock_drift)
        exact = [with_pseudorange_rate(row,
            PositionVelocityTime.calc_pseudorange_rate(user_pos, user_vel, clock_drift,
                row.position, row.velocity, row.clock_drift)) for row in rows]
        # From rest: hundreds of m/s and kilometres per second of drift in one step.
        solved, residuals =
            PositionVelocityTime.calc_user_velocity_and_clock_drift(exact, user_pos)
        @test norm(solved - truth) < 1e-6
        @test maximum(abs, residuals) < 1e-6
        # From a start a few m/s off, the previous epoch's solution, to rounding.
        start = truth + SVector(2.0, -1.5, 0.5, 0.3)
        solved, residuals =
            PositionVelocityTime.calc_user_velocity_and_clock_drift(exact, user_pos, start)
        @test norm(solved - truth) < 1e-10
        @test maximum(abs, residuals) < 1e-10
    end
end

# `calc_pvt` starts the velocity solve from the previous solution it trusts. The fix does
# not depend on that beyond the one step's precision, and a previous solution that is
# distrusted or not finite is not used.
@testset "calc_pvt warm-starts the velocity solve" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    groups = signal_group(gps_l1_states(0.0Hz))
    cold = calc_pvt(groups; kwargs...)
    warm = calc_pvt(groups, cold; kwargs...)
    rate_residuals(pvt) = [ustrip(m / s, info.rate_residual) for info in pvt.sats]
    @test norm(warm.velocity - cold.velocity) < 1e-9
    @test warm.relative_clock_drift ≈ cold.relative_clock_drift atol = 1e-9 / 3e8
    @test maximum(abs, rate_residuals(warm) - rate_residuals(cold)) < 1e-9

    # A previous solution with a non-finite velocity is the cold start.
    broken = PVTSolution(; position = cold.position, dop = cold.dop,
        velocity = ECEF(NaN, NaN, NaN), relative_clock_drift = NaN)
    @test calc_pvt(groups, broken; kwargs...).velocity ≈ cold.velocity atol = 1e-9
end
