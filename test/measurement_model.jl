# The modelled range rate must be the exact time derivative of the pseudorange the
# position solve models — Earth rotation over the travel time included — or the velocity
# solve inherits millimetres per second of bias from the model alone. The references
# below are derived from that pseudorange directly, by finite difference, so they share
# none of the model's algebra. `calc_pseudorange_rate` is that model for an estimator of its
# own (vector tracking), and the velocity solve fits the same model, so the tests also
# tie the two together.

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

@testset "calc_line_of_sight points receiver → satellite" begin
    # The textbook convention (Misra & Enge's 1ₖ, RTKLIB's `e`), pinned directly
    # because every H-based assertion is insensitive to it: `calc_H!` negates the
    # vector into the pseudorange's position partial, so a flipped convention with a
    # matching flip there would pass the whole suite unnoticed.
    user = [3.9074e6, 3.0684e5, 5.0150e6]
    sat = user + 2.0e7 * normalize([0.3, 0.1, 0.95])
    e = PositionVelocityTime.calc_line_of_sight(sat, [user; 0.0])
    @test norm(e) ≈ 1.0
    # Up to the ~30 m Sagnac rotation of the satellite position.
    @test dot(e, normalize(sat - user)) > 1 - 1e-10
    # And the design row carries its negative.
    H = PositionVelocityTime.calc_H(reshape(sat, 3, 1), [user; 0.0],
        PositionVelocityTime.BiasColumns([1], 1, [0], 0))
    @test H[1, 1:3] ≈ -e
end
