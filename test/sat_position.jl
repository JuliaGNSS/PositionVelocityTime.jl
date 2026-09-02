# Satellite orbit propagation: the velocity half of `calc_satellite_position_and_velocity`.
#
# The rates are analytic derivatives of the position expressions, so the reference here is
# a central difference of the position the same function returns. That also pins the true
# anomaly rate at perigee and apogee, where the textbook `sin E / sin ν` form of ν̇ is a
# numerically singular `0/0`.
@testset "Satellite velocity" begin
    base = first(galileo_e1b_states(0.0Hz)).decoder
    t_0e = Float64(base.data.t_0e)
    # Mean anomaly is `M_0 + n·(t − t_0e)`, so `M_0` places the satellite anywhere on the
    # orbit at `t = t_0e`; `M = 0`/`M = π` put it exactly at perigee/apogee.
    at_mean_anomaly(M_0) = let d = GNSSDecoder.GalileoINAVData(base.data; M_0 = M_0)
        GNSSDecoder.GNSSDecoderState(base; data = d, raw_data = d)
    end
    # Central difference of the propagated position, in BigFloat: `t ≈ 1.3e5 s` against a
    # ~3.6e6 m coordinate leaves too little Float64 headroom for a metre-level difference.
    function numeric_velocity(decoder, t; h = big(1.0))
        pos(x) = PositionVelocityTime.calc_satellite_position(decoder, x)
        (pos(big(t) + h) - pos(big(t) - h)) / (2h)
    end

    @testset "matches a central difference of the position at M_0 = $M_0" for M_0 in
                                                                             (0.0, 0.7, Float64(π), 4.2)
        decoder = at_mean_anomaly(M_0)
        velocity = PositionVelocityTime.calc_satellite_position_and_velocity(decoder, t_0e).velocity
        @test all(isfinite, velocity)
        @test velocity ≈ Float64.(numeric_velocity(decoder, t_0e)) rtol = 1e-6
    end

    # Perigee (`E = ν = 0`) is an exact `0/0` for the `sin E / sin ν` form; apogee
    # (`E = ν = π`) is saved from a hard NaN only by `sin(π) != 0` in floating point, but
    # is still ill-conditioned there. Sweep across both to make sure no sample degrades.
    @testset "stays finite sweeping through perigee and apogee" begin
        for M_0 in (0.0, Float64(π))
            decoder = at_mean_anomaly(M_0)
            velocities = [
                PositionVelocityTime.calc_satellite_position_and_velocity(decoder, t).velocity
                for t in range(t_0e - 1, t_0e + 1; length = 101)
            ]
            @test all(v -> all(isfinite, v), velocities)
            @test all(v -> norm(v) > 1000, velocities)   # a MEO satellite moves ~3.7 km/s
        end
    end
end

# The pseudorange differencing must survive the week wrap: seconds-of-week counts
# do not all wrap at once — within one constellation the wrap sweeps through the
# transmit-time spread, and in a mixed solve the +14 s BDT scale alignment keeps
# BeiDou times above 604800 for the first 14 s of every GPS week. Unfolded, every
# difference across the wrap is off by a whole week (1.8e14 m). The absolute
# level of the folded ranges is irrelevant — the receiver clock column absorbs a
# common constant, exactly as it does the unmodelled travel time — so what the
# tests pin is the *differences*.
@testset "calc_pseudo_ranges folds the week wrap" begin
    c = PositionVelocityTime.SPEEDOFLIGHT

    # No wrap in sight: plain differences against the latest transmit time.
    # Metre-level tolerances throughout: differencing ~6e5-magnitude counts in
    # Float64 leaves ~1e-10 s ≈ 0.03 m of representation noise, and the failure
    # being pinned is fourteen orders of magnitude above it.
    pr, t_ref = PositionVelocityTime.calc_pseudo_ranges([410000.001, 410000.005])
    @test t_ref == 410000.005
    @test pr[2] == 0.0
    @test pr[1] ≈ 0.004 * c atol = 1.0

    # One constellation, wrap inside the transmit-time spread: the wrapped
    # satellite transmitted 2 ms *later*, so its pseudorange must come out
    # 2 ms · c below the unwrapped one — not a week above it.
    pr, t_ref = PositionVelocityTime.calc_pseudo_ranges([604799.999, 0.001])
    @test t_ref == 604799.999
    @test pr[1] - pr[2] ≈ 0.002 * c atol = 1.0

    # The mixed-constellation straddle: a BeiDou time on the GPS count
    # (SOW + 14) just past 604800 against a GPS time just past 0 — the same
    # instant, so the same pseudorange.
    pr, _ = PositionVelocityTime.calc_pseudo_ranges([604803.0, 3.0])
    @test pr[1] ≈ pr[2] atol = 1.0
end
