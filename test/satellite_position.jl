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
