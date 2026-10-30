@testset "PVT satellites keyed by (signal, PRN): same PRN across signals coexists" begin
    gps = gps_l1_states(0.0Hz)
    gal = galileo_e1b_states(0.0Hz)
    healthy(s) = PositionVelocityTime.is_sat_healthy(s.decoder)

    # Pick a healthy satellite from each constellation and force the Galileo one to
    # share the GPS satellite's PRN, so the only thing telling them apart is the
    # signal tag in the key.
    gps_h = first(s for s in gps if healthy(s))
    gal_h_idx = findfirst(healthy, gal)
    shared_prn = gps_h.decoder.prn
    collided = with_prn(gal[gal_h_idx], shared_prn)
    groups = (
        signal_group(gps),
        signal_group([collided]),
        signal_group(gal[setdiff(1:length(gal), gal_h_idx)]),
    )

    baseline = calc_pvt((signal_group(gps), signal_group(gal)); approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    pvt = calc_pvt(groups; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false)

    # Both the GPS and the Galileo satellite with the shared PRN survive as
    # distinct entries — with a bare-Int PRN key one would have overwritten the
    # other and the used-satellite count would have dropped.
    @test haskey(pvt.observations, (:GPSL1CA, shared_prn))
    @test haskey(pvt.observations, (:GalileoE1B, shared_prn))
    @test pvt.observations[(:GPSL1CA, shared_prn)] !==
          pvt.observations[(:GalileoE1B, shared_prn)]
    @test length(pvt.observations) == length(baseline.observations)
end

# Course over ground: azimuth of the velocity vector in the local ENU frame, degrees
# clockwise from North. The recorded fixtures are stationary, so exercise the
# computation directly with known ENU velocities rotated into ECEF at a fixed position.
@testset "calc_course_over_ground is the ENU velocity azimuth clockwise from North" begin
    user = ECEF(ECEFfromLLA(wgs84)(LLA(50.1, 8.7, 120.0)))
    ecef_from_enu = ECEFfromENU(user, wgs84)
    # ENU velocity → ECEF velocity vector: rotate the tip and drop the origin translation.
    ecef_velocity(e, n, u) = ECEF(ecef_from_enu(ENU(e, n, u))) - user
    cog(e, n, u) = PositionVelocityTime.calc_course_over_ground(user, ecef_velocity(e, n, u))

    @test cog(0.0, 10.0, 0.0) ≈ 0° atol = 1e-6°      # North
    @test cog(10.0, 0.0, 0.0) ≈ 90° atol = 1e-6°     # East
    @test cog(0.0, -10.0, 0.0) ≈ 180° atol = 1e-6°   # South
    @test cog(-10.0, 0.0, 0.0) ≈ 270° atol = 1e-6°   # West
    @test cog(10.0, 10.0, 0.0) ≈ 45° atol = 1e-6°    # North-East
    # A stationary receiver has no horizontal velocity ⇒ course defined as 0°.
    @test PositionVelocityTime.calc_course_over_ground(user, ECEF(0.0, 0.0, 0.0)) == 0.0°
    # The vertical component does not affect the horizontal course.
    @test cog(10.0, 0.0, 7.0) ≈ cog(10.0, 0.0, 0.0)
    # Always wrapped to [0, 360)°.
    @test 0° ≤ cog(-3.0, -4.0, 0.0) < 360°
end
