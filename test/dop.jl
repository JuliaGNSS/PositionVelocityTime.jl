@testset "calc_DOP" begin
    # ECEF line-of-sight unit vector for a satellite at azimuth `az` (rad, from
    # North) and elevation `el` (rad), seen from `user`. Built by rotating the
    # local ENU direction back to ECEF via Geodesy, so the geometry matches what
    # the solver feeds calc_DOP (ECEF line-of-sight rows from calc_e).
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
    dop = PositionVelocityTime.calc_DOP(H_ecef, user)

    # Reference DOP computed directly in the local ENU frame (RTKLIB-style: rows
    # are the ENU line-of-sight components). calc_DOP must reproduce these from
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
    # DOP(-1, …) returned gracefully, not an exception.
    H_singular = repeat([1.0 0.0 0.0 1.0], 4)
    dop_singular = PositionVelocityTime.calc_DOP(H_singular, user)
    @test dop_singular.GDOP == -1
    @test dop_singular.PDOP == -1

    # Multi-GNSS: 3 position + 2 clock columns (3 sats per system). GDOP/HDOP are
    # independent of which clock is "primary"; TDOP reports the selected system's
    # clock variance, so it changes with primary_clock_index.
    azels6 = [azels; (deg2rad(300), deg2rad(35))]
    sys = [1, 1, 1, 2, 2, 2]
    H_multi = reduce(
        vcat,
        [
            [ecef_los(user, az, el)' (sys[i]==1 ? 1.0 : 0.0) (sys[i]==2 ? 1.0 : 0.0)] for
            (i, (az, el)) in enumerate(azels6)
        ],
    )
    d1 = PositionVelocityTime.calc_DOP(H_multi, user, 1)
    d2 = PositionVelocityTime.calc_DOP(H_multi, user, 2)
    @test d1.GDOP ≈ d2.GDOP
    @test d1.HDOP ≈ d2.HDOP
    @test d1.PDOP ≈ d2.PDOP
    @test !(d1.TDOP ≈ d2.TDOP)
    @test d1.GDOP > 0
end
