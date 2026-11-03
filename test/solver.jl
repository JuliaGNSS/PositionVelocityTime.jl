@testset "PVT Galileo E1B with frequency offset of $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    galileo_e1b = GalileoE1B()
    states = galileo_e1b_states(freq_offset)

    # Fixture data was recorded on 2021-05-31. Pin `approximate_year` so
    # this test stays valid as the wall clock drifts past 2030 (when the
    # default `now()` anchor would start picking the wrong rollover cycle
    # for these archived fixtures).
    pvt = calc_pvt(
        signal_group(states);
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    expected_lla = LLA(; lat = 50.77885207231635, lon = 6.065568145321566, alt = 289.09688064471146)
    @test pvt.position ≈ ECEFfromLLA(wgs84)(expected_lla) rtol = 1e-8
    @test TAIEpoch(pvt.time) ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.1183385390904732)
    @test pvt.velocity ≈ ECEF(0.0, 0.0, 0.0) atol = 9
    # The modelled range rate carries the light-time factor `1 − ρ̇/c` and the rate of
    # the Earth-rotation angle (see `calc_user_velocity_and_clock_drift`), which moved this by ~3 mm/s
    # (~0.01 Hz) from the value of the plain line-of-sight model; the factor also scales
    # with the Doppler, so the frequency offsets now differ by a few mHz.
    @test pvt.relative_clock_drift * get_center_frequency(galileo_e1b) ≈ -(1675.64Hz + freq_offset) atol = 0.01Hz

    warm_pvt = calc_pvt(
        signal_group(states),
        pvt;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    # Warm and cold starts agree to the position solve's precision, not bit for bit: it
    # stops once its step is below `x_tol` relative to `‖ξ‖`, ~2 µm with the ~2e7 m of
    # clock bias in `ξ`, so which side of that the two land on depends on last-bit
    # rounding — 2 µm apart on one GitHub runner. Geodesy's `≈` on `LLA` defaults to 1 µm.
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt) atol = 1e-5
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end
@testset "PVT GPS L1 with frequency offset of $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    gpsl1 = GPSL1CA()
    states = gps_l1_states(freq_offset)

    # Fixture data was recorded on 2021-05-31. See the note on the
    # Galileo testset above for why we pin `approximate_year`.
    pvt = calc_pvt(
        signal_group(states);
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    expected_lla = LLA(; lat = 50.778851781017025, lon = 6.065622611231713, alt = 291.96260731963366)
    @test pvt.position ≈ ECEFfromLLA(wgs84)(expected_lla) rtol = 1e-8
    @test TAIEpoch(pvt.time) ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.1491024351271335)
    @test pvt.velocity ≈ ECEF(0.0, 0.0, 0.0) atol = 2.5
    # The modelled range rate carries the light-time factor `1 − ρ̇/c` and the rate of
    # the Earth-rotation angle (see `calc_user_velocity_and_clock_drift`), which moved this by ~3 mm/s
    # (~0.01 Hz) from the value of the plain line-of-sight model; the factor also scales
    # with the Doppler, so the frequency offsets now differ by a few mHz.
    @test pvt.relative_clock_drift * get_center_frequency(gpsl1) ≈ -(1632.61Hz + freq_offset) atol = 0.01Hz

    warm_pvt = calc_pvt(
        signal_group(states),
        pvt;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt) atol = 1e-5  # see the first testset
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end

# Marks a GPS L1 C/A SatelliteState's decoder unhealthy: the six health bits are the
# subframe-1 SV health word, and a set MSB (nav data bad) is what `is_sat_healthy`
# rejects.
function with_health_bad(state)
    d = state.decoder
    unhealthy(data) = GNSSDecoder.GPSL1CAData(data; sv_health = 0b100000)
    new_decoder = GNSSDecoder.GNSSDecoderState(
        d; raw_data = unhealthy(d.raw_data), data = unhealthy(d.data))
    SatelliteState(;
        decoder = new_decoder,
        system = state.system,
        code_phase = state.code_phase,
        carrier_doppler = state.carrier_doppler,
        carrier_phase = state.carrier_phase,
    )
end

@testset "PVT GPS L1 + Galileo E1B combined, frequency offset $freq_offset" for freq_offset in (0.0Hz, 500Hz, -1000Hz)
    groups = (
        gps = signal_group(gps_l1_states(freq_offset)),
        galileo = signal_group(galileo_e1b_states(freq_offset)),
    )

    # Independent inter-system-bias solve: 8 GPS + 5 Galileo satellites, with one
    # clock bias per system (3 + 2 unknowns).
    pvt = calc_pvt(groups; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    expected_pos = ECEF(4.0186793226897363e6, 427033.09443239716, 4.918251247796992e6)
    @test pvt.position ≈ expected_pos rtol = 1e-8
    # Pinned since the range rate is modelled as the derivative of the
    # Earth-rotation-corrected pseudorange, which moved it by ~9 mm/s;
    # its light-time factor scales with the Doppler, so the offsets differ by < 0.6 mm/s.
    @test pvt.velocity ≈ ECEF(-1.4427421149590258, 0.5397343562431747, -2.127408455072271) atol = 1e-3
    @test TAIEpoch(pvt.time) ≈ TAIEpoch(2021, 5, 31, 12, 53, 14.285)

    # One clock bias per GNSS; the reference (GPS, most satellites) bias is
    # reported as time_correction, with Galileo's offset relative to it.
    @test pvt.reference_system == GPST()
    @test Set(keys(pvt.inter_system_biases)) == Set([GST()])

    # The combined fix agrees with each single-system fix and with the
    # GPS-anchored time scale.
    gps_only = calc_pvt(signal_group(gps_l1_states(freq_offset)); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    gal_only = calc_pvt(signal_group(galileo_e1b_states(freq_offset)); approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test norm(pvt.position - gps_only.position) < 5
    @test norm(pvt.position - gal_only.position) < 5
    @test pvt.time ≈ gps_only.time

    # A single receiver clock drift is estimated for the whole constellation
    # (one oscillator), so the combined drift agrees with each single-system drift
    # — to well within the inter-system spread, not split per system.
    @test pvt.relative_clock_drift ≈ gps_only.relative_clock_drift atol = 5e-8
    @test pvt.relative_clock_drift ≈ gal_only.relative_clock_drift atol = 5e-8

    # Satellites are keyed by (signal, PRN), so both constellations are represented
    # and nothing is merged or dropped across systems: the combined count equals the
    # sum of the single-system counts (PRNs here don't overlap).
    @test all(k -> k isa Tuple{Symbol,Int}, keys(pvt.sats))
    @test Set(first(k) for k in keys(pvt.sats)) == Set([:GPSL1CA, :GalileoE1B])
    @test length(pvt.sats) ==
          length(gps_only.sats) + length(gal_only.sats)
    let (sig, prn) = first(keys(pvt.sats))
        @test get_sat_info(pvt, sig, prn) === pvt.sats[(sig, prn)]
    end
    # A satellite not used in the fix returns `nothing` rather than throwing.
    @test get_sat_info(pvt, :GPSL1CA, 999) === nothing

    # Per-satellite post-fit residuals are exposed via SatInfo (measured − modeled
    # pseudorange, metres); finite and small for a good fix.
    resids = [info.residual for info in values(pvt.sats)]
    @test length(resids) == length(pvt.sats)
    @test all(isfinite, resids)
    @test maximum(abs, resids) < 10m

    # The rate-domain counterpart from the Doppler solve (measured − modeled range
    # rate, m/s): finite, and metre-per-second-level for a stationary receiver.
    rate_resids = [info.rate_residual for info in values(pvt.sats)]
    @test length(rate_resids) == length(pvt.sats)
    @test all(isfinite, rate_resids)
    @test maximum(abs, rate_resids) < 10.0m/s

    warm_pvt = calc_pvt(groups, pvt; approximate_year = 2021, enable_ionospheric_correction = false, enable_tropospheric_correction = false)
    @test get_LLA(warm_pvt) ≈ get_LLA(pvt) atol = 1e-5  # see the first testset
    @test warm_pvt.time ≈ pvt.time
    @test warm_pvt.velocity ≈ pvt.velocity atol = 1e-6
end

# Every unsolvable epoch takes the same exit — `prev_pvt`, never a throw — whether the
# shortfall is in the raw satellite count, in how many of them are usable, or in the layout
# the constellation demands.
@testset "an unsolvable epoch returns prev_pvt instead of throwing" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    gps = gps_l1_states(0.0Hz)
    reference = calc_pvt(signal_group(gps); kwargs...)
    @test reference.position != ECEF(0, 0, 0)

    # Too few states to hold 3 position + 1 clock unknown, down to none at all: with no
    # previous solution the default (origin) one comes back unchanged. The signal is
    # passed explicitly because `n = 0` leaves no satellite to read it off — a group
    # carries its ranging signal precisely so that an empty one is still well formed.
    for n in 0:3
        scarce = signal_group(gps[1].system, gps[1:n])
        pvt = calc_pvt(scarce; kwargs...)
        @test pvt.position == ECEF(0, 0, 0)
        @test isnothing(pvt.time)
        @test isempty(pvt.sats)
        # Given a previous fix, that fix is what is carried forward.
        @test calc_pvt(scarce, reference; kwargs...) === reference
    end

    # An unhealthy satellite is dropped before the count is taken, so 4 states carrying
    # only 3 usable ones exit exactly as 3 states do.
    unhealthy = [gps[1:3]; [with_health_bad(gps[4])]]
    @test length(unhealthy) == 4
    @test calc_pvt(signal_group(unhealthy), reference; kwargs...) === reference
end

# Regression: the position solve must run to its optimum, not stop at the first step
# LsqFit's *relative* `x_tol` calls small. `‖ξ‖` is dominated by the clock bias — the
# ~2e7 m of common range the pseudoranges are referred to — which puts the default
# tolerance at a ~0.2 m step, and with the default Levenberg-Marquardt damping a first
# step is only a fraction of the distance to the optimum. A warm start a metre out
# therefore took one ~10 cm step and reported convergence, returning a fix most of a
# metre away from the one the very same measurements give from cold. See the solver
# settings in `user_position`.
@testset "a warm start converges to the same fix as a cold start" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    states = gps_l1_states(0.0Hz)
    cold = calc_pvt(signal_group(states); kwargs...)

    # Previous solutions displaced by 1 m, 30 m and 100 km. With no position-dependent
    # correction active the measurements define one optimum, so every start must reach
    # it — to well under a millimetre, not to within the starting error.
    for offset in ([1.0, -1.0, 1.0], [30.0, 0.0, -20.0], [1e5, 1e5, -1e5])
        stale = PVTSolution(;
            position = ECEF((Vector(cold.position) .+ offset)...),
            time_correction = cold.time_correction,
            reference_system = cold.reference_system,
        )
        warm = calc_pvt(signal_group(states), stale; kwargs...)
        @test norm(warm.position - cold.position) < 1e-4
        @test warm.time ≈ cold.time
        # The residuals are those of the converged fit, not of the stale position.
        @test maximum(abs, [info.residual for info in values(warm.sats)]) ≈
              maximum(abs, [info.residual for info in values(cold.sats)]) rtol = 1e-6
    end
end

# Regression: a warm start is an optimisation, not an anchor. Seeded far enough away,
# Gauss-Newton converges to a spurious root of the range equations instead of the real
# fix — and once the previous solution *is* that root, every later epoch is seeded in the
# same basin and the output freezes there (measured on TEX-CUP: 2.7e8 m out, PDOP 587,
# byte-stable for 355 s, while a cold start of the same measurements solved to metres).
# Spurious roots flag themselves — the satellites collapse into a narrow cone there, so
# the fix's own DOP explodes while its residuals stay small — so `calc_pvt` refuses to
# seed from a previous solution whose DOP is implausible and solves the epoch from cold.
# Before the gate existed the poisoned position came straight back out of every epoch.
# Distrust affects the seed only: an unsolvable epoch still returns `prev_pvt` as
# passed, like every other unsolvable epoch does.
@testset "a distrusted previous fix is discarded rather than re-seeded" begin
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    states = gps_l1_states(0.0Hz)
    good = calc_pvt(signal_group(states); kwargs...)

    # A previous solution 2.7e8 m out, carrying the DOP such a fix earns (the observed
    # incident read PDOP 587) — exactly what feeding a locked-in epoch's output back in
    # as the next seed looks like.
    poisoned = PVTSolution(;
        position = ECEF((Vector(good.position) .+ [2.7e8, 0.0, 0.0])...),
        time_correction = good.time_correction,
        reference_system = good.reference_system,
        dop = PositionVelocityTime.DOP(613.0, 587.0, 400.0, 430.0, 176.6),
    )
    recovered = calc_pvt(signal_group(states), poisoned; kwargs...)

    # The seed is discarded whole, so the epoch is indistinguishable from one solved
    # without any previous solution at all — not merely closer to it.
    @test norm(recovered.position - good.position) < 1e-3
    @test recovered.time ≈ good.time
    @test recovered.velocity ≈ good.velocity atol = 1e-6
    # The reported geometry is the cold solve's, not the seed's.
    @test !isnothing(recovered.dop)
    @test recovered.dop.PDOP ≈ good.dop.PDOP
    @test recovered.dop.PDOP < PositionVelocityTime.MAX_TRUSTED_WARM_START_PDOP
    # Nothing of the poisoned seed survives into the reported satellite set.
    @test maximum(abs, [info.residual for info in values(recovered.sats)]) ≈
          maximum(abs, [info.residual for info in values(good.sats)]) rtol = 1e-6

    # An unsolvable epoch (too few satellites) returns the previous solution exactly
    # as passed, distrusted or not — the gate withholds it from the solve, not from
    # the caller.
    @test calc_pvt(signal_group(states[1:3]), poisoned; kwargs...) === poisoned
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
