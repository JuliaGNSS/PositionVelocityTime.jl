# `calc_satellite_clock_drift` must be the time derivative of the clock correction that
# `correct_clock` applies, so it has to be evaluated about the clock reference epoch
# `t_0c` — not about `t`, which would offset the `a_f2` term by `2·a_f2·t_0c`.
@testset "Satellite clock drift is the derivative of the clock correction" begin
    base = first(galileo_e1b_states(0.0Hz)).decoder
    with(; kwargs...) = let d = GNSSDecoder.GalileoINAVData(base.data; kwargs...)
        GNSSDecoder.GNSSDecoderState(base; data = d, raw_data = d)
    end
    # The applied correction, differentiated in BigFloat: the Float64 round trip through
    # `correct_clock` (t ≈ 6e5 s, correction ≈ 2e-4 s) loses the ~1e-13 s difference.
    function numeric_drift(decoder, t)
        applied(x) = x - PositionVelocityTime.correct_clock(ephemeris(decoder), GNSSDecoder.get_group_delay(decoder, GalileoE1B()), x)
        h = big(1.0)
        (applied(big(t) + h) - applied(big(t) - h)) / (2h)
    end

    @testset "a_f2 = 0 (the usual broadcast): drift is a_f1" begin
        decoder = with(; a_f1 = 1e-11, a_f2 = 0.0, t_0c = 132000.0)
        t = 132100.0
        @test PositionVelocityTime.calc_satellite_clock_drift(ephemeris(decoder), t) == 1e-11
        @test PositionVelocityTime.calc_satellite_clock_drift(ephemeris(decoder), t) ≈
              Float64(numeric_drift(decoder, t)) atol = 1e-13
    end

    @testset "non-zero a_f2 late in the week: no spurious 2·a_f2·t_0c offset" begin
        # |a_f2| at the edge of the broadcast range (8 bit, LSB 2^-55 s/s²) with a clock
        # reference epoch near the end of the week — the worst case for the offset.
        a_f2, t_0c = 3.5e-15, 604000.0
        decoder = with(; a_f1 = 1e-11, a_f2 = a_f2, t_0c = t_0c)
        t = t_0c + 100.0

        drift = PositionVelocityTime.calc_satellite_clock_drift(ephemeris(decoder), t)
        @test drift ≈ 1e-11 + 2 * a_f2 * 100.0
        # Within the neglected relativistic-rate term, well below the 4.2e-9 s/s
        # (≈1.3 m/s of range rate) the `t`-based polynomial would have been off by.
        @test drift ≈ Float64(numeric_drift(decoder, t)) atol = 1e-13
        @test !isapprox(drift, 1e-11 + 2 * a_f2 * t, atol = 1e-12)
    end
end

# The clock polynomial's argument is a time *difference*, so it must be folded modulo the
# week like every other one here. With `t_0c` late in the week and a transmit time just
# past the rollover, the raw `t - t_0c` reads −603900 s where the elapsed time is +900 s:
# a typical `a_f1` alone turns a ~9 ns correction into a ~6 µs one — 1.8 km of range —
# with the sign flipped, and the `a_f2` term is worse still.
@testset "Clock polynomial unwraps a week crossover" begin
    base = first(galileo_e1b_states(0.0Hz)).decoder
    with(; kwargs...) = let d = GNSSDecoder.GalileoINAVData(base.data; kwargs...)
        GNSSDecoder.GNSSDecoderState(base; data = d, raw_data = d)
    end

    a_f1, a_f2 = 1e-11, 3.5e-15
    t_0c, t = 604000.0, 100.0
    wrapped = t - t_0c + GNSSSignals.SECONDS_PER_WEEK
    @test wrapped == 900.0

    # The correction `correct_clock` applies, differenced against an otherwise identical
    # decoder with the drift terms zeroed: `a_f0`, the relativistic term and the group
    # delay cancel, leaving exactly `a_f1·Δt + a_f2·Δt²` — the part that sees the wrap.
    applied(decoder) = t - PositionVelocityTime.correct_clock(ephemeris(decoder), GNSSDecoder.get_group_delay(decoder, GalileoE1B()), t)
    drift_terms =
        applied(with(; a_f1, a_f2, t_0c)) - applied(with(; a_f1 = 0.0, a_f2 = 0.0, t_0c))
    @test drift_terms ≈ a_f1 * wrapped + a_f2 * wrapped^2 atol = 1e-13
    # The unwrapped argument would put a millisecond — hundreds of kilometres — here.
    @test !isapprox(drift_terms, a_f1 * (t - t_0c) + a_f2 * (t - t_0c)^2, atol = 1e-9)

    # `calc_satellite_clock_drift` is the derivative of that same polynomial about that
    # same epoch, so it has to unwrap identically.
    drift = PositionVelocityTime.calc_satellite_clock_drift(ephemeris(with(; a_f1, a_f2, t_0c)), t)
    @test drift ≈ a_f1 + 2 * a_f2 * wrapped
    @test !isapprox(drift, a_f1 + 2 * a_f2 * (t - t_0c), rtol = 1e-3)
end

# The carrier-phase term refines the transmit time by a fraction of a carrier cycle.
# `SatelliteState.carrier_phase` is in cycles (`Tracking.get_carrier_phase`), so it is
# divided by the centre frequency as it is.
@testset "Carrier-phase transmit-time term treats carrier_phase as cycles" begin
    signal = GPSL1CA()
    decoder = GNSSDecoder.GNSSDecoderState(
        GNSSDecoderState(signal, 1);
        num_bits_since_sync = 0,
        data = GNSSDecoder.GPSL1CAData(; TOW = 0),
        raw_data = GNSSDecoder.GPSL1CAData(; TOW = 0),
    )
    time_of(carrier_phase) = PositionVelocityTime.calc_uncorrected_transmit_time(
        SatelliteState(; decoder, signal, code_phase = 0.0, carrier_doppler = 0.0Hz,
            carrier_phase),
    )

    for cycles in (0.25, 0.5, -0.5)
        @test time_of(cycles) ≈ cycles / ustrip(Hz, get_center_frequency(signal))
    end
    # The wrapped carrier phase spans one cycle, so the term is bounded by half a carrier
    # wavelength of range — ~9.5 cm on L1, not the ~60 cm reading radians as cycles gives.
    half_wavelength =
        0.5 * PositionVelocityTime.SPEED_OF_LIGHT / ustrip(Hz, get_center_frequency(signal))
    @test half_wavelength ≈ 0.0952 atol = 1e-4
    @test maximum(
        abs(time_of(φ) * PositionVelocityTime.SPEED_OF_LIGHT) for
        φ in range(-0.5, 0.5; length = 21)
    ) ≈ half_wavelength
end
# The pseudorange differencing must survive the week wrap: seconds-of-week counts
# do not all wrap at once — within one constellation the wrap sweeps through the
# transmit-time spread, and in a mixed solve the +14 s BDT scale alignment keeps
# BeiDou times above 604800 for the first 14 s of every GPS week. Unfolded, every
# difference across the wrap is off by a whole week (1.8e14 m).
@testset "Transmit-time differences fold the week wrap" begin
    # Metre-level tolerances: differencing ~6e5-magnitude counts in Float64 leaves
    # ~1e-10 s ≈ 0.03 m of representation noise, and the failure being pinned is
    # fourteen orders of magnitude above it.
    c = PositionVelocityTime.SPEED_OF_LIGHT
    trail(reference, t) = PositionVelocityTime.calc_pseudorange(reference, t)

    # No wrap in sight: the plain difference.
    @test trail(410000.005, 410000.001) ≈ 0.004 * c atol = 1.0

    # One constellation, wrap inside the transmit-time spread: the wrapped satellite
    # transmitted 2 ms *later*, so it trails by −2 ms · c, not a week.
    @test trail(604799.999, 0.001) ≈ -0.002 * c atol = 1.0
    @test trail(0.001, 604799.999) ≈ 0.002 * c atol = 1.0

    # The mixed-constellation straddle: a BeiDou time on the GPS count (SOW + 14) just
    # past 604800 against a GPS time just past 0 — the same instant.
    @test trail(3.0, 604803.0) ≈ 0.0 atol = 1.0
    # The same, from the BeiDou satellite's own count and its leap-second correction.
    @test PositionVelocityTime.calc_pseudorange(3.0, 604789.0, 14.0) ≈ 0.0 atol = 1.0
end
