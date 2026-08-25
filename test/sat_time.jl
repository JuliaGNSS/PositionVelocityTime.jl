@testset "Week crossover" begin
    @test PositionVelocityTime.correct_week_crossovers(0) == 0

    @test PositionVelocityTime.correct_week_crossovers(350000) == 350000 - 604800

    @test PositionVelocityTime.correct_week_crossovers(-350000) == 604800 - 350000
end

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
        applied(x) = x - PositionVelocityTime.correct_clock(decoder, GalileoE1B(), x)
        h = big(1.0)
        (applied(big(t) + h) - applied(big(t) - h)) / (2h)
    end

    @testset "a_f2 = 0 (the usual broadcast): drift is a_f1" begin
        decoder = with(; a_f1 = 1e-11, a_f2 = 0.0, t_0c = 132000.0)
        t = 132100.0
        @test PositionVelocityTime.calc_satellite_clock_drift(decoder, t) == 1e-11
        @test PositionVelocityTime.calc_satellite_clock_drift(decoder, t) ≈
              Float64(numeric_drift(decoder, t)) atol = 1e-13
    end

    @testset "non-zero a_f2 late in the week: no spurious 2·a_f2·t_0c offset" begin
        # |a_f2| at the edge of the broadcast range (8 bit, LSB 2^-55 s/s²) with a clock
        # reference epoch near the end of the week — the worst case for the offset.
        a_f2, t_0c = 3.5e-15, 604000.0
        decoder = with(; a_f1 = 1e-11, a_f2 = a_f2, t_0c = t_0c)
        t = t_0c + 100.0

        drift = PositionVelocityTime.calc_satellite_clock_drift(decoder, t)
        @test drift ≈ 1e-11 + 2 * a_f2 * 100.0
        # Within the neglected relativistic-rate term, well below the 4.2e-9 s/s
        # (≈1.3 m/s of range rate) the `t`-based polynomial would have been off by.
        @test drift ≈ Float64(numeric_drift(decoder, t)) atol = 1e-13
        @test !isapprox(drift, 1e-11 + 2 * a_f2 * t, atol = 1e-12)
    end
end

# The carrier-phase term refines the transmit time by a fraction of a carrier cycle.
# `SatelliteState.carrier_phase` is in radians (`Tracking.get_carrier_phase`), so it must
# be converted to cycles before being divided by the centre frequency.
@testset "Carrier-phase transmit-time term treats carrier_phase as radians" begin
    system = GPSL1CA()
    decoder = GNSSDecoder.GNSSDecoderState(
        GNSSDecoderState(system, 1);
        num_bits_after_valid_syncro_sequence = 0,
        data = GNSSDecoder.GPSL1CAData(; TOW = 0),
        raw_data = GNSSDecoder.GPSL1CAData(; TOW = 0),
    )
    time_of(carrier_phase) = PositionVelocityTime.calc_uncorrected_time(
        SatelliteState(; decoder, system, code_phase = 0.0, carrier_doppler = 0.0Hz,
            carrier_phase),
    )

    for cycles in (0.25, 0.5, -0.5)
        @test time_of(cycles * 2π) ≈ cycles / ustrip(Hz, get_center_frequency(system))
    end
    # The wrapped carrier phase spans one cycle, so the term is bounded by half a carrier
    # wavelength of range — ~9.5 cm on L1, not the ~60 cm a radians-as-cycles reading gives.
    half_wavelength =
        0.5 * PositionVelocityTime.SPEEDOFLIGHT / ustrip(Hz, get_center_frequency(system))
    @test half_wavelength ≈ 0.0952 atol = 1e-4
    @test maximum(
        abs(time_of(φ) * PositionVelocityTime.SPEEDOFLIGHT) for
        φ in range(-Float64(π), Float64(π); length = 21)
    ) ≈ half_wavelength
end