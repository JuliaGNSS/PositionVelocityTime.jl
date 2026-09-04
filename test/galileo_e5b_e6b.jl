# Galileo E5b (I/NAV) and E6-B (C/NAV / HAS) support.
#
# These two additions are opposites. E5b broadcasts the *same* I/NAV message as
# E1-B — same container, same ephemeris, same clock — on a different carrier, so
# supporting it is almost entirely about the band: a second, independent
# measurement of the same satellite, with its own group-delay scaling and its own
# inter-frequency bias. E6-B broadcasts no ephemeris at all (its C/NAV carries HAS
# corrections *to other satellites*), so supporting it is about making sure it is
# recognised and excluded rather than propagated.

@testset "Galileo E5b (I/NAV) support" begin
    e1b_state = galileo_e1b_states(0.0Hz)[1]
    e1b = e1b_state.decoder
    d = e1b.data

    # E5b-I decodes into the very same `GalileoINAVData` as E1-B (GNSSDecoder splits
    # the I/NAV core out of `e1b.jl` precisely so the two can share it), so an E5b
    # decoder state is the E1B one with a different constants tag.
    as_e5b(data = d) = GNSSDecoder.GNSSDecoderState(
        GNSSDecoder.GalileoE5bDecoderState(e1b.prn);
        data = data,
        raw_data = data,
        num_bits_after_valid_syncro_sequence = e1b.num_bits_after_valid_syncro_sequence,
    )
    e5b = as_e5b()

    @testset "classification: own band and signal, shared time system and week" begin
        @test GNSSSignals.get_time_system(GalileoE5bI()) == GST()
        @test GNSSSignals.get_signal_id(GalileoE5bI()) == :GalileoE5bI
        # Its own band — this is what earns E5b a separate inter-frequency bias, and
        # what keeps it from being grouped with E1 the way E1B and E1C are.
        @test GNSSSignals.get_band_id(GalileoE5bI()) == :E5b
        @test GNSSSignals.get_band_id(GalileoE5bI()) != GNSSSignals.get_band_id(GalileoE1B())
        @test PositionVelocityTime.get_week(e5b) == d.WN
        @test PositionVelocityTime.system_start_epoch(GalileoE5bI()) ==
              PositionVelocityTime.system_start_epoch(GalileoE1B())
        # No GEO branch for Galileo.
        @test !PositionVelocityTime.is_geo_orbit(e5b)
    end

    @testset "orbit and clock are bit-for-bit the E1B ones" begin
        # Same container, same constants values, so the shared propagator must not
        # merely agree to a tolerance — it must return the identical numbers.
        for t in (132000.0, 132769.0, 135000.0)
            @test PositionVelocityTime.calc_satellite_position_and_velocity(e5b, t) ==
                  PositionVelocityTime.calc_satellite_position_and_velocity(e1b, t)
            @test PositionVelocityTime.calc_satellite_clock_drift(e5b, t) ==
                  PositionVelocityTime.calc_satellite_clock_drift(e1b, t)
        end
    end

    @testset "group delay: E1–E5b BGD scaled to the E5b carrier" begin
        # I/NAV's clock is referred to the E1/E5b ionosphere-free combination, so an
        # E1 range applies BGD(E1,E5b) as broadcast and an E5b range applies
        # `(f_E1/f_E5b)²·BGD(E1,E5b)` (OS SIS ICD Issue 2.2 §5.1.5).
        scaling = PositionVelocityTime.galileo_group_delay_scaling(GalileoE5bI())
        @test scaling ≈ (1575.42 / 1207.14)^2
        @test scaling > 1
        for signal in (GalileoE5bI(), GalileoE5bQ())
            @test PositionVelocityTime.correct_by_group_delay(e5b, signal, 100.0) ≈
                  100.0 - scaling * d.BGD_E1_E5b
        end
        # The same I/NAV message on its other band applies the bare BGD: the scaling is
        # a property of the ranging band, not of the message.
        @test PositionVelocityTime.correct_by_group_delay(e1b, GalileoE1B(), 100.0) ≈
              100.0 - d.BGD_E1_E5b
        # I/NAV also broadcasts BGD(E1,E5a), but no method reads it: an E5a range is
        # generated on E5a-I and so carries an F/NAV decoder, never an I/NAV one.
        # Cross-band pairings cannot arise, and are an error rather than a
        # plausible-looking wrong number.
        @test_throws MethodError PositionVelocityTime.correct_by_group_delay(
            e5b, GalileoE5aI(), 100.0)
        # E6 carries no BGD at all, so there is no right answer to give for it either.
        @test_throws MethodError PositionVelocityTime.correct_by_group_delay(
            e5b, GalileoE6B(), 100.0)
    end

    @testset "health reports the E5b facet, not the E1-B/C one" begin
        # Word type 5 carries both facets, so an I/NAV decoder always has both and the
        # signal layer picks. A satellite withdrawn on E1 but serving on E5b is usable
        # here and not there — the one decode-level difference between the two.
        e5b_only = GNSSDecoder.GalileoINAVData(
            d;
            E1B_SHS = GNSSDecoder.signal_out_of_service,
            E1B_DVS = GNSSDecoder.working_without_guarantee,
        )
        @test GNSSDecoder.is_sat_healthy(as_e5b(e5b_only))
        @test !GNSSDecoder.is_sat_healthy(
            GNSSDecoder.GNSSDecoderState(e1b; data = e5b_only, raw_data = e5b_only),
        )
    end

    @testset "calc_pvt: E1B + E5b is a two-band fix with its own IFB" begin
        kw = (;
            approximate_year = 2021,
            enable_ionospheric_correction = false,
            enable_tropospheric_correction = false,
        )
        e1b_states = galileo_e1b_states(0.0Hz)
        reference = calc_pvt(e1b_states; kw...)
        @test length(reference.sats) >= 4

        # Split a transmit time into the decoder's bit count plus a sub-symbol code
        # phase, the way `calc_uncorrected_time` reads them back (it reduces the code
        # phase modulo one data symbol, so the whole-symbol part must go in the count).
        function bits_and_code_phase(system, tow, target)
            datafreq = Float64(GNSSSignals.get_data_frequency(system) / Hz)
            codefreq = Float64(GNSSSignals.get_code_frequency(system) / Hz)
            elapsed = target - tow
            num_bits = floor(Int, elapsed * datafreq)
            (num_bits, (elapsed - num_bits / datafreq) * codefreq)
        end

        # An E5b copy of each satellite, transmit-time-consistent with its E1B original.
        # Both range off one I/NAV message, so the only thing that differs is the band:
        # its carrier (hence the code-phase/bit bookkeeping) and its group-delay scaling.
        function as_e5b_state(state)
            decoder = state.decoder
            # The E5b range applies a larger share of the BGD than the E1 range
            # (`s·BGD` against `BGD`), and the corrected transmit time carries `+BGD`,
            # so a copy that lands on the *same* corrected time must start `(1−s)·BGD`
            # away. Without this the test would read that difference out as a spurious
            # inter-frequency bias instead of the ~0 it is checking for.
            target =
                PositionVelocityTime.calc_uncorrected_time(state) +
                (1 - PositionVelocityTime.galileo_group_delay_scaling(GalileoE5bI())) *
                decoder.data.BGD_E1_E5b
            num_bits, code_phase =
                bits_and_code_phase(GalileoE5bI(), decoder.data.TOW, target)
            SatelliteState(;
                decoder = GNSSDecoder.GNSSDecoderState(
                    GNSSDecoder.GalileoE5bDecoderState(decoder.prn);
                    data = decoder.data,
                    raw_data = decoder.data,
                    num_bits_after_valid_syncro_sequence = num_bits,
                ),
                system = GalileoE5bI(),
                code_phase = code_phase,
                carrier_doppler = 0.0Hz,
                carrier_phase = 0.0,
            )
        end

        both = calc_pvt([e1b_states; map(as_e5b_state, e1b_states)]; kw...)
        @test length(both.sats) == 2 * length(reference.sats)
        # One time system (both are Galileo) but two bands, so exactly one IFB column,
        # on E5b against the L1 reference.
        @test both.reference_system == GST()
        @test collect(keys(both.inter_frequency_biases)) == [:E5b]
        @test both.inter_frequency_biases[:E5b].reference == :L1
        @test abs(both.inter_frequency_biases[:E5b].value) < 1e-2m
        @test norm(both.position - reference.position) < 1e-2
    end
end

@testset "Galileo E6-B (C/NAV) is recognised and excluded" begin
    kw = (;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )

    # A fully operational E6-B decoder: HAS status nominal, nothing wrong with it. It
    # still cannot contribute a pseudorange, because C/NAV broadcasts corrections to
    # *other* satellites' navigation data and no ephemeris of its own.
    e6b_data = GNSSDecoder.GalileoE6BData(; HAS_status = GNSSDecoder.has_operational_mode)
    e6b = GNSSDecoder.GNSSDecoderState(
        GNSSDecoder.GalileoE6BDecoderState(4);
        data = e6b_data,
        raw_data = e6b_data,
        num_bits_after_valid_syncro_sequence = 16,
    )
    e6b_state = SatelliteState(;
        decoder = e6b,
        system = GalileoE6B(),
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
        carrier_phase = 0.0,
    )

    @testset "never positioning-complete, whatever the service status" begin
        # This is what keeps it out of the solve, and it is a property of the message,
        # not of how much of it has been decoded.
        @test !GNSSDecoder.is_decoding_completed_for_positioning(e6b)
        @test GNSSDecoder.is_sat_healthy(e6b)          # the HAS *service* is nominal
        for status in (
            GNSSDecoder.has_test_mode,
            GNSSDecoder.has_status_reserved,
            GNSSDecoder.has_do_not_use,
        )
            data = GNSSDecoder.GalileoE6BData(; HAS_status = status)
            @test !GNSSDecoder.is_decoding_completed_for_positioning(
                GNSSDecoder.GNSSDecoderState(e6b; data = data, raw_data = data),
            )
        end
    end

    @testset "answers the per-decoder queries instead of erroring on absent fields" begin
        # C/NAV has no week number, no time-offset record and no ionospheric
        # coefficients, so the accessors that read those must not be dispatched onto
        # it. `GalileoE6BData` is an `AbstractGalileoData` but deliberately not an
        # `AbstractGalileoEphemerisData`, and keying on the narrower supertype is what
        # makes these answer "no" rather than raise a `FieldError` on a field that was
        # never broadcast.
        @test PositionVelocityTime.time_offset_available(e6b, GPST()) == false
        @test PositionVelocityTime.ntcm_g_params(e6b) === nothing
        @test PositionVelocityTime.klobuchar_params(e6b) === nothing
        # A week number genuinely does not exist for C/NAV, so this is a MethodError
        # rather than a wrong answer — the honest outcome for a question that has none.
        @test_throws MethodError PositionVelocityTime.get_week(e6b)
    end

    @testset "calc_pvt ignores it without touching its missing fields" begin
        ranging = galileo_e1b_states(0.0Hz)
        reference = calc_pvt(ranging; kw...)
        # Appending the E6-B state must change nothing: the health filter drops it
        # before anything reads `t_0e`, `a_f0` or any other field it does not have.
        # (Were it not filtered, this call would throw a `FieldError`, not misfix.)
        with_e6b = calc_pvt([ranging; [e6b_state]]; kw...)
        @test with_e6b.position == reference.position
        @test with_e6b.velocity == reference.velocity
        @test with_e6b.time == reference.time
        @test keys(with_e6b.sats) == keys(reference.sats)
        # E6 contributes no band either, so it cannot invent an inter-frequency bias.
        @test !haskey(with_e6b.inter_frequency_biases, :E6)
    end

    @testset "an E6-B-only epoch is unsolvable rather than an error" begin
        # `calc_pvt` returns `prev_pvt` unchanged when an epoch cannot be solved, so an
        # epoch of nothing but E6-B satellites is the empty solution, not a throw.
        unsolved = calc_pvt([e6b_state]; kw...)
        @test unsolved.position == PVTSolution().position
        @test isnothing(unsolved.time)
        @test isempty(unsolved.sats)
    end
end

# ---------------------------------------------------------------------------
# The group-delay method table, checked as a whole rather than case by case.
#
# `correct_by_group_delay` dispatches on (decoded message, ranging signal), and a
# missing cell is silent: it shows only as a MethodError from inside `calc_pvt`, on
# whatever pairing a user happens to bring. That is how four GPS cells stayed empty
# for years — every test paired a decoder with a signal of its own band, which is the
# only pairing that occurs, so nothing ever probed the rest of the table.
#
# The rule the table follows: a method exists for a ranging signal iff that signal
# shares a band with the data component the message was decoded from AND the message
# carries a correction term for it. The first half holds because a pseudorange is
# generated on the signal it was tracked on and the ephemeris comes from that band's
# data component, so cross-band and cross-constellation pairings cannot arise. The
# second half is why GPS LNAV serves only C/A: it carries `T_GD` and no ISCs, so it
# has nothing to say about the L1C components that share its band.
@testset "group-delay methods cover exactly the reachable pairings" begin
    # Per message: the ranging signals it can be asked about, and why the rest of its
    # own constellation is excluded.
    expected = [
        # LNAV: T_GD only, no ISCs -> its own signal, not the L1C components.
        GNSSDecoder.GPSL1CAData => [GPSL1CA()],
        # CNAV: ISC_L5I5/L5Q5/L2C. One message on two bands, so both are its own.
        # It also carries ISC_L1CA, which is cross-band and therefore unreachable.
        GNSSDecoder.GPSCNAVData => [GPSL5I(), GPSL5Q(), GPSL2CM(), GPSL2CL()],
        # CNAV-2: the full L1 ISC set. Its ISC_L2C/L5* are cross-band.
        GNSSDecoder.GPSL1C_DData => [GPSL1C_D(), GPSL1C_P(), GPSL1CA()],
        # I/NAV: one message on E1 and E5b. Carries BGD(E1,E5a) too, cross-band.
        GNSSDecoder.GalileoINAVData => [GalileoE1B(), GalileoE1C(),
            GalileoE1B_BOC11(), GalileoE1C_BOC11(), GalileoE5bI(), GalileoE5bQ()],
        # F/NAV: E5a only.
        GNSSDecoder.GalileoE5aData => [GalileoE5aI(), GalileoE5aQ(), GalileoE5aQP()],
        # D1/D2: one message on B1I and B3I. B3I needs no term at all.
        GNSSDecoder.BeiDouDNAVData => [BeiDouB1I(), BeiDouB3I()],
        # B-CNAV1/2: own band's data component and pilot. Each carries the other
        # band's pilot delay as well, which is cross-band.
        GNSSDecoder.BeiDouB1CData => [BeiDouB1C_D(), BeiDouB1C_P()],
        GNSSDecoder.BeiDouB2aData => [BeiDouB2aI(), BeiDouB2aQ()],
        # B-CNAV3: B2b_I has no pilot and no sibling ISC.
        GNSSDecoder.BeiDouB2bData => [BeiDouB2bI()],
    ]
    # Every civil ranging signal this package can be handed, by constellation. A new
    # signal added to GNSSSignals must be listed here and placed in `expected`.
    all_signals = Dict(
        :GPS => [GPSL1CA(), GPSL1C_D(), GPSL1C_P(), GPSL2CM(), GPSL2CL(), GPSL5I(), GPSL5Q()],
        :Galileo => [GalileoE1B(), GalileoE1C(), GalileoE1B_BOC11(), GalileoE1C_BOC11(),
            GalileoE5aI(), GalileoE5aQ(), GalileoE5aQP(),
            GalileoE5bI(), GalileoE5bQ(), GalileoE6B(), GalileoE6C()],
        :BeiDou => [BeiDouB1I(), BeiDouB3I(), BeiDouB1C_D(), BeiDouB1C_P(),
            BeiDouB2aI(), BeiDouB2aQ(), BeiDouB2bI()],
    )
    constellation(D) =
        D <: GNSSDecoder.AbstractGPSData ? :GPS :
        D <: GNSSDecoder.AbstractGalileoData ? :Galileo : :BeiDou

    for (D, servable) in expected
        served = Set(GNSSSignals.get_signal_id.(servable))
        for signal in all_signals[constellation(D)]
            id = GNSSSignals.get_signal_id(signal)
            has = hasmethod(
                PositionVelocityTime.correct_by_group_delay,
                Tuple{GNSSDecoder.GNSSDecoderState{<:D},typeof(signal),Float64},
            )
            # Compared as triples so a failure names the offending cell.
            @test (nameof(D), id, has) == (nameof(D), id, id in served)
        end
    end
end
