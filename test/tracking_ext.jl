using Tracking
using Tracking: TrackedSat, get_code_phase, get_carrier_doppler, get_carrier_phase

@testset "Tracking extension: SatelliteState from a Tracking sat state" begin
    # The extension is only available once Tracking is loaded.
    @test Base.get_extension(PositionVelocityTime, :PositionVelocityTimeTrackingExt) !==
          nothing

    gpsl1 = GPSL1CA()
    decoder = GNSSDecoderState(gpsl1, 1)
    sat_state = TrackedSat(gpsl1, 1, 123.0, 1234.0Hz; carrier_phase = 0.25)

    state = SatelliteState(decoder, gpsl1, sat_state)

    @test state isa SatelliteState
    @test state.decoder === decoder
    @test state.system === gpsl1
    # Each measurement must be wired from the matching Tracking getter.
    @test state.code_phase == get_code_phase(sat_state)
    @test state.carrier_doppler == get_carrier_doppler(sat_state)
    @test state.carrier_phase == get_carrier_phase(sat_state)
end

@testset "Tracking extension: signal groups from a TrackState" begin
    # `SignalGroup` below is Tracking's — this file does `using Tracking`, and
    # PositionVelocityTime deliberately does not export its own type of that name, so
    # the two coexist here with no clash and PVT's is written out in full.
    @test SignalGroup === Tracking.SignalGroup

    track_state = TrackState(; signals = (gps = (GPSL1CA(),), galileo = (GalileoE1B(),)))
    for (group, prn, code_phase) in
        ((:gps, 7, 12.0), (:gps, 8, 13.0), (:galileo, 2, 5.0))
        track_state = Tracking.add_satellite!(
            track_state; group, prn, carrier_doppler = 100.0Hz, code_phase)
    end
    # PRN 8 is tracked but has no decoder yet — the ordinary state of a satellite for
    # the first minute of its life — so it is skipped rather than erroring.
    decoders = (
        gps = Dictionary([7], [GNSSDecoderState(GPSL1CA(), 7)]),
        galileo = Dictionary([2], [GNSSDecoderState(GalileoE1B(), 2)]),
    )

    groups = PositionVelocityTime.signal_groups(track_state, decoders)
    # Both levels of key survive: the group names and, within each group, the PRNs.
    @test keys(groups) == keys(track_state.groups)
    @test collect(keys(groups.gps.satellites)) == [7]
    @test collect(keys(groups.galileo.satellites)) == [2]
    # The ranging signal defaults to the group's first signal, and it is what the
    # satellite states are built on.
    @test groups.gps.signal === first(track_state.groups.gps.signals)
    @test groups.gps.satellites[7].system === groups.gps.signal
    @test groups.gps.satellites[7].decoder === decoders.gps[7]
    # The measurements are wired from the matching Tracking getters, as in the
    # single-satellite constructor above.
    tracked = Tracking.get_sat_state(track_state, :gps, 7)
    @test groups.gps.satellites[7].code_phase == get_code_phase(tracked)
    @test groups.gps.satellites[7].carrier_doppler == get_carrier_doppler(tracked)

    # Each group's satellites are concretely typed, which is the whole point, and the
    # conversion is fully inferable — unlike `signal_groups` over a flat vector, whose
    # NamedTuple type is only known at runtime.
    @test isconcretetype(eltype(groups.gps.satellites))
    @test inferred_return_type(
        PositionVelocityTime.signal_groups, Tuple{typeof(track_state), typeof(decoders)}
    ) === typeof(groups)

    # The single-group constructor, and its `signal` keyword — the escape hatch for a
    # range generated on a signal of the group other than the driver.
    one_group = PositionVelocityTime.SignalGroup(track_state.groups.gps, decoders.gps)
    @test collect(keys(one_group.satellites)) == [7]
    @test one_group.signal === groups.gps.signal
    ranging_on = GPSL1CA()
    @test PositionVelocityTime.SignalGroup(
        track_state.groups.gps, decoders.gps; signal = ranging_on,
    ).signal === ranging_on
    # A group whose satellites have no decoders at all converts to an empty group
    # rather than failing.
    @test isempty(
        PositionVelocityTime.SignalGroup(
            track_state.groups.gps, Dictionary{Int,GNSSDecoderState}(),
        ).satellites,
    )
end
