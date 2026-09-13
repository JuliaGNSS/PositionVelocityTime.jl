# Real-world satellite states used as inputs to `calc_pvt` benchmarks so that
# timings reflect realistic decoder data and a converging least-squares geometry.
# The fixture functions live in test/fixtures.jl, shared with the test suite, so
# the benchmark inputs cannot drift from the test fixtures (that file also carries
# the GNSSDecoder major-version shim AirspeedVelocity needs to run this script
# against the base revision).

using Dictionaries
using GNSSDecoder
using GNSSSignals
using PositionVelocityTime
using PositionVelocityTime: SignalGroup
using Unitful: Hz

include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

# `calc_pvt` takes signal groups. The satellites are stored in a `Dictionary` keyed by
# PRN, which is the shape a receiver carries (and what the Tracking extension builds),
# rather than the plain vector the fixture functions return.
signal_group(signal, states) =
    SignalGroup(signal, Dictionary([state.decoder.prn for state in states], states))

"5 Galileo E1B satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_galileo_group() = signal_group(GalileoE1B(), galileo_e1b_states(0.0Hz))

"9 GPS L1 C/A satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_gps_group() = signal_group(GPSL1CA(), gps_l1_states(0.0Hz))

"The first `n` satellites of `group`, as a group on the same ranging signal."
take_satellites(group, n) =
    SignalGroup(group.signal, Dictionary(collect(keys(group.satellites))[1:n],
        collect(group.satellites)[1:n]))
