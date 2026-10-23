# Real-world satellite states used as inputs to `calc_pvt` benchmarks so that
# timings reflect realistic decoder data and a converging least-squares geometry.
# The fixture functions live in test/fixtures.jl, shared with the test suite, so
# the benchmark inputs cannot drift from the test fixtures.

using Dictionaries
using GNSSDecoder
using GNSSSignals
using PositionVelocityTime
using Unitful: Hz

include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

"""
    pvt_input(signal, states)

The satellites of one ranging signal as the `SignalGroup` `calc_pvt` takes, keyed by
PRN — the shape a receiver carries, and what Tracking.jl builds.
"""
pvt_input(signal, states) = PositionVelocityTime.SignalGroup(
    signal, Dictionary([state.decoder.prn for state in states], states))

"The first `n` satellites of one per-signal input, in the same shape."
take_satellites(input, n) = PositionVelocityTime.SignalGroup(
    input.signal,
    Dictionary(collect(keys(input.satellites))[1:n], collect(input.satellites)[1:n]),
)

"5 Galileo E1B satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_galileo_input() = pvt_input(GalileoE1B(), galileo_e1b_states(0.0Hz))

"9 GPS L1 C/A satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_gps_input() = pvt_input(GPSL1CA(), gps_l1_states(0.0Hz))
