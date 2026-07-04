# Real-world satellite states used as inputs to `calc_pvt` benchmarks so that
# timings reflect realistic decoder data and a converging least-squares geometry.
# The fixture functions live in test/fixtures.jl, shared with the test suite, so
# the benchmark inputs cannot drift from the test fixtures (that file also carries
# the GNSSDecoder major-version shim AirspeedVelocity needs to run this script
# against the base revision).

using GNSSDecoder
using GNSSSignals
using PositionVelocityTime
using Unitful: Hz

include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

"5 Galileo E1B satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_galileo_states() = galileo_e1b_states(0.0Hz)

"9 GPS L1 C/A satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_gps_states() = gps_l1_states(0.0Hz)
