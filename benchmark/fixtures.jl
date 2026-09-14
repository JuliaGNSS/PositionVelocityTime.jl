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
using Unitful: Hz

include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

# AirspeedVelocity runs *this* script against both the PR and the base revision, and
# until the base revision is 6.0 its `calc_pvt` still takes a flat vector of
# `SatelliteState`s where 6.0 takes signal groups. So every input below is built through
# the shape-neutral helpers here: on 6.0 they produce groups, on 5.x they hand the
# vectors straight back, and the two revisions stay comparable instead of the base run
# failing to load. Drop them once the base revision of the comparison is 6.0 — the same
# treatment, and the same deletion condition, as the GNSSDecoder major-version shim in
# `test/fixtures.jl`.
#
# `getfield` rather than a plain `PositionVelocityTime.SignalGroup`: on the base
# revision that binding does not exist, and naming it — even in a branch that never
# runs — makes the module's import of it undeclared at load time.
const HAS_SIGNAL_GROUPS = isdefined(PositionVelocityTime, :SignalGroup)

signal_group(signal, satellites) =
    getfield(PositionVelocityTime, :SignalGroup)(signal, satellites)

"""
    pvt_input(signal, states)

The satellites of one ranging signal, in whatever shape this revision's `calc_pvt`
takes. On 6.0 that is a `SignalGroup` whose satellites are a `Dictionary` keyed by PRN —
the shape a receiver carries, and what the Tracking extension builds.
"""
pvt_input(signal, states) =
    HAS_SIGNAL_GROUPS ?
    signal_group(signal, Dictionary([state.decoder.prn for state in states], states)) :
    states

"""
    combine_inputs(inputs::NamedTuple)

One mixed-constellation epoch from several per-signal inputs: a NamedTuple of groups on
6.0, the pooled flat vector — whose element type is abstract, which is the whole point —
on 5.x.
"""
combine_inputs(inputs::NamedTuple) =
    HAS_SIGNAL_GROUPS ? inputs : reduce(vcat, values(inputs))

"Number of satellites in one per-signal input."
num_satellites(input) = HAS_SIGNAL_GROUPS ? length(input.satellites) : length(input)

"The first `n` satellites of one per-signal input, in the same shape."
function take_satellites(input, n)
    HAS_SIGNAL_GROUPS || return input[1:n]
    satellites = input.satellites
    signal_group(
        input.signal,
        Dictionary(collect(keys(satellites))[1:n], collect(satellites)[1:n]),
    )
end

"5 Galileo E1B satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_galileo_input() = pvt_input(GalileoE1B(), galileo_e1b_states(0.0Hz))

"9 GPS L1 C/A satellites over Aachen, 2021-05-31 (from test/fixtures.jl)."
make_gps_input() = pvt_input(GPSL1CA(), gps_l1_states(0.0Hz))
