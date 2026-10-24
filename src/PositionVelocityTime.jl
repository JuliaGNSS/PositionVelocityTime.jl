module PositionVelocityTime
using CoordinateTransformations,
    Geodesy,
    GNSSDecoder,
    GNSSSignals,
    LinearAlgebra,
    StaticArrays,
    Unitful,
    Dates

using Unitful: s, Hz, m, °, ustrip
using Dictionaries: Dictionary, IndexError, set!

include("tai_time.jl")
using .TAITimes: TAITime
# Unexported but documented decoder vocabulary this package shares: the week
# length, and the two message-family supertypes the propagator dispatches on.
using GNSSDecoder: SECONDS_PER_WEEK, AbstractGPSCNAVData, AbstractBeiDouCNAVData

export calc_pvt,
    calc_pvt!,
    PVTWorkspace,
    PVTSolution,
    TAITime,
    SatelliteInfo,
    InterFrequencyBias,
    SatelliteState,
    get_lla,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    calc_elevation_and_azimuth

# The pieces `calc_pvt` is built from — for a consumer running its own estimator over the
# same measurement model, such as a navigation filter — stay unexported but documented
# ("The Measurement-Model Surface" in the API reference); bind them with
# `using PositionVelocityTime: calc_transmit_time, …`.

# Naming conventions. Names a user meets — exports, documented functions and their
# arguments, struct fields — are spelled out ("satellite", not "sat"); internals (locals,
# private helpers' arguments, `PVTWorkspace`'s fields) may keep short forms. `get_` reads
# or converts what its argument holds (`get_week`, `get_lla`), `calc_` computes, and
# physical models and named quantities keep noun names (`ionospheric_delay`,
# `day_of_year`). `ξ` and `H` are the textbook least-squares state and design matrix.

"""
    SPEED_OF_LIGHT

The speed of light (m/s), of every range/time conversion here.
"""
const SPEED_OF_LIGHT = 299792458.0

"""
    EARTH_ROTATION_RATE

The Earth's rotation rate ω_e (rad/s), the WGS-84 value GPS and Galileo (GTRF) both use,
which turns the ECEF frame over a signal's travel time.
"""
const EARTH_ROTATION_RATE = 7.2921151467e-5

# In dependency order: a file annotates only with types the files before it define.
include("collection.jl")
include("time_systems.jl")
include("solution.jl")
# The trim-safe least-squares fit behind `estimate_position`.
include("levenberg_marquardt.jl")
using .LevenbergMarquardt: curve_fit!, LMWorkspace, grown
include("bias_layout.jl")
include("measurement_model.jl")
include("estimation.jl")
include("satellite_time.jl")
include("satellite_position.jl")
include("atmosphere/ionosphere.jl")
include("atmosphere/klobuchar.jl")
include("atmosphere/ntcm_g.jl")
include("atmosphere/bdgim.jl")
include("atmosphere/troposphere.jl")
include("atmosphere/delays.jl")
include("solver.jl")
# The precompile workload solves real fixtures, so it must come after every
# method it dispatches into.
include("precompile/states.jl")
include("precompile/workload.jl")
end
