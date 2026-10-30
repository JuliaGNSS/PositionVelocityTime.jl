module PositionVelocityTime
using CoordinateTransformations,
    Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, StaticArrays, Unitful, Dates

using Unitful: s, Hz, m, °, ustrip
using Dictionaries: Dictionary, IndexError, set!

# GNSSSignals' week arithmetic, documented but unexported there.
using GNSSSignals: fold_week_crossover

export calc_pvt,
    calc_pvt!,
    PVTWorkspace,
    PVTSolution,
    PVTResult,
    PVTStatus,
    ObservationResult,
    InterFrequencyBias,
    BroadcastClockConstraint,
    DOP,
    SupportedTimeSystem,
    SignalObservation,
    get_lla,
    calc_satellite_position_and_velocity,
    calc_elevation_and_azimuth

# The pieces `calc_pvt` is built from — for a consumer running its own estimator over the
# same measurement model, such as a navigation filter — stay unexported but documented
# ("The Measurement-Model Surface" in the API reference); bind them with
# `using PositionVelocityTime: calc_transmit_time, …`.

# Naming conventions. Names a user writes — exports, documented functions, keyword
# arguments, struct fields — are spelled out ("satellite", not "sat"), except the bias
# layout's IFB names (`ifb_layout`, `ifb_columns`, `num_ifbs`); positional arguments and
# internals (locals, private helpers, `PVTWorkspace`'s fields) may use short forms (`sat`,
# `pos`, `vel`, `ifb`, `mjd`). `get_` reads
# or converts what its argument holds (`get_lla`), `calc_` computes, and
# physical models and named quantities keep noun names (`ionospheric_delay`,
# `tropospheric_delay`). `ξ` and `H` are the textbook least-squares state and design
# matrix.

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
