module PositionVelocityTime
using CoordinateTransformations,
    DocStringExtensions,
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
    SatInfo,
    InterFrequencyBias,
    SatelliteState,
    get_LLA,
    get_sat_info,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    get_sat_enu

# The measurement-model surface. `calc_pvt` above is the whole scalar solver,
# but a consumer running its own estimator over the same measurement model — a
# navigation filter closing tracking loops through its own Kalman update —
# consumes the model in pieces: the per-satellite times and clock rates, the
# bias layout and design-matrix columns, the predicted ranges and geometry, the
# atmospheric and time-system corrections. Those pieces stay unexported — they
# are solver internals, not names every `using` should carry — but they are
# documented as a surface (see "The Measurement-Model Surface" in the API
# reference), and a consumer binds them explicitly with
# `using PositionVelocityTime: calc_corrected_time, …`, which declares the
# dependency at a single site.

"""
    SPEED_OF_LIGHT

The speed of light (m/s), shared by every range/time conversion in this
package and its consumers.
"""
const SPEED_OF_LIGHT = 299792458.0

"""
    EARTH_ROTATION_RATE

The Earth's rotation rate ω_e (rad/s), the WGS-84 value GPS and Galileo (GTRF) both
use. It sets the Earth's turn over a signal's travel time, which carries the satellite
from the ECEF frame of transmission into the one of reception.
"""
const EARTH_ROTATION_RATE = 7.2921151467e-5

# In dependency order: a file annotates only with types the files before it define.
include("collection.jl")
include("time_systems.jl")
include("solution.jl")
# The trim-safe least-squares fit behind `user_position`.
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
