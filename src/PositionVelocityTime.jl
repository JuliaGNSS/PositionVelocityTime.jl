module PositionVelocityTime
using CoordinateTransformations,
    Geodesy, GNSSDecoder, GNSSSignals, LinearAlgebra, StaticArrays, Unitful, Dates

using Unitful: s, Hz, m, °, ustrip
using Dictionaries: Dictionary, set!

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

# Julia 1.11+ records the supported qualified API without importing it into every
# caller's namespace. Build the expression at load time to retain Julia 1.10 support.
@static if VERSION >= v"1.11"
    eval(
        Expr(
            :public,
            :PreparedObservation,
            :prepare_observations,
            :prepare_observations!,
            :SUPPORTED_TIME_SYSTEMS,
            :SPEED_OF_LIGHT,
            :calc_course_over_ground,
            :calc_transmit_time,
            :calc_uncorrected_transmit_time,
            :calc_satellite_clock_drift,
            :calc_time_scale_offset,
            :calc_pseudorange,
            :BiasColumns,
            :num_parameters,
            :ClockConstraint,
            :decide_bias_layout,
            :decide_bias_layout!,
            :BiasLayout,
            :fit_bias_columns,
            :dop_bias_columns,
            :select_reference_system,
            :predict_pseudoranges!,
            :calc_design_matrix,
            :calc_design_matrix!,
            :calc_fit_covariance!,
            :calc_dop,
            :calc_dop!,
            :calc_pseudorange_rate,
            :broadcast_time_offset,
            :has_broadcast_time_offset,
            :broadcast_inter_system_bias,
            :predict_atmospheric_delays,
            :predict_atmospheric_delays!,
            :IonosphericModel,
            :select_ionospheric_parameters,
            :ionospheric_delay,
            :tropospheric_delay,
            :ifb_layout,
            :solved,
            :insufficient_observations,
            :degenerate_geometry,
        ),
    )
end

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
