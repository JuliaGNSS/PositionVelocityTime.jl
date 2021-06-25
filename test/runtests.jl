
using Test, PositionVelocityTime, GNSSDecoder
import PositionVelocityTime: GNSSDecoderState, GPSData, GPSL1Constants, ECEF, SatelliteState


include("test_data.jl")
include("sat_position.jl")
include("user_position.jl")
include("position_velocity_time.jl")
include("sv_time.jl")