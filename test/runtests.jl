
using Test, PositionVelocityTime, GNSSDecoder, AstroTime, GNSSSignals, Geodesy, Dates, LinearAlgebra
using Unitful: Hz, m, s, °, ustrip
# `Dictionaries.Dictionary` is how GNSSDecoder keys its paged records (BeiDou B1C's
# per-GNSS BGTO sets, the almanacs); the BeiDou tests build one directly.
using Dictionaries: Dictionary

include("aqua.jl")
include("fixtures.jl")
include("sat_time.jl")
include("sat_position.jl")
include("pvt.jl")
include("dop.jl")
include("cnav.jl")
include("gps_l2c.jl")
include("galileo_e5a.jl")
include("galileo_e5b_e6b.jl")
include("beidou.jl")
include("inter_frequency_bias.jl")
include("pvt_iono_tropo.jl")
include("get_week.jl")
include("tracking_ext.jl")
include("ionosphere.jl")
include("troposphere.jl")
include("pvt_integration.jl")
