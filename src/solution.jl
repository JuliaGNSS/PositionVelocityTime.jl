# What `calc_pvt` returns, and what is derived from or done to a solution.

"""
    DOP

Dilution of precision of the geometry a PVT solution was solved from.

# Fields
- `GDOP::Float64`: geometric DOP (all parameters)
- `PDOP::Float64`: position DOP
- `VDOP::Float64`: vertical DOP
- `HDOP::Float64`: horizontal DOP
- `TDOP::Float64`: time DOP
"""
struct DOP
    GDOP::Float64
    PDOP::Float64
    VDOP::Float64
    HDOP::Float64
    TDOP::Float64
end

"""
    SatelliteInfo

One satellite's entry in a [`PVTSolution`](@ref).

# Fields
- `position::ECEF{Float64}`: ECEF position at transmit time (m).
- `transmit_time::typeof(1.0s)`: transmit time, as a time of week on the satellite's own
  system.
- `residual::typeof(1.0m)`: post-fit pseudorange residual of the (atmosphere-corrected)
  pseudorange — a fit-quality and outlier indicator.
- `rate_residual::typeof(1.0m/s)`: post-fit range-rate residual of the carrier-Doppler
  velocity solve; it flags a satellite whose Doppler disagrees with the velocity fix
  (cycle slips, dynamics) independently of its pseudorange.

Both are *measured minus modelled* ("observed minus computed"), as RTKLIB's `rescode` and
`resdop` are, and GNSS-SDR's and PocketSDR's through it. The range rate is taken in the
geometric sense, positive while the satellite recedes, as `resdop`'s is; the
Doppler-signed rate a tracking loop works in, `λ · carrier_doppler`, is its negative.
This is the convention every residual and rate in the package follows.
"""
struct SatelliteInfo
    position::ECEF{Float64}
    transmit_time::typeof(1.0s)
    residual::typeof(1.0m)
    rate_residual::typeof(1.0m/s)
end

"""
    InterFrequencyBias

A receiver inter-frequency bias in a [`PVTSolution`](@ref): the differential delay of one
band's RF chain.

# Fields
- `value::typeof(1.0m)`: how much longer this band's receiver delay is than
  `reference_band`'s (m).
- `reference_band::Symbol`: the band (e.g. `:L1`; see `get_band_id`) whose delay is folded
  into the clock biases. Chosen per coverage component (see
  [`inter_frequency_bias_layout`](@ref)), so a disconnected solution can carry several.
"""
struct InterFrequencyBias
    value::typeof(1.0m)
    reference_band::Symbol
end

"""
    PVTSolution

Position, velocity and time solved from GNSS measurements.

# Fields
- `position::ECEF{Float64}`: user position (m).
- `velocity::ECEF{Float64}`: user velocity (m/s).
- `course_over_ground::typeof(1.0°)`: direction of travel, the azimuth of the horizontal
  velocity in the local East-North-Up frame, clockwise from true North in `[0, 360)°`.
  It is the direction of motion, not the vehicle heading a single antenna cannot observe;
  `0°` when the horizontal velocity is zero.
- `clock_bias::typeof(1.0m)`: receiver clock bias of `reference_system` (m); another
  system's is `clock_bias + inter_system_biases[system]`.
- `time::Union{TAITime, Nothing}`: the epoch, on TAI (see [`TAITime`](@ref)).
- `clock_drift::typeof(1.0m/s)`: receiver clock drift, the rate of `clock_bias` (m/s);
  divided by the speed of light it is the relative drift (s/s).
- `dop::Union{DOP, Nothing}`: dilution of precision.
- `satellites::Dictionary{Tuple{Symbol, Int}, SatelliteInfo}`: per satellite, keyed by
  `(signal, PRN)` with the signal id of `get_signal_id` (e.g. `:GPSL1CA`), which keeps a
  PRN apart across constellations and across the signals of one satellite (see
  [`SatelliteInfo`](@ref)).
- `reference_system::Union{SupportedTimeSystem, Nothing}`: the time system (e.g. `GPST()`)
  that `time` and `clock_bias` refer to.
- `inter_system_biases::Dictionary{SupportedTimeSystem, typeof(1.0m)}`: for each other time
  system, its clock bias minus the reference system's (m). It is a system effect, not a
  receiver one — the offset between the two system times, so the Galileo entry is
  `−c · Δt` with `Δt = GST − GPST` (the GGTO), and the BeiDou one carries BDT's steering
  residual, never the defined 14 s (see [`leap_second_offset_to_gpst`](@ref)). It is
  estimated from the geometry wherever observable; a broadcast GGTO or BGTO fixes it only
  where the satellites cannot (see [`decide_bias_layout`](@ref)). The difference of two
  entries is the bias between those two systems; empty for a single-system solution.
- `inter_frequency_biases::Dictionary{Symbol, InterFrequencyBias}`: for each band beyond a
  reference band (keyed by `get_band_id`, e.g. `:L5`), the receiver inter-frequency bias,
  shared by every constellation on that band. Empty for a single-band solution.

Immutable. [`calc_pvt!`](@ref) returns a new solution that takes over the containers of
the one handed to it; [`calc_pvt`](@ref) builds new ones and never modifies `prev_pvt`.
"""
@kwdef struct PVTSolution
    position::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    velocity::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    course_over_ground::typeof(1.0°) = 0.0°
    clock_bias::typeof(1.0m) = 0.0m
    time::Union{TAITime,Nothing} = nothing
    clock_drift::typeof(1.0m/s) = 0.0m/s
    dop::Union{DOP,Nothing} = nothing
    satellites::Dictionary{Tuple{Symbol,Int},SatelliteInfo} =
        Dictionary{Tuple{Symbol,Int},SatelliteInfo}()
    reference_system::Union{SupportedTimeSystem,Nothing} = nothing
    inter_system_biases::Dictionary{SupportedTimeSystem,typeof(1.0m)} =
        Dictionary{SupportedTimeSystem,typeof(1.0m)}()
    inter_frequency_biases::Dictionary{Symbol,InterFrequencyBias} =
        Dictionary{Symbol,InterFrequencyBias}()
end

"""
    calc_course_over_ground(position::ECEF, velocity::ECEF) -> typeof(1.0°)

Course over ground: see [`PVTSolution`](@ref)'s `course_over_ground`. `ENUfromECEF` maps
`position` to its origin, so applying it to `position + velocity` rotates the velocity
into the ENU frame without translating it.
"""
function calc_course_over_ground(position::ECEF, velocity::ECEF)
    enu_velocity = ENUfromECEF(position, wgs84)(ECEF(position + velocity))
    rad2deg(mod2pi(atan(enu_velocity[1], enu_velocity[2]))) * °
end

"""
    copy(pvt::PVTSolution) -> PVTSolution

A copy of `pvt` that shares no container with it, so a later [`calc_pvt!`](@ref) into
`pvt`'s containers leaves it untouched: how to keep a solution a `calc_pvt!` loop reuses.
"""
Base.copy(pvt::PVTSolution) = PVTSolution(
    pvt.position,
    pvt.velocity,
    pvt.course_over_ground,
    pvt.clock_bias,
    pvt.time,
    pvt.clock_drift,
    pvt.dop,
    copy(pvt.satellites),
    pvt.reference_system,
    copy(pvt.inter_system_biases),
    copy(pvt.inter_frequency_biases),
)

# A copy of `source` in the containers of `solution`, or `source` itself where they share
# them: copying a container onto itself would empty it.
function copy_into(solution::PVTSolution, source::PVTSolution)
    solution.satellites === source.satellites &&
        solution.inter_system_biases === source.inter_system_biases &&
        solution.inter_frequency_biases === source.inter_frequency_biases &&
        return source
    copy_container!(solution.satellites, source.satellites)
    copy_container!(solution.inter_system_biases, source.inter_system_biases)
    copy_container!(solution.inter_frequency_biases, source.inter_frequency_biases)
    PVTSolution(
        source.position,
        source.velocity,
        source.course_over_ground,
        source.clock_bias,
        source.time,
        source.clock_drift,
        source.dop,
        solution.satellites,
        source.reference_system,
        solution.inter_system_biases,
        solution.inter_frequency_biases,
    )
end

function copy_container!(destination::Dictionary, source::Dictionary)
    destination === source && return destination
    empty_keeping_capacity!(destination)
    for (key, value) in pairs(source)
        set!(destination, key, value)
    end
    destination
end

"""
    empty_keeping_capacity!(dict::Dictionary) -> dict
    empty_keeping_capacity!(solution::PVTSolution) -> solution

Empty a `Dictionary` keeping its storage, so a refill of as many entries allocates
nothing (`empty!` hands it a fresh, minimal hash table). The [`PVTSolution`](@ref) form
empties all three of a solution's containers, for a consumer that builds its solutions
into the containers of one it keeps, as [`calc_pvt!`](@ref) does.
"""
function empty_keeping_capacity! end

# The empty state Dictionaries' own `empty!` builds, with the slot table kept at its
# current (power-of-two) size. It reaches into `Dictionaries.Indices` (0.4), so a test in
# `test/allocation_free.jl` pins it.
function empty_keeping_capacity!(dict::Dictionary)
    indices = keys(dict)
    fill!(getfield(indices, :slots), 0)
    empty!(getfield(indices, :hashes))
    empty!(getfield(indices, :values))
    setfield!(indices, :holes, 0)
    empty!(getfield(dict, :values))
    dict
end

function empty_keeping_capacity!(solution::PVTSolution)
    empty_keeping_capacity!(solution.satellites)
    empty_keeping_capacity!(solution.inter_system_biases)
    empty_keeping_capacity!(solution.inter_frequency_biases)
    solution
end

"""
    get_lla(pvt::PVTSolution) -> LLA

The solution's position as WGS84 geodetic latitude, longitude and altitude.
"""
function get_lla(pvt::PVTSolution)
    LLAfromECEF(wgs84)(pvt.position)
end
