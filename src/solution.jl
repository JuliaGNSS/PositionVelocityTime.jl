# What `calc_pvt` returns, and what is derived from or done to a solution.

"""
    DOP

Dilution of precision of the geometry a PVT solution was solved from, as
[`PVTSolution`](@ref)'s `dop` (e.g. `pvt.dop.GDOP`). `NaN` throughout for a geometry that
fails the rank test.

# Fields
- `GDOP::Float64`: geometric DOP (all parameters)
- `PDOP::Float64`: position DOP
- `HDOP::Float64`: horizontal DOP
- `VDOP::Float64`: vertical DOP
- `TDOP::Float64`: time DOP
"""
struct DOP
    GDOP::Float64
    PDOP::Float64
    HDOP::Float64
    VDOP::Float64
    TDOP::Float64
end

"""
    ObservationResult

One signal observation's result in a [`PVTSolution`](@ref), keyed by signal and PRN.

# Fields
- `position::ECEF{Float64}`: ECEF position at transmit time (m).
- `transmit_time::typeof(1.0s)`: transmit time, as a time of week on the satellite's own
  system.
- `residual::typeof(1.0m)`: post-fit pseudorange residual of the (atmosphere-corrected)
  pseudorange — a fit-quality and outlier indicator.
- `rate_residual::typeof(1.0m/s)`: post-fit range-rate residual of the carrier-Doppler
  velocity solve; it flags a satellite whose Doppler disagrees with the velocity fix
  (cycle slips, dynamics) independently of its pseudorange.

`residual` and `rate_residual` are *measured minus modelled* ("observed minus
computed"), as RTKLIB's `rescode` and `resdop` are, and so GNSS-SDR's and PocketSDR's. The range rate is taken in the
geometric sense, positive while the satellite recedes, as `resdop`'s is; the
Doppler-signed rate a tracking loop works in, `λ · carrier_doppler`, is its negative.
This is the convention every residual and rate in the package follows.
"""
struct ObservationResult
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
  [`ifb_layout`](@ref)), so a disconnected solution can carry several.
"""
struct InterFrequencyBias
    value::typeof(1.0m)
    reference_band::Symbol
end

"""
    BroadcastClockConstraint

A receiver-clock difference a [`PVTSolution`](@ref) constrained with a broadcast GNSS
time offset (GGTO, BGTO), where the observations could not separate the two receiver
clocks (see [`decide_bias_layout`](@ref)). This assumes any differential receiver delay
between those clocks is captured by the shared per-band delay model. Unmodelled receiver
delays and errors in the broadcast offset therefore enter the constrained result.
Each constraint fixes an unknown: a fix estimates `3 + (1 + length(inter_system_biases)) +
length(inter_frequency_biases) − length(broadcast_clock_constraints)` parameters.

# Fields
- `system::SupportedTimeSystem`: the time system whose broadcast offset was used.
- `target::SupportedTimeSystem`: the time system the offset is toward.
- `value::typeof(1.0m)`: the bias of `system` against `target`, as
  [`broadcast_inter_system_bias`](@ref) evaluates it. The solution meets it: the
  inter-system bias of `system` minus that of `target`, the reference system's counting as
  zero, equals it to rounding.
"""
struct BroadcastClockConstraint
    system::SupportedTimeSystem
    target::SupportedTimeSystem
    value::typeof(1.0m)
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
- `time::Union{TAITime, Nothing}`: the epoch, on TAI (GNSSSignals' `TAITime`).
- `clock_bias::typeof(1.0m)`: receiver clock bias of `reference_system` (m); another
  system's is `clock_bias + inter_system_biases[system]`.
- `clock_drift::typeof(1.0m/s)`: receiver clock drift, the rate of `clock_bias` (m/s);
  divided by the speed of light it is the relative drift (s/s).
- `reference_system::Union{SupportedTimeSystem, Nothing}`: the time system (e.g. `GPST()`)
  that `time` and `clock_bias` refer to.
- `inter_system_biases::Dictionary{SupportedTimeSystem, typeof(1.0m)}`: for each other time
  system, its effective receiver clock bias minus the reference system's (m), under a
  model with one receiver delay shared by all systems on each band. It includes the
  system-time offset (`−c · Δt`, with `Δt = GST − GPST` for Galileo) and receiver delays
  not separated by that model. When constellation × band coverage is disconnected,
  different components' reference-band delays also enter this difference. The defined
  14 s between BDT and GPST is already removed by `get_time_scale_offset`. The bias is
  estimated wherever observable; a broadcast GGTO or BGTO constrains it only where the
  observations cannot (see [`decide_bias_layout`](@ref)). The difference of two entries
  is their effective clock-bias difference; empty for a single-system solution.
- `broadcast_clock_constraints::Vector{BroadcastClockConstraint}`: the
  clock differences constrained by a broadcast time offset, one per
  pair of clocks it tied (see [`BroadcastClockConstraint`](@ref)). Empty whenever the
  satellites separate every clock, the usual case.
- `inter_frequency_biases::Dictionary{Symbol, InterFrequencyBias}`: for each band beyond a
  reference band (keyed by `get_band_id`, e.g. `:L5`), the receiver inter-frequency bias,
  shared by every constellation on that band. Empty for a single-band solution.
- `dop::Union{DOP, Nothing}`: dilution of precision.
- `observations::Dictionary{Tuple{Symbol, Int}, ObservationResult}`: per observation, keyed by
  `(signal, PRN)` with the signal id of `get_signal_id` (e.g. `:GPSL1CA`), which keeps a
  PRN apart across constellations and across the signals of one satellite (see
  [`ObservationResult`](@ref)).

Immutable. [`calc_pvt!`](@ref) returns a new solution that takes over the containers of
the one handed to it; [`calc_pvt`](@ref) builds new ones and never modifies `prev_pvt`.
"""
@kwdef struct PVTSolution
    position::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    velocity::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    course_over_ground::typeof(1.0°) = 0.0°
    time::Union{TAITime,Nothing} = nothing
    clock_bias::typeof(1.0m) = 0.0m
    clock_drift::typeof(1.0m/s) = 0.0m/s
    reference_system::Union{SupportedTimeSystem,Nothing} = nothing
    inter_system_biases::Dictionary{SupportedTimeSystem,typeof(1.0m)} =
        Dictionary{SupportedTimeSystem,typeof(1.0m)}()
    broadcast_clock_constraints::Vector{BroadcastClockConstraint} =
        BroadcastClockConstraint[]
    inter_frequency_biases::Dictionary{Symbol,InterFrequencyBias} =
        Dictionary{Symbol,InterFrequencyBias}()
    dop::Union{DOP,Nothing} = nothing
    observations::Dictionary{Tuple{Symbol,Int},ObservationResult} =
        Dictionary{Tuple{Symbol,Int},ObservationResult}()
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
    pvt.time,
    pvt.clock_bias,
    pvt.clock_drift,
    pvt.reference_system,
    copy(pvt.inter_system_biases),
    copy(pvt.broadcast_clock_constraints),
    copy(pvt.inter_frequency_biases),
    pvt.dop,
    copy(pvt.observations),
)

# A copy of `source` in the containers of `solution`, or `source` itself where they share
# them: copying a container onto itself would empty it.
function copy_into(solution::PVTSolution, source::PVTSolution)
    solution.observations === source.observations &&
        solution.inter_system_biases === source.inter_system_biases &&
        solution.inter_frequency_biases === source.inter_frequency_biases &&
        solution.broadcast_clock_constraints ===
        source.broadcast_clock_constraints &&
        return source
    copy_container!(solution.observations, source.observations)
    copy_container!(solution.inter_system_biases, source.inter_system_biases)
    copy_container!(solution.inter_frequency_biases, source.inter_frequency_biases)
    copy_container!(solution.broadcast_clock_constraints,
        source.broadcast_clock_constraints)
    PVTSolution(
        source.position,
        source.velocity,
        source.course_over_ground,
        source.time,
        source.clock_bias,
        source.clock_drift,
        source.reference_system,
        solution.inter_system_biases,
        solution.broadcast_clock_constraints,
        solution.inter_frequency_biases,
        source.dop,
        solution.observations,
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
copy_container!(destination::Vector, source::Vector) =
    destination === source ? destination : copy!(destination, source)

"""
    empty_keeping_capacity!(dict::Dictionary) -> dict
    empty_keeping_capacity!(solution::PVTSolution) -> solution

Empty a `Dictionary` keeping its storage, so a refill of as many entries allocates
nothing (`empty!` hands it a fresh, minimal hash table). The [`PVTSolution`](@ref) form
empties all four of a solution's containers (the `Vector` by `empty!`, which keeps its
capacity), for a consumer that builds its solutions
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
    empty_keeping_capacity!(solution.observations)
    empty_keeping_capacity!(solution.inter_system_biases)
    empty_keeping_capacity!(solution.inter_frequency_biases)
    empty!(solution.broadcast_clock_constraints)
    solution
end

"""
    get_lla(pvt::PVTSolution) -> LLA

The solution's position as WGS84 geodetic latitude, longitude and altitude.
"""
function get_lla(pvt::PVTSolution)
    LLAfromECEF(wgs84)(pvt.position)
end
