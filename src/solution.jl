# What `calc_pvt` returns: the solution, its per-satellite and per-bias entries, and what
# is derived from or done to a solution.

"""
    DOP

Dilution of Precision (DOP) values describing the geometric quality of the satellite
constellation used for a PVT solution.

# Fields
- `GDOP::Float64`: Geometric DOP (overall quality)
- `PDOP::Float64`: Position DOP (3D position quality)
- `VDOP::Float64`: Vertical DOP
- `HDOP::Float64`: Horizontal DOP
- `TDOP::Float64`: Time DOP
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

Per-satellite information attached to a [`PVTSolution`](@ref) (one entry per
satellite used in the fix).

# Fields
- `position::ECEF{Float64}`: Satellite ECEF position at transmit time (metres).
- `transmit_time::Float64`: Satellite transmit time (system time of week, seconds).
- `residual::typeof(1.0m)`: Post-fit least-squares pseudorange residual (metres) — the
  (atmosphere-corrected) measured minus the modeled pseudorange. A per-satellite
  fit-quality / outlier indicator.
- `rate_residual::typeof(1.0m/s)`: Post-fit least-squares range-rate residual (metres per
  second) — the measured minus the modeled range rate of the carrier-Doppler velocity and
  clock-drift solve. The rate-domain counterpart of `residual`: it flags a satellite whose
  Doppler disagrees with the velocity fix (cycle slips, dynamics) independently of its
  pseudorange.

Both are *measured minus modeled* ("observed minus computed"), the orientation GNSS
software reports observation residuals in — RTKLIB's `rescode` / `resdop`, and GNSS-SDR
and PocketSDR through it. Note that `rate_residual` follows `resdop` in the geometric
range-rate sense, positive while the satellite recedes; a receiver forming the same
residual from its tracking loops' `λ · carrier_doppler` works in the opposite sign (see
`estimate_velocity` and `calc_pseudorange_rate`, which set that sign).
"""
struct SatelliteInfo
    position::ECEF{Float64}
    transmit_time::Float64
    residual::typeof(1.0m)
    rate_residual::typeof(1.0m/s)
end

"""
    InterFrequencyBias

A receiver inter-frequency bias attached to a [`PVTSolution`](@ref): the differential
hardware delay of one frequency band's RF chain, together with the band it is measured
against.

# Fields
- `value::typeof(1.0m)`: The bias (metres), relative to `reference` — how much longer
  this band's receiver chain delay is than the reference band's.
- `reference_band::Symbol`: The reference band (e.g. `:L1`; see `get_band_id`) whose delay is
  folded into the clock biases and against which `value` is measured. Chosen per
  coverage component (see [`ifb_layout`](@ref)), so different biases in a
  disconnected solution can in principle carry different references.
"""
struct InterFrequencyBias
    value::typeof(1.0m)
    reference_band::Symbol
end

"""
    PVTSolution

Complete Position, Velocity, and Time solution from GNSS measurements.

# Fields
- `position::ECEF{Float64}`: User position in ECEF coordinates (meters)
- `velocity::ECEF{Float64}`: User velocity in ECEF coordinates (m/s)
- `course_over_ground::typeof(1.0°)`: Horizontal direction of travel (degrees), the
  azimuth of the velocity vector in the local East-North-Up frame at `position` —
  measured clockwise from true North and wrapped to `[0, 360)°`, following the GNSS
  course-over-ground (COG) convention. Derived from the velocity alone, so this is the
  direction of motion, not vehicle heading (which a single-antenna receiver cannot
  observe). The vertical (Up) velocity component is ignored. `0°` when the horizontal
  velocity is zero (stationary or purely vertical), where course is undefined.
- `clock_bias::typeof(1.0m)`: Estimated receiver clock bias of the reference GNSS
  (meters). For a multi-GNSS solution this is the bias of `reference_system`;
  other systems' biases are `clock_bias + inter_system_biases[system]`.
- `time::Union{TAITime, Nothing}`: Estimated time as a TAI epoch (GNSSSignals' `TAITime`;
  with AstroTime loaded, `TAIEpoch(time)` converts it)
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `satellites::Dictionary{Tuple{Symbol, Int}, SatelliteInfo}`: Maps `(signal, PRN)` to satellite
  info — position, transmit time, and the post-fit pseudorange and range-rate
  residuals (see [`SatelliteInfo`](@ref)). The signal tag (e.g. `:GPSL1CA`,
  `:GalileoE1B`; see `get_signal_id`) keeps the
  same PRN apart both across constellations (GPS PRN 5 vs Galileo E05) and across
  signals of one constellation (a satellite tracked on GPS L1 C/A and L5 yields two
  entries sharing a PRN). Receiver-clock grouping is by time system, not signal —
  see `reference_system`.
- `reference_system::Union{SupportedTimeSystem, Nothing}`: GNSS time system (e.g.
  `GNSSSignals.GPST()`, `GST()`) that `time` and `clock_bias` are referenced to.
- `inter_system_biases::Dict{SupportedTimeSystem, typeof(1.0m)}`: For each GNSS time system
  other than `reference_system`, the offset of that system's time scale relative to the reference
  system's (meters) — the inter-system bias. This is a **system / space-segment** effect,
  not a receiver one (this receiver has no inter-system hardware bias): it is the GNSS
  system-time offset, so the Galileo entry equals `−c · Δt_systems`, `Δt_systems = GST −
  GPST` (the GGTO), and the BeiDou entry likewise carries the BDT steering residual (the
  BGTO's subject) — never the defined 14 s count offset, which is a convention already
  removed from the measurements (see [`calc_time_scale_offsets`](@ref)). It is
  **estimated directly from the geometry whenever observable** (no broadcast error); the
  broadcast GGTO/BGTO is used to derive it only as a fallback (the offset-aided
  collapse), as that broadcast value may be erroneous. Reference-independent
  (the difference of two entries is the offset between those two systems); empty for a
  single-system solution.
- `inter_frequency_biases::Dict{Symbol, InterFrequencyBias}`: For each frequency band
  other than the reference band, the receiver inter-frequency bias relative to it — the
  differential hardware delay of that band's RF chain, estimated as an extra unknown
  when satellites are processed on more than one band. The key is the band
  (e.g. `:L1`, `:L5`; see `get_band_id`) and the value is shared across all
  constellations on that band. Each [`InterFrequencyBias`](@ref) carries both the bias
  (in metres) and the reference band it is measured against — the anchor whose delay is
  folded into the clock biases — so each reported bias's anchor is explicit rather than
  implicit. The reference is chosen per coverage component (see
  [`ifb_layout`](@ref)), so a disconnected constellation can in principle yield
  several references. Empty for a single-band solution. A solution then needs
  `n ≥ 3 + M + B` satellites for `M` time systems and `B` extra bands.

Immutable: a solution's fields cannot be reassigned once it is built. The allocation-free
[`calc_pvt!`](@ref) returns a new solution too; what it reuses are the `satellites`,
`inter_system_biases` and `inter_frequency_biases` containers of the solution it is
explicitly handed as its output, which the returned solution takes over.
[`calc_pvt`](@ref) builds new containers every time and never modifies the `prev_pvt`
it reads.
"""
@kwdef struct PVTSolution
    position::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    velocity::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    course_over_ground::typeof(1.0°) = 0.0°
    clock_bias::typeof(1.0m) = 0.0m
    time::Union{TAITime,Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    satellites::Dictionary{Tuple{Symbol,Int},SatelliteInfo} =
        Dictionary{Tuple{Symbol,Int},SatelliteInfo}()
    reference_system::Union{SupportedTimeSystem,Nothing} = nothing
    inter_system_biases::Dict{SupportedTimeSystem,typeof(1.0m)} =
        Dict{SupportedTimeSystem,typeof(1.0m)}()
    inter_frequency_biases::Dict{Symbol,InterFrequencyBias} =
        Dict{Symbol,InterFrequencyBias}()
end

"""
    calc_course_over_ground(position::ECEF, velocity::ECEF) -> typeof(1.0°)

Course over ground (degrees): the azimuth of `velocity` in the local East-North-Up
frame at `position`, measured clockwise from true North and wrapped to `[0, 360)°`,
following the GNSS COG convention. This is the direction of motion derived from the
velocity, not vehicle heading. The vertical (Up) component is ignored. Returns `0°`
when the horizontal velocity is zero (stationary or purely vertical), where course is
undefined.

`ENUfromECEF` is affine with `position` as its origin (which it maps to zero), so
applying it to `position + velocity` rotates the velocity vector into the ENU frame
without introducing any translation.
"""
function calc_course_over_ground(position::ECEF, velocity::ECEF)
    enu_velocity = ENUfromECEF(position, wgs84)(ECEF(position + velocity))
    rad2deg(mod2pi(atan(enu_velocity[1], enu_velocity[2]))) * °
end

"""
    copy(pvt::PVTSolution) -> PVTSolution

A copy of `pvt` that shares no container with it: `satellites`, `inter_system_biases` and
`inter_frequency_biases` are copied, every other field is an immutable value. Their
entries are immutable too, so the copy is independent of `pvt` — a later
[`calc_pvt!`](@ref) into `pvt`'s containers leaves it untouched. This is how to keep a
solution that a `calc_pvt!` loop will reuse.
"""
Base.copy(pvt::PVTSolution) = PVTSolution(
    pvt.position,
    pvt.velocity,
    pvt.course_over_ground,
    pvt.clock_bias,
    pvt.time,
    pvt.relative_clock_drift,
    pvt.dop,
    copy(pvt.satellites),
    pvt.reference_system,
    copy(pvt.inter_system_biases),
    copy(pvt.inter_frequency_biases),
)

# A copy of `source` in the containers of `solution`. `source` itself where the two
# share their containers (above all, where they are the same solution): there is
# nothing to copy, and copying a container onto itself would empty it.
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
        source.relative_clock_drift,
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

function copy_container!(destination::Dict, source::Dict)
    destination === source && return destination
    empty!(destination)
    for (key, value) in source
        destination[key] = value
    end
    destination
end

"""
    empty_keeping_capacity!(dict::Dictionary) -> dict
    empty_keeping_capacity!(solution::PVTSolution) -> solution

Empty a `Dictionary` without giving up its storage, so refilling it with as many entries
as before allocates nothing. `empty!(::Dictionary)` hands the dictionary a fresh, minimal
hash table instead, which a refill then reallocates and regrows.

The [`PVTSolution`](@ref) form empties all three of a solution's containers that way —
`satellites`, `inter_system_biases` and `inter_frequency_biases` (the last two Base `Dict`s,
whose own `empty!` already keeps their capacity) — for a consumer that builds its own
solutions into the containers of one it keeps, as [`calc_pvt!`](@ref) does.
"""
function empty_keeping_capacity! end

# An all-zero slot table with no hashes, keys, values or holes is exactly the empty state
# Dictionaries' own `empty!` builds, only with the slot table at its current size — which
# is a power of two, as the hashing requires, because only Dictionaries itself ever sized
# it.
#
# This reaches into the fields of `Dictionaries.Indices` (0.4), which is why it is
# pinned by a test of its own in `test/allocation_free.jl`.
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
    empty!(solution.inter_system_biases)
    empty!(solution.inter_frequency_biases)
    solution
end

"""
    get_lla(pvt::PVTSolution) -> LLA

Convert the ECEF position in the PVT solution to geodetic coordinates
(latitude, longitude, altitude) using the WGS84 ellipsoid.
"""
function get_lla(pvt::PVTSolution)
    LLAfromECEF(wgs84)(pvt.position)
end
