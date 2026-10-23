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
    get_satellite_info,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    calc_elevation_and_azimuth

# The measurement-model surface. `calc_pvt` above is the whole scalar solver,
# but a consumer running its own estimator over the same measurement model — a
# navigation filter closing tracking loops through its own Kalman update —
# consumes the model in pieces: the per-satellite times and clock rates, the
# bias layout and design-matrix columns, the predicted ranges and geometry, the
# atmospheric and time-system corrections. Those pieces stay unexported — they
# are solver internals, not names every `using` should carry — but they are
# documented as a surface (see "The Measurement-Model Surface" in the API
# reference), and a consumer binds them explicitly with
# `using PositionVelocityTime: calc_transmit_time, …`, which declares the
# dependency at a single site.

# Naming conventions. Names a user meets — exports, documented functions and their
# arguments, struct fields — are spelled out ("satellite", not "sat"); internals (locals,
# private helpers' arguments, `PVTWorkspace`'s fields) may keep short forms. `get_` reads
# or converts what its argument holds (`get_week`, `get_lla`), `calc_` computes, and
# physical models and named quantities keep noun names (`ionospheric_delay`,
# `day_of_year`). `ξ` and `H` are the textbook least-squares state and design matrix.

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

# PDOP above which a previous solution is distrusted as a warm-start seed and
# discarded (see the warm-start seeding in `_solve_pvt!`). Genuine geometries this
# side of useless stay well under it — a PDOP of 20 is already an unusable fix —
# while the spurious far-away roots the gate exists to catch show hundreds (the
# observed incident: 587). Discarding a fix that honestly earned a high DOP is
# harmless — the cold solve lands in the same place — so the exact value only
# sets how often that happens.
const MAX_TRUSTED_WARM_START_PDOP = 50.0

"""
    ClockConstraint

One known difference between two receiver clock biases, fitted as an extra observation
row: `ξ[3 + column] − ξ[3 + target_column] = offset`. It is how a broadcast GNSS time
offset (GGTO, BGTO) ties the clocks of two time systems together; see
[`decide_bias_layout`](@ref).

# Fields
- `column::Int`: the clock column of the system whose broadcast offset is used.
- `target_column::Int`: the clock column of the system the offset is toward.
- `offset::Float64`: the known difference in metres
  ([`broadcast_inter_system_bias`](@ref)).
"""
struct ClockConstraint
    column::Int
    target_column::Int
    offset::Float64
end

"""
    BiasColumns

Per-satellite assignment of the estimated bias columns of the least-squares design matrix,
shared by [`predict_pseudoranges!`](@ref), [`calc_design_matrix!`](@ref) and
[`estimate_position`](@ref). The state vector is `[x, y, z, tc₁, …, tc_num_clock_biases,
ifb₁, …, ifb_num_ifb]` with [`num_lsq_params`](@ref)`(bias_columns)` entries. The two column
kinds have different physical sources: a clock column is the receiver clock for one GNSS
time system (the spacing *between* systems is a system/space-segment effect — the broadcast
GNSS time offset, GGTO or BGTO), whereas an inter-frequency-bias column is the receiver's
per-band RF-chain delay. Known per-satellite corrections (atmosphere) are applied to the
pseudoranges in [`calc_pvt`](@ref), not carried here.

# Fields
- `clock_bias_indices::Vector{Int}`: per satellite, the clock column (1…`num_clock_biases`)
  of its GNSS time system; the design-matrix `1.0` lands at `3 + clock_bias_indices[j]`.
- `num_clock_biases::Int`: number of clock biases (also the offset of the IFB block).
- `ifb_indices::Vector{Int}`: per satellite, the inter-frequency-bias column
  (1…`num_ifb`) of its frequency band, or `0` for the reference band; the `1.0` lands
  at `3 + num_clock_biases + ifb_indices[j]`.
- `num_ifb::Int`: number of inter-frequency biases (frequency bands beyond the reference).
- `constraints::Vector{ClockConstraint}`: the clock-difference rows, appended after the
  satellite rows. Without the argument, an empty vector shared by every `BiasColumns`
  built that way and by [`dop_columns`](@ref), which must not be modified.
"""
struct BiasColumns
    clock_bias_indices::Vector{Int}
    num_clock_biases::Int
    ifb_indices::Vector{Int}
    num_ifb::Int
    constraints::Vector{ClockConstraint}
end

BiasColumns(clock_bias_indices, num_clock_biases, ifb_indices, num_ifb) = BiasColumns(
    clock_bias_indices, num_clock_biases, ifb_indices, num_ifb, NO_CLOCK_CONSTRAINTS)

# The constraints of every `BiasColumns` without any: shared, so building one allocates
# nothing for them. Never pushed to.
const NO_CLOCK_CONSTRAINTS = ClockConstraint[]

"""
    num_lsq_params(bias_columns::BiasColumns) -> Int

Length of the least-squares state vector for `bias_columns`:
`3 + num_clock_biases + num_ifb`.
"""
num_lsq_params(bias_columns::BiasColumns) =
    3 + bias_columns.num_clock_biases + bias_columns.num_ifb

"""
    SatelliteState{CP<:Real,D<:GNSSDecoderState,S<:AbstractGNSSSignal}

Combines the GNSS decoder state with code and carrier phase measurements for a single satellite.

# Fields
- `decoder::GNSSDecoderState`: GNSS decoder state containing decoded navigation data
- `system::AbstractGNSSSignal`: GNSS system (e.g., `GPSL1CA()`, `GalileoE1B()`)
- `code_phase::CP`: Code phase measurement in chips
- `carrier_doppler`: Carrier Doppler frequency in Hz
- `carrier_phase::CP`: Carrier phase measurement in radians, matching
  `Tracking.get_carrier_phase` (default: `0.0`)

# Constructors
    SatelliteState(; decoder, system, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, system, tracked_satellite)

The second constructor extracts code phase, carrier Doppler, and carrier phase from
`tracked_satellite`, a `Tracking.TrackedSat`. Tracking.jl defines it (from
Tracking 9.1 on), so it is available once `Tracking` is loaded; this package does not
depend on Tracking.
"""
@kwdef struct SatelliteState{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    system::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
end

# The signal-group container and the flat per-satellite measurement row every function
# below works on, plus the collection pass that turns the former into the latter.
# Included here rather than with the rest at the bottom of the file because the solver
# methods further down annotate their arguments with these types.
include("measurement.jl")
# The trim-safe least-squares fit behind `estimate_position`, included this early because
# the solver's workspace below holds its buffers.
include("levenberg_marquardt.jl")
using .LevenbergMarquardt: curve_fit!, LMWorkspace, grown

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
- `value::typeof(1.0m)`: The bias (metres), relative to `reference_band` — how much
  longer this band's receiver chain delay is than the reference band's.
- `reference_band::Symbol`: The reference band (e.g. `:L1`; see `get_band_id`) whose delay is
  folded into the clock biases and against which `value` is measured. Chosen per
  coverage component (see [`band_ifb_layout`](@ref)), so different biases in a
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
- `time::Union{TAITime, Nothing}`: Estimated time as a TAI epoch (see [`TAITime`](@ref);
  with AstroTime loaded, `TAIEpoch(time)` converts it)
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `satellites::Dictionary{Tuple{Symbol, Int}, SatelliteInfo}`: Maps `(signal, PRN)` to
  satellite info — position, transmit time, and the post-fit pseudorange and range-rate
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
  broadcast GGTO/BGTO fixes it only where the satellites leave it unobservable (see
  [`decide_bias_layout`](@ref)), as that broadcast value may be erroneous. Reference-independent
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
  [`band_ifb_layout`](@ref)), so a disconnected constellation can in principle yield
  several references. Empty for a single-band solution. A solution then needs
  `n ≥ 3 + M + B` satellites for `M` time systems and `B` extra bands, less one per
  broadcast offset that ties two clocks together (see [`decide_bias_layout`](@ref)).

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
    get_satellite_info(pvt_solution::PVTSolution, signal::Symbol, prn::Integer)
        -> Union{SatelliteInfo,Nothing}

Return the [`SatelliteInfo`](@ref) (position, transmit time and the post-fit
pseudorange/range-rate residuals) of the
satellite with the given `prn` on GNSS `signal` (e.g. `:GPSL1CA`, `:GalileoE1B`; see
`get_signal_id`), or `nothing` if that satellite was not used in the fix. The
signal tag is required because the same PRN can belong to different constellations or
be tracked on several signals; see the `satellites` field of [`PVTSolution`](@ref).
"""
function get_satellite_info(pvt_solution::PVTSolution, signal::Symbol, prn::Integer)
    get(pvt_solution.satellites, (signal, Int(prn)), nothing)
end

"""
    calc_elevation_and_azimuth(user_position::ECEF, satellite_position)
        -> (elevation, azimuth)
    calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, satellite_position)
        -> (elevation, azimuth)

Elevation and azimuth (radians) of the satellite at `satellite_position` (ECEF) as seen
from the user position, in the local East-North-Up frame there. Azimuth is measured
clockwise from North and wrapped to `[0, 2π)`, the GNSS convention.

Pass a precomputed `enu_from_ecef = ENUfromECEF(user_position, wgs84)` to reuse one
transform across all satellites that share a user position (one geodetic conversion per
epoch instead of per satellite).
"""
function calc_elevation_and_azimuth(enu_from_ecef::ENUfromECEF, satellite_position)
    enu = SphericalFromCartesian()(enu_from_ecef(ECEF(satellite_position)))
    elevation = enu.ϕ
    # `SphericalFromCartesian` measures θ counter-clockwise from East (the ENU +x
    # axis); azimuth is measured clockwise from North, i.e. π/2 − θ.
    azimuth = mod2pi(π / 2 - enu.θ)
    return elevation, azimuth
end

calc_elevation_and_azimuth(user_position::ECEF, satellite_position) =
    calc_elevation_and_azimuth(ENUfromECEF(user_position, wgs84), satellite_position)

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
    band_ifb_layout(system_keys, bands)
        -> (ifb_indices, extra_bands, reference_bands, num_components)

Lay out the receiver inter-frequency biases from the (constellation × band) coverage
graph. `system_keys[j]` and `bands[j]` are satellite `j`'s constellation and frequency
band; the constellation enters only through equality and first appearance, so any
per-satellite key that separates the constellations will do — [`decide_bias_layout`](@ref)
passes the *clock group* of each satellite (its clock column, unless a
[`ClockConstraint`](@ref) ties its system to another), which is concretely typed, unlike
the `TimeSystem` singletons themselves.

Two bands share a *coverage component* iff some constellation is tracked on
both (directly or through a chain of shared constellations); within a component the
clock/IFB split has a single gauge freedom, so one reference band per component (the
most-populated, ties by first appearance) is fixed and an IFB column is created only
for the *other* bands of that component. This estimates exactly the observable IFBs —
a band that is the sole band of its component (its constellation lives only there)
gets none, its delay folding into that constellation's clock — so the resulting design
matrix is always full rank.

`ifb_indices[j]` is satellite `j`'s IFB column (1…`length(extra_bands)`), or `0` for a
per-component reference band; `extra_bands[i]` is the band of IFB column `i` (ordered
deterministically) and `reference_bands[i]` is the reference band of that column's
coverage component (the anchor its IFB is measured against); `num_components` is the
number of coverage components (`1` ⇔ the graph is connected).
"""
function band_ifb_layout(system_keys, bands)
    ifb_indices = Int[]
    extra_bands = Vector{eltype(bands)}()
    reference_bands = Vector{eltype(bands)}()
    num_components = band_ifb_layout!(ifb_indices, extra_bands, reference_bands,
        BandLayoutScratch{eltype(bands)}(), system_keys, bands)
    return ifb_indices, extra_bands, reference_bands, num_components
end

# The scratch vectors of `band_ifb_layout!`. Bands are identified by their position in
# `unique_bands` (first appearance), so the union-find, the counts and the per-component
# references are all plain `Int` vectors — the same bookkeeping a `Dict` keyed by band
# would hold, without a hash table to allocate per epoch.
struct BandLayoutScratch{B}
    unique_bands::Vector{B}
    band_index::Vector{Int}   # per satellite: its band's position in `unique_bands`
    parent::Vector{Int}       # per band: union-find parent
    band_count::Vector{Int}   # per band: satellites on it
    reference::Vector{Int}    # per union-find root: the component's reference band, or 0
    column::Vector{Int}       # per band: its IFB column, or 0 for a reference band
end

BandLayoutScratch{B}() where {B} = BandLayoutScratch{B}(B[], Int[], Int[], Int[], Int[], Int[])

"""
    band_ifb_layout!(ifb_indices, extra_bands, reference_bands, scratch, system_keys, bands)
        -> num_components

[`band_ifb_layout`](@ref) writing its three vectors into the given ones (emptied first)
and working in `scratch`, so that it allocates nothing once they have grown to the epoch.
The result is identical: the same components, the same references and tie-breaks, and
the same column order.
"""
function band_ifb_layout!(
    ifb_indices,
    extra_bands,
    reference_bands,
    scratch::BandLayoutScratch,
    system_keys,
    bands,
)
    (; unique_bands, band_index, parent, band_count, reference, column) = scratch
    num_sats = length(bands)
    empty!(unique_bands)
    resize!(band_index, num_sats)
    for j in 1:num_sats
        u = findfirst(isequal(bands[j]), unique_bands)
        if isnothing(u)
            push!(unique_bands, bands[j])
            u = length(unique_bands)
        end
        band_index[j] = u
    end
    num_bands = length(unique_bands)
    # Union-find over bands: union the bands a single constellation is tracked on, by
    # linking every satellite's band to the band of its constellation's first satellite.
    resize!(parent, num_bands)
    parent .= 1:num_bands
    for j in 1:num_sats
        first_of_system = something(findfirst(i -> isequal(system_keys[i], system_keys[j]), 1:j))
        a = find_root!(parent, band_index[first_of_system])
        c = find_root!(parent, band_index[j])
        a == c || (parent[a] = c)
    end
    resize!(band_count, num_bands)
    fill!(band_count, 0)
    for j in 1:num_sats
        band_count[band_index[j]] += 1
    end
    # Reference band per component = most-populated in the component (ties: first seen).
    resize!(reference, num_bands)
    fill!(reference, 0)
    for u in 1:num_bands
        r = find_root!(parent, u)
        if reference[r] == 0 || band_count[u] > band_count[reference[r]]
            reference[r] = u
        end
    end
    empty!(extra_bands)
    empty!(reference_bands)
    resize!(column, num_bands)
    for u in 1:num_bands
        reference_band = reference[find_root!(parent, u)]
        if reference_band == u
            column[u] = 0
        else
            push!(extra_bands, unique_bands[u])
            push!(reference_bands, unique_bands[reference_band])
            column[u] = length(extra_bands)
        end
    end
    resize!(ifb_indices, num_sats)
    for j in 1:num_sats
        ifb_indices[j] = column[band_index[j]]
    end
    # Every root, and only a root, holds its component's reference.
    return count(!iszero, reference)
end

# The union-find root of `b` in `parent`, compressing the path to it. A loop rather
# than the recursive closure it could be written as: a closure that calls itself is
# boxed, which makes every call through it a dynamic dispatch.
function find_root!(parent, b)
    root = b
    while parent[root] != root
        root = parent[root]
    end
    while parent[b] != root
        b, parent[b] = parent[b], root
    end
    return root
end

"""
    BiasLayout()

The least-squares bias layout of one epoch, as [`decide_bias_layout!`](@ref) decides it.
Its vectors are refilled in place, so a layout reused across epochs allocates nothing
once it has grown to the largest one.

# Fields
- `clock_bias_indices::Vector{Int}`: per satellite, the clock column of its time system.
- `clock_systems::Vector{SupportedTimeSystem}`: the time system of each clock column, in
  order of first appearance.
- `constraints::Vector{ClockConstraint}`: the broadcast-offset rows; usually empty.
- `clock_group_indices::Vector{Int}`: per satellite, its *clock group*, in which the systems
  `constraints` tie together count as one.
- `num_clock_groups::Int`: the number of clock groups.
- `ifb_indices::Vector{Int}`, `extra_bands::Vector{Symbol}`,
  `reference_bands::Vector{Symbol}`: the inter-frequency-bias columns, as
  [`band_ifb_layout`](@ref) lays them out over the clock groups — `extra_bands[i]` is
  the band of IFB column `i` and `reference_bands[i]` the reference band it is measured
  against.

[`fit_columns`](@ref) is the design the position solve fits, [`dop_columns`](@ref) the one
its DOP is taken from.
"""
mutable struct BiasLayout
    clock_bias_indices::Vector{Int}
    clock_systems::Vector{SupportedTimeSystem}
    constraints::Vector{ClockConstraint}
    clock_group_indices::Vector{Int}
    num_clock_groups::Int
    ifb_indices::Vector{Int}
    extra_bands::Vector{Symbol}
    reference_bands::Vector{Symbol}
    # Scratch.
    bands::Vector{Symbol}
    satellite_counts::Vector{Int}            # per clock column: distinct satellites
    candidates::Vector{ClockConstraint}      # one per pair of systems with an offset
    system_parent::Vector{Int}               # per clock column: union-find parent
    clock_group_of_root::Vector{Int}         # per union-find root: its group, or 0
    band_scratch::BandLayoutScratch{Symbol}
end

BiasLayout() = BiasLayout(Int[], SupportedTimeSystem[], ClockConstraint[], Int[], 0, Int[],
    Symbol[], Symbol[], Symbol[], Int[], ClockConstraint[], Int[], Int[],
    BandLayoutScratch{Symbol}())

"""
    fit_columns(layout::BiasLayout) -> BiasColumns

The [`BiasColumns`](@ref) the position solve fits: one clock column per time system, and
the layout's clock constraints as extra rows. Its vectors **are** the layout's.
"""
fit_columns(layout::BiasLayout) = BiasColumns(layout.clock_bias_indices,
    length(layout.clock_systems), layout.ifb_indices, length(layout.extra_bands),
    layout.constraints)

"""
    dop_columns(layout::BiasLayout) -> BiasColumns

The [`BiasColumns`](@ref) the DOP is taken from: one clock column per clock group and no
constraint rows, i.e. the constraints taken as exact. Its vectors **are** the layout's,
except `constraints`: an empty vector shared by every call, which must not be modified.
"""
dop_columns(layout::BiasLayout) = BiasColumns(layout.clock_group_indices,
    layout.num_clock_groups, layout.ifb_indices, length(layout.extra_bands),
    NO_CLOCK_CONSTRAINTS)

"""
    decide_bias_layout(measurements) -> Union{BiasLayout,Nothing}

Decide the least-squares bias layout (see [`BiasLayout`](@ref)) for the flat
[`SatelliteMeasurement`](@ref) rows of one epoch, or `nothing` when the epoch cannot be
solved.

Every time system present gets its own clock column. A broadcast time offset between two
of them — Galileo's GGTO toward GPS Time, GPS's CNAV/CNAV-2 GGTO toward Galileo System
Time, BeiDou's BGTO toward either (see [`time_offset_available`](@ref)) — can tie them as
a [`ClockConstraint`](@ref) row, but only where the satellites cannot separate them:

1. Without any constraint, the layout is taken if the satellites suffice —
   `n ≥ 3 + M + B` measurements for `M` time systems and `B` inter-frequency biases, of
   which `3 + M` from *distinct* satellites (a satellite tracked on several bands supplies
   one line of sight) — and the constellation × band coverage graph is connected, so
   every inter-frequency bias is separable from the inter-system biases.
2. Otherwise the smallest set of constraints that achieves both is added. Systems a
   constraint ties count as one in both conditions, and a clock group's bands are
   connected through it. Among sets of one size, those tying the systems with the fewest
   distinct satellites are preferred (ties: the better-observed partner, then first
   appearance).
3. If no set connects the coverage graph, the smallest set that satisfies the counts is
   used (none, if the satellites already suffice). No IFB column is created for a band
   stranded on its own constellation, so its inter-frequency bias folds into that
   constellation's clock and the reported inter-system bias reads `GGTO + (IFB
   difference)`. Position, residuals and DOP are unaffected.

Failing all three, the epoch is unsolvable. Each constraint pins an otherwise
unobservable direction, so the fit meets it exactly and the constrained system's
inter-system bias reads as the broadcast `−c·Δt`. Between GPS and Galileo, Galileo's GGTO
(specified to 20 ns, 95 %) is preferred over GPS's. A constraint is evaluated at the
transmit time of the first satellite carrying it ([`broadcast_inter_system_bias`](@ref)).

The conditions are structural, and so necessary but not sufficient: the geometry can
still be degenerate, which `calc_pvt` catches from the DOP of the solved design.
"""
function decide_bias_layout(measurements)::Union{Nothing,BiasLayout}
    layout = BiasLayout()
    decide_bias_layout!(layout, measurements) ? layout : nothing
end

"""
    decide_bias_layout!(layout::BiasLayout, measurements) -> Bool

[`decide_bias_layout`](@ref) into `layout`: `true` and the layout stored there, or `false`
when the epoch cannot be solved. Allocates nothing once `layout` has grown to the epoch.

After `false`, `layout` is still a consistent layout of the epoch, just not a solvable one.
It ties every pair of systems that a broadcast offset can tie, through the most preferred
offsets that tie a system to one it is not tied to yet (a spanning forest over the
systems). Its `constraints` are those offsets, or none if no offset is available, and its
groups and inter-frequency-bias columns go with them. That is for a consumer, such as a
navigation filter, that keeps the clocks tied through an epoch it cannot solve.
"""
function decide_bias_layout!(layout::BiasLayout, measurements)
    num_sats = length(measurements)
    bands = resize!(layout.bands, num_sats)
    for j in 1:num_sats
        bands[j] = measurements[j].band_id
    end
    clock_systems = unique_time_systems!(layout.clock_systems,
        measurement.time_system for measurement in measurements)
    num_systems = length(clock_systems)
    clock_bias_indices = resize!(layout.clock_bias_indices, num_sats)
    satellite_counts = fill!(resize!(layout.satellite_counts, num_systems), 0)
    for j in 1:num_sats
        # `something`: every system is in `clock_systems` by construction, and saying so
        # keeps the columns a `Vector{Int}` rather than widening them with `Nothing`.
        column = something(time_system_index(clock_systems, measurements[j].time_system))
        clock_bias_indices[j] = column
        is_new_satellite(measurements, j) && (satellite_counts[column] += 1)
    end

    # Most epochs need no constraint, so that is tried before any offset is evaluated.
    candidates = empty!(layout.candidates)
    solvable, connected = group_clocks!(layout, measurements, 0)
    solvable && connected && return set_constraints!(layout, 0)

    # One candidate per pair of systems with a broadcast offset between them, from the
    # first satellite that decoded it (the offset is constellation-wide). Galileo's GGTO
    # is tried before GPS's.
    for a in 1:num_systems, b in (a+1):num_systems
        own, other = clock_systems[b] === GST() ? (b, a) : (a, b)
        push_offset_candidate!(candidates, measurements, clock_systems, own, other) ||
            push_offset_candidate!(candidates, measurements, clock_systems, other, own)
    end
    sort!(candidates; alg = InsertionSort, by = candidate -> (
        min(satellite_counts[candidate.column], satellite_counts[candidate.target_column]),
        -max(satellite_counts[candidate.column], satellite_counts[candidate.target_column]),
    ))

    # Sets of candidates as bit masks, candidate `i` at bit `i − 1`. Among sets of one
    # size, running the masks upwards visits first those that leave the least preferred
    # candidates out. A set never needs more constraints than there are systems to join.
    num_candidates = length(candidates)
    fallback = solvable ? 0 : -1
    for num_constraints in 1:min(num_candidates, num_systems - 1),
        mask in 0:(1 << num_candidates - 1)

        count_ones(mask) == num_constraints || continue
        solvable, connected = group_clocks!(layout, measurements, mask)
        solvable || continue
        connected && return set_constraints!(layout, mask)
        fallback < 0 && (fallback = mask)
    end
    if fallback >= 0
        group_clocks!(layout, measurements, fallback)
        return set_constraints!(layout, fallback)
    end
    # Unsolvable: every offset that ties two systems not tied yet, for a consumer that
    # keeps the clocks tied through the epoch.
    mask = spanning_constraints(layout)
    group_clocks!(layout, measurements, mask)
    set_constraints!(layout, mask)
    false
end

is_selected(mask, i) = (mask >> (i - 1)) & 1 == 1

# The candidates of `layout`, most preferred first, that each tie two systems no earlier
# one has tied: the largest set without a redundant offset, as a mask.
function spanning_constraints(layout::BiasLayout)
    parent = resize!(layout.system_parent, length(layout.clock_systems))
    parent .= eachindex(parent)
    mask = 0
    for (i, candidate) in enumerate(layout.candidates)
        a = find_root!(parent, candidate.column)
        b = find_root!(parent, candidate.target_column)
        a == b && continue
        parent[a] = b
        mask |= 1 << (i - 1)
    end
    mask
end

# Tie the systems of the candidates in `mask` into clock groups and lay out the IFB columns
# over them, into `layout`. Returns `(solvable, connected)`: whether the satellites suffice
# for the remaining unknowns, and whether the coverage graph is connected. A set that closes
# a loop holds a redundant offset, and is never solvable.
function group_clocks!(layout::BiasLayout, measurements, mask)
    num_systems = length(layout.clock_systems)
    parent = resize!(layout.system_parent, num_systems)
    parent .= 1:num_systems
    for (i, candidate) in enumerate(layout.candidates)
        is_selected(mask, i) || continue
        a = find_root!(parent, candidate.column)
        b = find_root!(parent, candidate.target_column)
        a == b && return false, false
        parent[a] = b
    end
    clock_group_of_root = fill!(resize!(layout.clock_group_of_root, num_systems), 0)
    clock_group_indices = resize!(layout.clock_group_indices, length(measurements))
    num_clock_groups = 0
    for j in eachindex(clock_group_indices)
        root = find_root!(parent, layout.clock_bias_indices[j])
        if clock_group_of_root[root] == 0
            num_clock_groups += 1
            clock_group_of_root[root] = num_clock_groups
        end
        clock_group_indices[j] = clock_group_of_root[root]
    end
    layout.num_clock_groups = num_clock_groups
    # Keyed by clock group, so the tied systems count as one constellation.
    num_components = band_ifb_layout!(layout.ifb_indices, layout.extra_bands,
        layout.reference_bands, layout.band_scratch, clock_group_indices, layout.bands)
    # Both are necessary for a full-rank design, neither is sufficient: the geometry itself
    # can still be degenerate, which `calc_pvt` screens for once the design exists.
    solvable =
        length(measurements) >= 3 + num_clock_groups + length(layout.extra_bands) &&
        sum(layout.satellite_counts) >= 3 + num_clock_groups
    return solvable, num_components == 1
end

# The candidates in `mask`, as the constraints of `layout`.
function set_constraints!(layout::BiasLayout, mask)
    constraints = empty!(layout.constraints)
    for (i, candidate) in enumerate(layout.candidates)
        is_selected(mask, i) && push!(constraints, candidate)
    end
    true
end

# Push the constraint of clock column `column` toward `target_column`, from the first
# row of the former's system that broadcasts that offset; `false` if no row does.
function push_offset_candidate!(candidates, measurements, clock_systems, column,
    target_column)
    row = offset_row(measurements, clock_systems[column], clock_systems[target_column])
    row == 0 && return false
    offset = broadcast_inter_system_bias(measurements[row], clock_systems[target_column])
    push!(candidates, ClockConstraint(column, target_column, offset))
    true
end

# The first row of `system` that carries a broadcast offset toward `target`, or 0.
function offset_row(measurements, system, target)
    for j in eachindex(measurements)
        measurements[j].time_system === system &&
            time_offset_available(measurements[j], target) && return j
    end
    0
end

# Identity-based helpers over the `time_system` field. Every
# `GNSSSignals.TimeSystem` is a singleton, so `===` is both the exactly right
# comparison and the one that compiles to a pointer test — cheaper than `==` or a
# `Dict` lookup, which branch on the field's `Union` (see `SupportedTimeSystem`) before
# they compare anything.

"""
    unique_time_systems(systems) -> Vector{SupportedTimeSystem}

The distinct GNSS time systems of `systems` (an iterable of `TimeSystem`s), in order of
first appearance — the order that fixes the clock columns of [`BiasColumns`](@ref).
"""
unique_time_systems(systems) = unique_time_systems!(SupportedTimeSystem[], systems)

# `unique_time_systems` into `unique_systems`, emptied first.
function unique_time_systems!(unique_systems, systems)
    empty!(unique_systems)
    for system in systems
        any(other -> other === system, unique_systems) || push!(unique_systems, system)
    end
    unique_systems
end

"""
    time_system_index(systems, system) -> Union{Int,Nothing}

Position of `system` in `systems`, or `nothing` if absent.
"""
time_system_index(systems, system) = findfirst(other -> other === system, systems)

# Whether row `j` is the first of its physical satellite, identified by `(time system,
# PRN)` — a PRN is unique only within its GNSS, and a satellite tracked on several bands
# contributes several rows but one line of sight.
#
# A scan of the rows already in hand, rather than a set of keys: a
# `(SupportedTimeSystem, Int)` key is not a concrete type, so collecting them would
# heap-allocate one box per satellite on every epoch. Quadratic in the satellite count,
# which for the dozens of rows an epoch holds is the cheaper of the two.
function is_new_satellite(measurements, j)
    measurement = measurements[j]
    for i in firstindex(measurements):(j-1)
        measurements[i].prn == measurement.prn &&
            measurements[i].time_system === measurement.time_system && return false
    end
    true
end

"""
    select_reference_system(layout::BiasLayout) -> system

The reference time system of the epoch `layout` was decided for (by
[`decide_bias_layout!`](@ref), whether or not it found the epoch solvable): the one whose
clock bias, reference time, week and start epoch define the reported time
([`PVTSolution`](@ref)'s `reference_system`). It is the system with the most distinct
satellites, ties going to the one that appears first. Throws an `ArgumentError` for a
layout without a time system. Allocates nothing.

Satellites rather than rows, because a satellite tracked on several bands is one line of
sight: its extra rows mostly determine inter-frequency biases, not its system's clock. It
is also the count `decide_bias_layout!` prefers to constrain the systems by, so the
reference system is never the one whose clock only a broadcast offset determines.
"""
function select_reference_system(layout::BiasLayout)
    isempty(layout.clock_systems) &&
        throw(ArgumentError("no time system to select a reference system from"))
    layout.clock_systems[reference_clock_column(layout)]
end

# The clock column of `select_reference_system`. `argmax` returns the first maximum, and
# the columns are in order of first appearance.
reference_clock_column(layout::BiasLayout) = argmax(layout.satellite_counts)

"""
    time_scale_offset_to_gpst(time_system) -> Float64

Signed offset of a GNSS time system's *count* against GPS Time's for the same
instant — `get_tai_offset(GPST) − get_tai_offset(time_system)`, so **negative**
where the system's count reads lower (BDT: `19 − 33 = −14.0`), and `0.0` for a
system that counts alike. `calc_time_scale_offsets` adds the *difference of two
of these* to a transmit time, which is what puts the +14 s onto a BeiDou time
in a GPS-referenced solve.

This is structural, not a bias: it follows from the time scales' definitions, not
from either system's steering. Both GPST and GST are `TAI − 19 s`, so they count
identically and this is `0.0` for every GPS and Galileo satellite. BDT is
`TAI − 33 s`, so a BeiDou second-of-week reads 14 s lower than the GPS
time-of-week for the same instant — BDT week 0 second 0 *is* GPS week 1356 time
of week 14, the 14 leap seconds that had accrued between the two epochs.

Derived from `GNSSSignals.get_tai_offset` rather than tabulated, so a constellation
added later is covered without touching this.

!!! note "Why this is a measurement correction and not a time-scale shift"

    It would be tempting to fold the 14 s into `get_time_of_week` and be done. That
    would be wrong: the satellite position is propagated from the *same* transmit
    time against the message's own `t_0e`, which is on the broadcasting system's
    scale, so shifting the reported time would move the 14 s into the ephemeris —
    about 55 km of along-track error at BeiDou MEO velocities. The transmit time
    must stay on its own scale and the correction must land on the pseudorange.
"""
time_scale_offset_to_gpst(time_system::GNSSSignals.TimeSystem) =
    ustrip(s, get_tai_offset(GPST()) - get_tai_offset(time_system))

"""
    calc_time_scale_offsets(measurements, reference_system) -> Vector{Float64}

Seconds to add to each satellite's transmit time to express it in
`reference_system`'s count, so that `calc_pseudoranges` may difference them.
`0.0` for every satellite of a system that counts alike — which is every GPS and
Galileo satellite, and all of them in a single-constellation epoch.

Relative to the *reference* system rather than to GPS Time, because the reference
system's count is also what dates the reported epoch: `reference_time` comes out
of the same differencing, and is combined with the reference system's week and
`system_start_epoch`. Anchoring on GPST instead leaves a BeiDou-referenced solve
reporting a GPS-count time of week dated from the BDT epoch — a 14 s error in
`PVTSolution.time` traded for the one this removes from the pseudoranges.

This is applied to the times *handed to the differencing*, and nowhere else. The
satellite position must still be propagated from the untouched transmit time
against the message's own `t_0e`, which is on the broadcasting system's scale, so
shifting the time itself would move the offset into the ephemeris — about 55 km
of along-track error at BeiDou MEO velocities. `SatelliteInfo.transmit_time` and
[`calc_steering_offset`](@ref) likewise keep the unconverted value.
"""
calc_time_scale_offsets(measurements, reference_system) = calc_time_scale_offsets!(
    Vector{Float64}(undef, length(measurements)), measurements, reference_system)

# `calc_time_scale_offsets` into `offsets`, resized to the satellite count.
function calc_time_scale_offsets!(offsets, measurements, reference_system)
    offsets = resize!(offsets, length(measurements))
    reference = time_scale_offset_to_gpst(reference_system)
    # Each row already carries its own anchor (`count_offset_to_gpst`), precomputed by
    # `collect_measurements`, so this reads a number rather than branching on the
    # `time_system` field per satellite.
    for (j, measurement) in enumerate(measurements)
        offsets[j] = reference - measurement.count_offset_to_gpst
    end
    offsets
end

"""
    broadcast_inter_system_bias(row, target) -> Float64

The inter-system bias (metres) of `row`'s time system against the time system `target`
that the broadcast offset predicts at `row`'s transmit time: `−c · Δt_systems`, with
`Δt_systems` the broadcast steering offset `row` carries toward `target` (see
[`calc_steering_offset`](@ref), which evaluates it at any other time). It is the receiver clock difference `tc_row − tc_target`
a [`ClockConstraint`](@ref) fixes, signed like [`PVTSolution`](@ref)'s
`inter_system_biases`. Zero when `row` is of `target` itself, `NaN` where `row` carries
no usable offset toward it (see [`time_offset_available`](@ref)).
"""
function broadcast_inter_system_bias(row::SatelliteMeasurement, target::SupportedTimeSystem)
    row.time_system === target && return 0.0
    -SPEED_OF_LIGHT *
    calc_steering_offset(row.time_offsets[time_offset_index(target)], row.transmit_time)
end

"""
    predict_atmospheric_delays(ξ, measurements, correction,
                               reference_time, doy, enable_tropospheric_correction) -> Vector{Float64}

Per-satellite ionospheric + tropospheric delay (metres), to be subtracted from the
pseudoranges. The user position is the first three elements of the least-squares
state vector `ξ = [x, y, z, tc₁, …]` (ECEF, metres); the remaining clock-bias
components are ignored. `correction` is the constellation-wide ionospheric model
from [`select_ionospheric_correction`](@ref) (`nothing` skips the ionosphere); the
troposphere uses the blind Saastamoinen zenith delays mapped by the Niell mapping
functions, whose seasonal term takes the day of year `doy` (see
[`tropospheric_delay`](@ref)), unless `enable_tropospheric_correction` is `false`.

A single corrected solve is enough: the delays depend on position only through the
satellite elevation/azimuth (and, for the troposphere, the user height), and
∂delay/∂position is negligible over the metre-level position uncertainty (a 15 m
shift moves the elevation by ~1e-5°), so delays predicted at a nearby position are
accurate to well under a millimetre — no iterate-to-convergence needed. The user
geodetic coordinates and the ENU transform depend only on `ξ`, so they are
computed once and reused across satellites.
"""
Base.@nospecializeinfer function predict_atmospheric_delays(
    ξ,
    measurements,
    @nospecialize(correction),
    reference_time,
    doy,
    enable_tropospheric_correction,
)::Vector{Float64}
    # Taken `@nospecialize`d, so that with a single method a caller holding `correction`
    # as `Any` resolves statically to this one unspecialised body, where an `Any`
    # argument to a specialising function would be a dynamic dispatch (which a
    # `juliac --trim` build rejects). The narrowing to each model's concrete type is
    # spelled out here instead; the solver itself does not come through here, it passes
    # an `IonosphericModel` straight to `predict_atmospheric_delays!`.
    #
    # The return annotation is load-bearing too: without it the `Any` argument makes the
    # return type `Any`.
    predict(model) = predict_atmospheric_delays!(
        Vector{Float64}(undef, length(measurements)), ξ, measurements,
        IonosphericModel(model), reference_time, doy, enable_tropospheric_correction)
    correction isa KlobucharParams && return predict(correction)
    correction isa BeiDouKlobucharParams && return predict(correction)
    correction isa NTCMGParams && return predict(correction)
    correction isa BDGIMParams && return predict(correction)
    isnothing(correction) && return predict(nothing)
    throw(ArgumentError("not an ionospheric correction: $(typeof(correction))"))
end

"""
    predict_atmospheric_delays!(delays, ξ, measurements, model::IonosphericModel,
                                reference_time, doy, enable_tropospheric_correction)
        -> delays

[`predict_atmospheric_delays`](@ref) into `delays` (resized to the satellite count), with
the ionospheric correction carried by an [`IonosphericModel`](@ref). Allocates nothing.
"""
function predict_atmospheric_delays!(
    delays,
    ξ,
    measurements,
    model,
    reference_time,
    doy,
    enable_tropospheric_correction,
)
    # The per-satellite loop specialises once per ionospheric model, behind this barrier:
    # the field is a five-member `Union`, one more than Julia splits by itself, so the
    # split is spelled out.
    predict(correction) = _predict_atmospheric_delays!(delays, ξ, measurements,
        correction, reference_time, doy, enable_tropospheric_correction)
    correction = model.correction
    correction isa KlobucharParams && return predict(correction)
    correction isa BeiDouKlobucharParams && return predict(correction)
    correction isa NTCMGParams && return predict(correction)
    correction isa BDGIMParams && return predict(correction)
    return predict(nothing)
end

function _predict_atmospheric_delays!(
    delays,
    ξ,
    measurements,
    correction,
    reference_time,
    doy,
    enable_tropospheric_correction,
)
    delays = resize!(delays, length(measurements))
    user_pos = ECEF(ξ[1], ξ[2], ξ[3])
    user_lla = LLAfromECEF(wgs84)(user_pos)
    enu_from_ecef = ENUfromECEF(user_pos, wgs84)
    for (j, measurement) in enumerate(measurements)
        elevation, azimuth = calc_elevation_and_azimuth(enu_from_ecef, measurement.position)
        iono = ionospheric_delay(
            correction,
            measurement.center_frequency,
            elevation,
            azimuth,
            user_lla,
            reference_time,
        )
        tropo =
            enable_tropospheric_correction ? tropospheric_delay(elevation, user_lla, doy) :
            0.0
        delays[j] = iono + tropo
    end
    delays
end

"""
    PVTWorkspace()

The reusable scratch storage of [`calc_pvt!`](@ref), in the order the solve uses it: the
flat measurement rows and their bias layout, the pseudoranges and the warm-start seed, the
satellite positions and the atmosphere-corrected pseudoranges the position solve fits, its
least-squares buffers and residuals, the design and normal-equations matrices of the DOP,
and the velocity solve's design rows and residuals. Every buffer grows to the largest epoch
it has seen and is reused from then on; its contents after a solve are not part of any
interface.

Satellite positions are kept as a `Vector{SVector{3,Float64}}` — a fixed-size quantity
per satellite, in a vector that is resized in place. The least-squares state, whose
length `3 + M + B` depends on the epoch's time systems and bands, is a plain vector
instead of a static one, so the solver compiles once rather than once per layout; the
matrices sized by it are kept at their largest and used through views.
"""
mutable struct PVTWorkspace
    measurements::Vector{SatelliteMeasurement}
    layout::BiasLayout
    pseudoranges::Vector{Float64}
    prev_ξ::Vector{Float64}
    sat_positions::Vector{SVector{3,Float64}}
    atmospheric_delays::Vector{Float64}
    corrected_ranges::Vector{Float64}
    least_squares::LMWorkspace
    residuals::Vector{Float64}
    design_matrix::Matrix{Float64}
    normal_matrix::Matrix{Float64}
    velocity_design::Vector{SVector{4,Float64}}
    rate_residuals::Vector{Float64}
end

PVTWorkspace() = PVTWorkspace(
    SatelliteMeasurement[],
    BiasLayout(),
    Float64[],
    Float64[],
    SVector{3,Float64}[],
    Float64[],
    Float64[],
    LMWorkspace(),
    Float64[],
    Matrix{Float64}(undef, 0, 0),
    Matrix{Float64}(undef, 0, 0),
    SVector{4,Float64}[],
    Float64[],
)

"""
    calc_pvt(groups, prev_pvt::PVTSolution = PVTSolution();
             approximate_year::Integer = year(now(UTC)),
             enable_ionospheric_correction::Bool = true,
             enable_tropospheric_correction::Bool = true) -> PVTSolution

Calculate Position, Velocity, and Time (PVT) from GNSS satellite measurements.

Satellites from different GNSS and frequency bands may be combined. Because each
constellation references its broadcasts to its own system time, one receiver clock
bias is estimated per GNSS time system; and because the receiver's RF chain delays
each band differently, one receiver inter-frequency bias is estimated per frequency
band beyond a reference band (shared across constellations on that band; see
`get_band_id`). The state vector is therefore
`[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]` for `M` distinct time systems and `B`
extra bands. Position and time are found by least squares; velocity and clock drift
are solved from carrier Doppler.

A solution requires `n ≥ 3 + M + B` healthy satellite measurements (each system needs at
least one satellite, and a system contributing a single satellite spends it entirely
on that system's clock bias), of which `3 + M` must come from *distinct* satellites: a
satellite tracked on several bands supplies one line of sight, and its extra
measurements constrain the inter-frequency biases rather than the geometry. When either
condition fails, a broadcast offset between two tracked systems — Galileo's GGTO
(Galileo–GPS Time Offset), BeiDou's BGTO (BDT–GNSS Time Offset, toward GPS or Galileo),
or GPS's CNAV/CNAV-2 GGTO — ties their clock biases together as an extra observation,
which makes a 4-satellite mixed fix possible. Offsets are used only where the satellites
cannot separate the clocks (see [`decide_bias_layout`](@ref)), so their own error never
shifts what the geometry determines.

Unless disabled via `enable_ionospheric_correction`, the ionospheric delay is
corrected automatically using only the coefficients decoded from the navigation
messages. A single model is chosen for the whole solve and applied to every
satellite, preferring the more accurate global TEC models: NTCM-G if Galileo
Effective Ionisation Level coefficients have been decoded, else BDGIM if a BDS-3
B-CNAV (B1C/B2a/B2b) coefficient set has, else the GPS Klobuchar model, else
BeiDou's own Klobuchar variant from its legacy B1I/B3I message, else no
correction. See
[`select_ionospheric_correction`](@ref) and [`ionospheric_delay`](@ref).

Unless disabled via `enable_tropospheric_correction`, the tropospheric delay is
corrected with a blind model (no broadcast coefficients needed): Saastamoinen
zenith delays mapped to the line of sight by the Niell mapping functions.
See [`tropospheric_delay`](@ref).

# Arguments
- `groups`: this epoch's measurements, as [`SignalGroups`](@ref) — a `NamedTuple` of
  [`SignalGroup`](@ref)s, one per ranging signal, each holding that signal's
  [`SatelliteState`](@ref)s. A bare `Tuple` of groups is accepted (and numbered
  `group1`, `group2`, …), as is a single `SignalGroup` (which becomes the `:default`
  group), so a one-constellation solve stays a one-liner:

  ```julia
  using PositionVelocityTime: SignalGroup
  calc_pvt(SignalGroup(GPSL1CA(), gps_states))
  calc_pvt((gps = SignalGroup(GPSL1CA(), gps_states),
            galileo = SignalGroup(GalileoE1B(), galileo_states)))
  ```

  Grouping is what makes the solve type-stable: within a group every satellite shares
  one concrete state type, so [`collect_measurements`](@ref) dispatches statically, and
  the solver behind it compiles once for every constellation mix. Build the groups where
  the satellites are tracked — with Tracking.jl 9.1 or later loaded,
  [`signal_groups`](@ref)`(track_state, decoders)` builds them from a whole `TrackState`.
  A pooled `Vector{SatelliteState}` is refused with an error saying what to build
  instead.

  Each `(signal, PRN)` pair must appear at most once — a receiver produces one
  measurement per signal per satellite, and a duplicate would enter the least-squares
  solve twice. A `Dictionary`-backed group prevents this within itself; two groups that
  share a ranging signal can still collide, and are not checked for.

  **Order is significant**: the flat measurement order is group order × within-group
  order, and it fixes the reference system's tie-break ("most satellites, ties by first
  appearance"), [`band_ifb_layout`](@ref)'s reference-band tie-break, and the insertion
  order of `PVTSolution.satellites`.
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin). A
  previous solution whose own DOP is implausible (`GDOP < 0`, or `PDOP` above
  `MAX_TRUSTED_WARM_START_PDOP` = 50) is not used as a seed and the epoch is
  solved from cold, so one spurious fix cannot re-seed itself epoch after epoch.

# Keyword Arguments
- `approximate_year`: Calendar year of the observation, used to resolve the
  GPS L1 C/A 1024-week rollover ambiguity (legacy LNAV broadcasts only a
  10-bit week number, so the receiver needs external information to
  determine which 1024-week cycle the recording is in). Anything within
  ±9 years of the actual observation date works. Defaults to the current
  UTC year, which is correct for live signals; for processing archived
  recordings, pass the rough year of the recording.
- `enable_ionospheric_correction`: when `true` (default), apply the automatic
  ionospheric correction described above. Set to `false` to skip it entirely
  and solve from the raw pseudoranges (e.g. for diagnostics or when an external
  correction is applied elsewhere).
- `enable_tropospheric_correction`: when `true` (default), apply the Saastamoinen
  tropospheric correction. Set to `false` to skip it.

# Returns
A new [`PVTSolution`](@ref) containing position, velocity, time, DOP values, and
satellite information; `prev_pvt` is never modified. [`calc_pvt!`](@ref) is the same
solve reusing the containers of a solution passed to it, without allocating. Returns
`prev_pvt` if the epoch cannot be solved: too few healthy satellites to solve the
constellation even with the broadcast-offset constraints (counted as distinct satellites
— a satellite tracked on several bands supplies one line of sight, so measurements alone
are not enough), or a geometry whose solved design matrix is rank deficient (reported as
a negative GDOP). None of these throw, so a receiver can pass
whatever it currently tracks each epoch and carry `prev_pvt` forward. Distrust of
`prev_pvt` (see above) affects only the seed: an unsolvable epoch still returns
`prev_pvt` exactly as passed.
"""
function calc_pvt(
    groups,
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    # A fresh solution lends the solve its (empty) containers; on an unsolvable epoch the
    # solver hands `prev_pvt` itself back.
    _, solution = _calc_pvt!(PVTSolution(), PVTWorkspace(), groups, prev_pvt,
        approximate_year, enable_ionospheric_correction, enable_tropospheric_correction)
    solution
end

"""
    calc_pvt!(solution::PVTSolution, workspace::PVTWorkspace, groups, prev_pvt::PVTSolution;
              approximate_year::Integer = year(now(UTC)),
              enable_ionospheric_correction::Bool = true,
              enable_tropospheric_correction::Bool = true) -> PVTSolution

[`calc_pvt`](@ref), reusing the containers of `solution` and the buffers of `workspace`,
so that a steady stream of epochs is solved without allocating. The solve is the same —
same arguments, same keywords, same numbers — and so is the result it returns.

[`PVTSolution`](@ref) is immutable, so the fix is a new solution, returned; keep the
return value:

```julia
pvt = PVTSolution()
workspace = PVTWorkspace()
for groups in epochs
    pvt = calc_pvt!(pvt, workspace, groups, pvt)
end
```

**`solution` is the output argument, and the only thing overwritten**: the returned
solution takes over its `satellites`, `inter_system_biases` and `inter_frequency_biases`
containers, which are emptied and refilled with this epoch's values. `solution` itself
is therefore spent — its own fields still describe the epoch it was computed for, but
its containers now hold the new one — so use the returned solution from here on, and
keep a [`copy`](@ref Base.copy(::PVTSolution)) of `solution` if you need it after the
call (a `Dictionary` made by `map` over `solution.satellites` shares its keys, and sees the
change too). `prev_pvt` is only read. Passing the same solution as both is allowed and
is the usual receiver loop above: every value the seed needs is read before anything is
written.

Where [`calc_pvt`](@ref) returns `prev_pvt` — an epoch that cannot be solved — so does
this: `prev_pvt` itself when it is `solution` (or shares its containers), and otherwise a
copy of it in `solution`'s containers, so that `prev_pvt` stays untouched.

It allocates nothing once `workspace` and `solution` have held an epoch with at least as
many satellites and as many estimated biases; a larger epoch grows them, once. That
covers everything from the collection pass to the returned solution, but not the
`groups` themselves, which the caller builds (and can reuse too — a `SignalGroup`'s
satellites may be any vector).

The `workspace` holds no state between epochs that affects the result — only buffers —
but it must not be used by two solves at once; give each task its own.
"""
function calc_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    groups,
    prev_pvt::PVTSolution;
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    solved, result = _calc_pvt!(solution, workspace, groups, prev_pvt, approximate_year,
        enable_ionospheric_correction, enable_tropospheric_correction)
    solved ? result : copy_into(solution, prev_pvt)
end

# The collection pass and the solve: `(true, fix)` with the fix in `solution`'s
# containers, or `(false, prev_pvt)` for an epoch that cannot be solved, with `solution`
# untouched. Specialises on the group shape, like `collect_measurements`, and is small;
# the solver behind it compiles once.
function _calc_pvt!(
    solution,
    workspace,
    groups,
    prev_pvt,
    approximate_year,
    enable_ionospheric_correction,
    enable_tropospheric_correction,
)
    # The function barrier. `_collect_measurements!` specialises on the group shape and
    # reduces it to flat rows; `_solve_pvt!` below is one compiled body for every
    # constellation mix there is.
    ionospheric_correction =
        _collect_measurements!(workspace.measurements, groups, approximate_year)
    _solve_pvt!(
        solution,
        workspace,
        IonosphericModel(enable_ionospheric_correction ? ionospheric_correction : nothing),
        prev_pvt,
        enable_tropospheric_correction,
    )
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

# The solver: ~200 lines that must compile exactly once, so every argument is of a
# concrete type — the ionospheric model included, which is wrapped in an
# `IonosphericModel` rather than passed as the `Union` it is (see there). The inner
# barrier that does specialise on the model, once per model, where the work actually is,
# sits behind `predict_atmospheric_delays!`.
#
# Returns `(true, fix)`, the fix a new solution holding `solution`'s containers, or
# `(false, prev_pvt)` for an epoch that cannot be solved. The containers are written only
# once the epoch is known to be solvable, and only after everything is read from
# `prev_pvt`, which may be the same solution. A tuple rather than a
# `Union{Nothing,PVTSolution}`: the solution holds references, so such a `Union` would
# be returned boxed — an allocation per epoch.
function _solve_pvt!(
    solution::PVTSolution,
    workspace::PVTWorkspace,
    ionospheric_model,
    prev_pvt::PVTSolution,
    enable_tropospheric_correction::Bool,
)
    measurements = workspace.measurements
    num_sats = length(measurements)

    # Every satellite that was not decoded-complete-and-healthy has already been
    # dropped by the collection pass, and each surviving row carries the GNSSSignals
    # keys that drive the solution: `time_system` (`GPST()`/`GST()`/`BDT()`) groups the
    # receiver clock bias — one per time system, ordered by first appearance —
    # and `band_id` (`:L1`, `:L5`, …) groups the receiver inter-frequency bias — one
    # per band beyond a per-coverage-component reference. (`signal_id`, `:GPSL1CA` …, is
    # the per-signal `satellites` identity used below, not a grouping key.)
    #
    # Solvability is decided here: `decide_bias_layout!` keeps only observable IFBs and
    # adds broadcast-offset constraints where the geometry is disconnected or
    # under-determined. A degenerate geometry, which no count can see, is caught after
    # the solve by the DOP.
    layout = workspace.layout
    decide_bias_layout!(layout, measurements) || return (false, prev_pvt)
    (; clock_systems, extra_bands, reference_bands, constraints) = layout
    bias_columns = fit_columns(layout)
    (; clock_bias_indices, num_clock_biases) = bias_columns
    num_params = num_lsq_params(bias_columns)

    # Reference system — its clock bias, reference time, week and start epoch
    # define the reported time (see `select_reference_system`).
    reference_clock_index = reference_clock_column(layout)
    reference_system = clock_systems[reference_clock_index]
    # `something`, not a bare `findfirst`: every clock column has a satellite, so this
    # never fails, and saying so keeps the index an `Int` rather than a
    # `Union{Int,Nothing}` that would widen the DOP call below.
    reference_index = something(findfirst(==(reference_clock_index), clock_bias_indices))

    # The common reference cancels out of the reported time (the reference clock
    # bias absorbs it), so any latest-transmit-time reference works — but the times
    # must first be put on one count (`calc_time_scale_offsets!`).
    # A BeiDou second-of-week reads 14 s below a GPS time of week for the same instant,
    # so differencing them raw hands every BeiDou measurement 4.2e9 m of structural
    # offset: not merely a biased BeiDou clock column, but a parameter nine orders of
    # magnitude above the others in the normal equations. Zero for GPS and Galileo, and
    # for any single-constellation epoch.
    # The offsets are written straight into the pseudorange buffer, and the transmit
    # times added onto them, so no separate offset buffer is needed.
    pseudoranges =
        calc_time_scale_offsets!(workspace.pseudoranges, measurements, reference_system)
    for (j, measurement) in enumerate(measurements)
        pseudoranges[j] += measurement.transmit_time
    end
    _, reference_time = calc_pseudoranges!(pseudoranges)
    # The clock constraints' known offsets are observations too, after the satellites'.
    resize!(pseudoranges, num_sats + length(constraints))
    for (i, constraint) in enumerate(constraints)
        pseudoranges[num_sats+i] = constraint.offset
    end

    # The reference system's week and start epoch date the epoch absolutely: the day of
    # year feeds the tropospheric mapping's seasonal term here, and week/start epoch
    # date the reported time after the solve.
    reference_measurement = measurements[reference_index]
    week = reference_measurement.week
    start_time = reference_measurement.system_start_epoch
    doy = day_of_year(reference_measurement.system_start_time, week, reference_time)

    # Gauss-Newton converges within its seed's basin, so a spurious far-away root
    # would re-seed itself epoch after epoch. Such roots flag themselves with an
    # impossible geometry (negative GDOP, exploded PDOP): solve from a cold seed
    # instead. Only the seed is affected — an unsolvable epoch still returns
    # `prev_pvt` unchanged.
    distrusted =
        !isnothing(prev_pvt.dop) &&
        (prev_pvt.dop.GDOP < 0 || prev_pvt.dop.PDOP > MAX_TRUSTED_WARM_START_PDOP)

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias. A
    # distrusted previous solution seeds nothing: the state stays all zeros, the cold
    # start.
    #
    # A system the previous solution has no bias for (it was not tracked then) is seeded
    # with the previous reference clock. That is off by the inter-system bias, typically
    # metres; the alternative, zero, is off by the whole receiver clock bias, often
    # orders of magnitude more.
    prev_abs_bias(sys) =
        sys === prev_pvt.reference_system ? prev_pvt.clock_bias :
        prev_pvt.clock_bias + get(prev_pvt.inter_system_biases, sys, 0.0m)
    prev_ξ = fill!(resize!(workspace.prev_ξ, num_params), 0.0)
    if !distrusted
        prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
        for (k, sys) in enumerate(clock_systems)
            prev_ξ[3+k] = ustrip(m, prev_abs_bias(sys))
        end
        # Only warm-start an IFB column from the previous solution when that band's bias
        # was measured against the same reference band; a reference change (the anchor
        # differs across epochs) makes the stored value refer to a different quantity,
        # so seed at 0 instead. The IFB enters the design linearly, so a stale seed only
        # costs iterations, but this keeps the starting point meaningful.
        #
        # `haskey` and an index, not `get(…, nothing)`: an `InterFrequencyBias` holds a
        # `Symbol`, so a `Union` of it with `nothing` is returned boxed.
        prev_ifbs = prev_pvt.inter_frequency_biases
        for (i, band) in enumerate(extra_bands)
            prev_ξ[3+num_clock_biases+i] =
                haskey(prev_ifbs, band) && prev_ifbs[band].reference_band == reference_bands[i] ?
                ustrip(m, prev_ifbs[band].value) : 0.0
        end
    end

    # The satellite positions at transmission, from the collection pass.
    sat_positions = resize!(workspace.sat_positions, num_sats)
    for j in 1:num_sats
        sat_positions[j] = measurements[j].position
    end

    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model was chosen for the whole solve by the
    # collection pass. The prediction is `predict_atmospheric_delays!`, a function
    # barrier that specialises on the model this body deliberately does not.
    #
    # Is any atmospheric correction active at all? When neither the ionosphere (no
    # model selected) nor the troposphere contributes, skip the prediction entirely and
    # the solve runs on the raw pseudoranges.
    correct_atmosphere =
        !isnothing(ionospheric_model.correction) || enable_tropospheric_correction
    function corrected_ranges(ξ_prediction)
        delays = predict_atmospheric_delays!(
            workspace.atmospheric_delays, ξ_prediction, measurements, ionospheric_model,
            reference_time, doy, enable_tropospheric_correction)
        ranges = resize!(workspace.corrected_ranges, length(pseudoranges)) .= pseudoranges
        view(ranges, 1:num_sats) .-= delays
        ranges
    end
    least_squares = workspace.least_squares

    # `ξ` is the least-squares workspace's own parameter vector, so it stays valid until
    # the next solve on that workspace — which is after everything below.
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, residuals_uncorrected = estimate_position!(least_squares,
            workspace.residuals, sat_positions, pseudoranges, bias_columns, prev_ξ)
        if correct_atmosphere
            estimate_position!(least_squares, workspace.residuals, sat_positions,
                corrected_ranges(ξ_uncorrected), bias_columns, ξ_uncorrected)
        else
            (ξ_uncorrected, residuals_uncorrected)
        end
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction. With no
        # active correction, solve directly on the raw pseudoranges.
        estimate_position!(least_squares, workspace.residuals, sat_positions,
            correct_atmosphere ? corrected_ranges(prev_ξ) : pseudoranges, bias_columns,
            prev_ξ)
    end
    position = ECEF(ξ[1], ξ[2], ξ[3])

    # Check the geometry at the converged position — the DOP is reported to the caller, and
    # a rank deficiency here (negative GDOP) means `ξ` is meaningless, so nothing should be
    # derived from it. The satellite conditions in `decide_bias_layout!` already reject the
    # count-shaped degeneracies before the solve; this catches what a count cannot see.
    #
    # It also stands ahead of the velocity solve, whose design is this one's but for a few
    # 1e-6 relative; that solve tests its own design against the same tolerance as a
    # backstop and returns `NaN` should it fail, an unsolvable epoch as well (see
    # `estimate_velocity`).
    #
    # The DOP takes the constraints as exact (see `dop_columns`).
    dop_bias_columns = dop_columns(layout)
    num_dop_params = num_lsq_params(dop_bias_columns)
    workspace.design_matrix = grown(workspace.design_matrix, num_sats, num_dop_params)
    H = calc_design_matrix!(view(workspace.design_matrix, 1:num_sats, 1:num_dop_params),
        sat_positions, ξ, dop_bias_columns)
    workspace.normal_matrix = grown(workspace.normal_matrix, num_dop_params, num_dop_params)
    dop = calc_dop!(view(workspace.normal_matrix, 1:num_dop_params, 1:num_dop_params), H,
        position, layout.clock_group_indices[reference_index])

    dop.GDOP < 0 && return (false, prev_pvt)

    # The velocity solve takes one Gauss–Newton step, which is the least-squares fit to
    # rounding from the previous epoch's velocity and clock drift, and to ~1e-8 m/s from
    # rest (see `estimate_velocity`). So it starts from the previous
    # solution whenever the position solve trusts that as its seed, and from rest
    # otherwise: on a cold start (the default `prev_pvt` is at rest anyway), after a
    # distrusted solution, or should the previous one carry a non-finite value.
    previous_velocity_and_drift = SVector(prev_pvt.velocity[1], prev_pvt.velocity[2],
        prev_pvt.velocity[3], prev_pvt.relative_clock_drift * SPEED_OF_LIGHT)
    velocity_start =
        !distrusted && all(isfinite, previous_velocity_and_drift) ?
        previous_velocity_and_drift : zero(SVector{4,Float64})
    user_velocity_and_clock_drift, rate_residuals =
        estimate_velocity!(workspace.rate_residuals,
            workspace.velocity_design, measurements, position, velocity_start)
    all(isfinite, user_velocity_and_clock_drift) || return (false, prev_pvt)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    clock_bias = ξ[3+reference_clock_index]
    # A receiver clock running ahead reads the epoch late by its bias.
    corrected_reference_time = reference_time - clock_bias / SPEED_OF_LIGHT

    # Everything is read from `prev_pvt` by now, so the containers of `solution` —
    # possibly `prev_pvt`'s own — are written from here on.
    #
    # Per-satellite `satellites` key: (signal id, PRN) — signal-level (not time system), so
    # a satellite tracked on two signals of one constellation stays distinct; the
    # receiver-clock grouping is separate, by time system.
    #
    # A duplicate key is refused, as `Dictionary(keys, values)` refused it before. It is
    # checked before any container is written, so the throw leaves `solution` — possibly
    # `prev_pvt` itself — as it was, not half-refilled. Not with `insert!`, which formats
    # the key into its message — a dynamic `show` a `juliac --trim` build rejects — so
    # the message here is fixed.
    for j in 2:num_sats, k in 1:j-1
        measurements[j].signal_id === measurements[k].signal_id &&
            measurements[j].prn == measurements[k].prn && throw(IndexError(
                "a (signal, PRN) pair appears twice in one epoch; each must appear once"))
    end
    satellites = empty_keeping_capacity!(solution.satellites)
    for (j, measurement) in enumerate(measurements)
        key = (measurement.signal_id, measurement.prn)
        set!(satellites, key, SatelliteInfo(measurement.position, measurement.transmit_time,
            residuals[j] * m, rate_residuals[j] * (m / s)))
    end

    # Inter-system biases relative to the reference system's clock, in
    # meters. The reference is omitted (its bias is `clock_bias`).
    inter_system_biases = empty!(solution.inter_system_biases)
    for (k, sys) in enumerate(clock_systems)
        sys === reference_system && continue
        inter_system_biases[sys] = (ξ[3+k] - clock_bias) * m
    end

    # Receiver inter-frequency biases relative to the reference band, in meters
    # (the reference band is omitted; its bias is folded into the clock biases).
    inter_frequency_biases = empty!(solution.inter_frequency_biases)
    for (i, band) in enumerate(extra_bands)
        inter_frequency_biases[band] =
            InterFrequencyBias(ξ[3+num_clock_biases+i] * m, reference_bands[i])
    end

    return true,
    PVTSolution(
        position,
        velocity,
        calc_course_over_ground(position, velocity),
        clock_bias * m,
        # Assumes `start_time.fraction == 0` (true for GPS/Galileo: integer-second
        # origins).
        TAITime(
            week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) +
            start_time.second,
            corrected_reference_time - floor(Int, corrected_reference_time),
        ),
        user_velocity_and_clock_drift[4] / SPEED_OF_LIGHT,
        dop,
        satellites,
        reference_system,
        inter_system_biases,
        inter_frequency_biases,
    )
end

"""
    system_start_epoch(system) -> TAITime

Absolute TAI epoch of a ranging signal's GNSS time-scale origin (week 0, time of
week 0), from GNSSSignals' `get_tai_system_start_time` — the epoch already
labelled on the atomic scale (GPS `1980-01-06T00:00:19` TAI, Galileo
`1999-08-22T00:00:19` TAI, BeiDou `2006-01-01T00:00:33` TAI), so no leap-second
table is consulted, which is also what makes the value safe to derive at precompile
time.
"""
system_start_epoch(system) = TAITime(get_tai_system_start_time(system))

"""
    get_week(decoder::GNSSDecoderState{<:GPSL1CAData}; approximate_year)

Return the absolute GPS week number for a GPSL1 decoder, resolving the
1024-week rollover ambiguity using `approximate_year` as a calendar
anchor.

The legacy GPS L1 C/A LNAV message broadcasts only a 10-bit week number
(0–1023) modulo 1024, so the receiver cannot determine which 1024-week
cycle the recording is in from the data alone (IS-GPS-200, §20.3.3.3).
Each cycle is ~19.6 years, so any anchor within ±9 years of the true
observation date selects the correct cycle.

GPS week 0 is 1980-01-06; cycle boundaries fall on 1999-08-22,
2019-04-07, 2038-11-21, 2058-07-08, …

For Galileo, the broadcast WN is 12 bits and does not need this
treatment in any practical operational scenario.
"""
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.GPSL1CAData};
    approximate_year::Integer = year(now(UTC)),
)
    # GPS week 0 begins 1980-01-06. Compute the integer week count from
    # there to mid-`approximate_year`, then choose the cycle base such
    # that `cycle_base + WN` is closest to that anchor.
    days_at_anchor = Date(approximate_year, 6, 30) - Date(1980, 1, 6)
    weeks_at_anchor = Dates.value(days_at_anchor) ÷ 7
    n_cycles = round(Int, (weeks_at_anchor - decoder.data.WN) / 1024)
    return n_cycles * 1024 + decoder.data.WN
end

# Every other navigation message broadcasts the absolute week number — Galileo I/NAV
# and F/NAV 12-bit, GPS CNAV (L5, L2C) and CNAV-2 (L1C) 13-bit, BeiDou D1/D2 and
# B-CNAV 13-bit — so there is no rollover to resolve, `approximate_year` is unused,
# and one method serves all of them. GPS L1 C/A above is the sole exception, and the
# reason this function takes `approximate_year` at all.
#
# The Galileo half of the union is `AbstractGalileoEphemerisData`, not the wider
# `AbstractGalileoData`: `GalileoE6BData` is Galileo data with no week number at all
# (C/NAV stamps a time of hour instead), so the wider bound would put a method on it
# that raises a `FieldError` on the field it goes looking for. That is precisely the
# distinction GNSSDecoder introduced the narrower supertype to draw.
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{
        <:Union{
            GNSSDecoder.AbstractGalileoEphemerisData,
            GNSSDecoder.AbstractBeiDouData,
            AbstractGPSCNAVData,
        },
    };
    approximate_year::Integer = year(now(UTC)),
)
    decoder.data.WN
end

"""
    get_lla(pvt::PVTSolution) -> LLA

Convert the ECEF position in the PVT solution to geodetic coordinates
(latitude, longitude, altitude) using the WGS84 ellipsoid.
"""
function get_lla(pvt::PVTSolution)
    LLAfromECEF(wgs84)(pvt.position)
end

include("user_position.jl")
include("sat_time.jl")
include("sat_position.jl")
include("ionosphere.jl")
include("troposphere.jl")
# The precompile workload solves real fixtures, so it must come after every
# method it dispatches into.
include("precompile.jl")
end
