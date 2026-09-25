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

# PDOP above which a previous solution is distrusted as a warm-start seed and
# discarded (see the gate at the top of `calc_pvt`). Genuine geometries this
# side of useless stay well under it — a PDOP of 20 is already an unusable fix —
# while the spurious far-away roots the gate exists to catch show hundreds (the
# observed incident: 587). Discarding a fix that honestly earned a high DOP is
# harmless — the cold solve lands in the same place — so the exact value only
# sets how often that happens.
const MAX_TRUSTED_WARM_START_PDOP = 50.0

"""
    BiasColumns

Per-satellite assignment of the estimated bias columns of the least-squares design
matrix, shared by [`calc_ρ_hat!`](@ref), [`calc_H!`](@ref) and [`user_position`](@ref).
The state vector is `[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifb]` with
[`num_lsq_params`](@ref)`(bias_columns)` entries. The two column kinds have different
physical sources: a clock column is the receiver clock for one GNSS time system (the
spacing *between* systems is a system/space-segment effect — the broadcast GNSS time
offset, GGTO or BGTO), whereas an inter-frequency-bias column is the receiver's per-band
RF-chain delay. Known per-satellite corrections (atmosphere, the broadcast time-system
offset) are applied to the pseudoranges in [`calc_pvt`](@ref), not carried here.

# Fields
- `clock_bias_indices::Vector{Int}`: per satellite, the clock column (1…`num_clock_biases`)
  of its GNSS time system; the design-matrix `1.0` lands at `3 + clock_bias_indices[j]`.
- `num_clock_biases::Int`: number of clock biases (also the offset of the IFB block).
- `ifb_indices::Vector{Int}`: per satellite, the inter-frequency-bias column
  (1…`num_ifb`) of its frequency band, or `0` for the reference band; the `1.0` lands
  at `3 + num_clock_biases + ifb_indices[j]`.
- `num_ifb::Int`: number of inter-frequency biases (frequency bands beyond the reference).
"""
struct BiasColumns
    clock_bias_indices::Vector{Int}
    num_clock_biases::Int
    ifb_indices::Vector{Int}
    num_ifb::Int
end

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
    SatelliteState(decoder, system, sat_state)

The second constructor extracts code phase, carrier Doppler, and carrier phase from a
`Tracking` satellite state (`Tracking.TrackedSat`). It is provided by a package extension
that is loaded automatically once `Tracking` is available, so `Tracking` is only a weak
dependency of this package.
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
# The trim-safe least-squares fit behind `user_position`, included this early because the
# solver's workspace below holds its buffers.
include("levenberg_marquardt.jl")
using .LevenbergMarquardt: curve_fit, curve_fit!, LMWorkspace, grown

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
    SatInfo

Per-satellite information attached to a [`PVTSolution`](@ref) (one entry per
satellite used in the fix).

# Fields
- `position::ECEF{Float64}`: Satellite ECEF position at transmit time (metres).
- `time::Float64`: Satellite transmit time (system time of week, seconds).
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
`calc_user_velocity_and_clock_drift`, whose `yⱼ` sets that sign).
"""
struct SatInfo
    position::ECEF{Float64}
    time::Float64
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
- `reference::Symbol`: The reference band (e.g. `:L1`; see `get_band_id`) whose delay is
  folded into the clock biases and against which `value` is measured. Chosen per
  coverage component (see [`band_ifb_layout`](@ref)), so different biases in a
  disconnected solution can in principle carry different references.
"""
struct InterFrequencyBias
    value::typeof(1.0m)
    reference::Symbol
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
- `time_correction::typeof(1.0m)`: Estimated receiver clock bias of the reference GNSS
  (meters). For a multi-GNSS solution this is the bias of `reference_system`;
  other systems' biases are `time_correction + inter_system_biases[system]`.
- `time::Union{TAITime, Nothing}`: Estimated time as a TAI epoch (see [`TAITime`](@ref);
  with AstroTime loaded, `TAIEpoch(time)` converts it)
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `sats::Dictionary{Tuple{Symbol, Int}, SatInfo}`: Maps `(signal, PRN)` to satellite
  info — position, transmit time, and the post-fit pseudorange and range-rate
  residuals (see [`SatInfo`](@ref)). The signal tag (e.g. `:GPSL1CA`,
  `:GalileoE1B`; see `get_signal_id`) keeps the
  same PRN apart both across constellations (GPS PRN 5 vs Galileo E05) and across
  signals of one constellation (a satellite tracked on GPS L1 C/A and L5 yields two
  entries sharing a PRN). Receiver-clock grouping is by time system, not signal —
  see `reference_system`.
- `reference_system::Union{SupportedTimeSystem, Nothing}`: GNSS time system (e.g.
  `GNSSSignals.GPST()`, `GST()`) that `time` and `time_correction` are referenced to.
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
  [`band_ifb_layout`](@ref)), so a disconnected constellation can in principle yield
  several references. Empty for a single-band solution. A solution then needs
  `n ≥ 3 + M + B` satellites for `M` time systems and `B` extra bands.

Immutable: a solution's fields cannot be reassigned once it is built. The allocation-free
[`calc_pvt!`](@ref) returns a new solution too; what it reuses are the `sats`,
`inter_system_biases` and `inter_frequency_biases` containers of the solution it is
explicitly handed as its output, which the returned solution takes over.
[`calc_pvt`](@ref) builds new containers every time and never modifies the `prev_pvt`
it reads.
"""
@kwdef struct PVTSolution
    position::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    velocity::ECEF{Float64} = ECEF(0.0, 0.0, 0.0)
    course_over_ground::typeof(1.0°) = 0.0°
    time_correction::typeof(1.0m) = 0.0m
    time::Union{TAITime,Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    sats::Dictionary{Tuple{Symbol,Int},SatInfo} = Dictionary{Tuple{Symbol,Int},SatInfo}()
    reference_system::Union{SupportedTimeSystem,Nothing} = nothing
    inter_system_biases::Dict{SupportedTimeSystem,typeof(1.0m)} =
        Dict{SupportedTimeSystem,typeof(1.0m)}()
    inter_frequency_biases::Dict{Symbol,InterFrequencyBias} =
        Dict{Symbol,InterFrequencyBias}()
end

"""
    get_sat_info(pvt_solution::PVTSolution, signal::Symbol, prn::Integer) -> Union{SatInfo,Nothing}

Return the [`SatInfo`](@ref) (position, transmit time and the post-fit
pseudorange/range-rate residuals) of the
satellite with the given `prn` on GNSS `signal` (e.g. `:GPSL1CA`, `:GalileoE1B`; see
`get_signal_id`), or `nothing` if that satellite was not used in the fix. The
signal tag is required because the same PRN can belong to different constellations or
be tracked on several signals; see the `sats` field of [`PVTSolution`](@ref).
"""
function get_sat_info(pvt_solution::PVTSolution, signal::Symbol, prn::Integer)
    get(pvt_solution.sats, (signal, Int(prn)), nothing)
end

"""
    get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF) -> Spherical
    get_sat_enu(enu_from_ecef::ENUfromECEF, sat_pos_ecef::ECEF) -> Spherical

Convert satellite position to East-North-Up (ENU) spherical coordinates (azimuth and
elevation) relative to the user position.

# Arguments
- `user_pos_ecef::ECEF`: User position in ECEF coordinates
- `enu_from_ecef::ENUfromECEF`: A precomputed `ENUfromECEF(user_pos_ecef, wgs84)`
  transform. Pass this form to reuse one transform across all satellites that
  share a user position (one geodetic conversion per epoch instead of per satellite).
- `sat_pos_ecef::ECEF`: Satellite position in ECEF coordinates

# Returns
Spherical coordinates containing azimuth and elevation of the satellite as seen from
the user position.
"""
function get_sat_enu(enu_from_ecef::ENUfromECEF, sat_pos_ecef::ECEF)
    SphericalFromCartesian()(enu_from_ecef(sat_pos_ecef))
end

get_sat_enu(user_pos_ecef::ECEF, sat_pos_ecef::ECEF) =
    get_sat_enu(ENUfromECEF(user_pos_ecef, wgs84), sat_pos_ecef)

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
passes the *clock column* of each satellite, which is exactly that and, unlike the
`TimeSystem` singletons themselves, concretely typed.

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
    BiasLayout

The full least-squares bias layout of one epoch, as [`decide_bias_layout`](@ref) returns
it: a `NamedTuple` type spelled out so the decision has a concrete return type. That
matters beyond tidiness — the layout is destructured straight into the solver, which is
compiled once for every constellation mix, so an unresolved field type here would become
a dynamic dispatch on every epoch. See [`decide_bias_layout`](@ref) for what each field
means.
"""
const BiasLayout = @NamedTuple{
    bias_columns::BiasColumns,
    extra_bands::Vector{Symbol},
    reference_bands::Vector{Symbol},
    hub_system::Union{Nothing,SupportedTimeSystem},
    hub_rows::Vector{SatelliteMeasurement},
}

"""
    decide_bias_layout(measurements) -> Union{BiasLayout,Nothing}

Decide the full least-squares bias layout — one clock column per GNSS time system plus
the per-band inter-frequency-bias columns from [`band_ifb_layout`](@ref) — for the flat
[`SatelliteMeasurement`](@ref) rows of one epoch, and return it as
a [`BiasLayout`](@ref),
or `nothing` when the constellation cannot be solved. Every classification the decision
needs — each satellite's `time_system`, `band_id`, `prn` and `time_offsets` — is a
field of the row:

- `bias_columns::BiasColumns`: the per-satellite column assignment (see
  [`BiasColumns`](@ref)).
- `extra_bands::Vector{Symbol}`: band of each inter-frequency-bias column, so
  `extra_bands[i]` belongs to column `i` of the IFB block.
- `reference_bands::Vector{Symbol}`: per IFB column, the reference band of its coverage
  component — the anchor that column's bias is measured against (see
  [`band_ifb_layout`](@ref)).
- `hub_system::Union{Nothing,SupportedTimeSystem}`: the system the collapsed
  clocks were merged onto, or `nothing` for a layout that estimates every clock bias
  independently, which is the common case.
- `hub_rows::Vector{SatelliteMeasurement}`: one row per collapsed time system — the
  satellite whose broadcast offset to `hub_system` converts that whole system's
  measurements. Each row carries its own `time_system`, so it is its own key. The same
  vocabulary as [`time_offset_available`](@ref), which decides membership, and
  [`calc_steering_offset`](@ref), which evaluates one entry (see
  [`calc_hub_range_offsets`](@ref)). Empty exactly when `hub_system` is `nothing`.

Only the *classification* fields of each row are read, never its transmit time: the
counts, the coverage graph and the availability flags are all the decision needs. The
offsets its `hub_rows` enable are evaluated afterwards, at each satellite's own
transmit time.

The decision is observability-driven, not merely count-driven:

- When the coverage graph is connected and the satellites suffice for the unknowns
  (`n ≥ 3 + num_systems + num_ifb` measurements, of which `3 + num_systems` must come
  from *distinct* satellites — extra bands of an already-tracked satellite add
  inter-frequency-bias information, not geometry), estimate everything independently: the
  inter-system offset and the receiver inter-frequency biases are observed directly from
  the geometry, so neither inherits the broadcast-GGTO error (the satellite group delays
  are already removed per satellite upstream, so the per-band column carries the receiver
  chain).
- Otherwise merge every clock that can be merged onto a *hub* system, using a
  broadcast offset to that hub (see [`time_offset_available`](@ref)). Hubs are tried
  in the fixed order GPST, GST, BDT, and the first whose merged layout has enough
  satellites wins, so the choice is deterministic and a GPS-anchored collapse behaves
  exactly as it always did. Toward GPS Time the offsets on the air are Galileo's GGTO
  and BeiDou's BGTO; toward Galileo System Time they are BeiDou's BGTO variant and
  the GGTO GPS itself broadcasts on CNAV/CNAV-2 — which is what lets a GPS-free
  Galileo + BeiDou epoch, or a GPS + Galileo epoch whose Galileo satellites have not
  decoded their GGTO yet, still collapse. The merge removes a clock unknown per
  merged system (the scarce-satellite case) and reconnects a disjoint band split (the
  disconnected case — where a band's IFB column would otherwise be collinear with the
  stranded constellation's clock), making the inter-frequency bias observable again
  (it then carries the broadcast-offset error, alongside the offset-based
  inter-system bias). A system whose satellites carry no offset to the hub keeps its
  own clock column, so a mixed epoch can collapse Galileo and leave BeiDou
  independent, or the reverse.
- Failing that, fall back to the (already observability-restricted) independent layout
  if the satellite count allows. No IFB column is created for a band stranded on its own
  constellation, so its inter-frequency bias folds into that constellation's clock and
  the reported inter-system bias reads out as `GGTO + (IFB difference)` — the two are no
  longer separable. Position, residuals and DOP are still finite and non-degenerate;
  only the bias decomposition is ambiguous. Else return `nothing`.

Because `band_ifb_layout` never creates an unobservable IFB column, every returned
layout is structurally sound — the degenerate disjoint-band case is removed by
construction, not caught after the fact. The satellite conditions above are structural
too, and thus necessary but not sufficient: a returned layout can still have degenerate
*geometry* (lines of sight that span too little), which no satellite count can see. That
is left to the checks `calc_pvt` makes on the solved geometry — the DOP's positive-definite
test and the velocity solve's own — rather than pre-screened.
"""
function decide_bias_layout(measurements)::Union{Nothing,BiasLayout}
    workspace = BiasLayoutWorkspace()
    decide_bias_layout!(workspace, measurements) ? bias_layout(workspace) : nothing
end

"""
    BiasLayoutWorkspace()

The storage [`decide_bias_layout!`](@ref) decides a layout into: the layout's own
vectors, which [`bias_layout`](@ref) hands out as a [`BiasLayout`](@ref), and the scratch
the decision works in. Every vector is emptied and refilled per epoch, so a workspace
reused across epochs allocates nothing once it has grown to the largest one.
"""
mutable struct BiasLayoutWorkspace
    # The decided layout.
    clock_bias_indices::Vector{Int}
    num_clock_biases::Int
    ifb_indices::Vector{Int}
    extra_bands::Vector{Symbol}
    reference_bands::Vector{Symbol}
    hub_system::Union{Nothing,SupportedTimeSystem}
    hub_rows::Vector{SatelliteMeasurement}
    # Scratch.
    bands::Vector{Symbol}
    effective_systems::Vector{SupportedTimeSystem}
    unique_effective_systems::Vector{SupportedTimeSystem}
    band_scratch::BandLayoutScratch{Symbol}
end

BiasLayoutWorkspace() = BiasLayoutWorkspace(Int[], 0, Int[], Symbol[], Symbol[], nothing,
    SatelliteMeasurement[], Symbol[], SupportedTimeSystem[], SupportedTimeSystem[],
    BandLayoutScratch{Symbol}())

"""
    bias_layout(workspace::BiasLayoutWorkspace) -> BiasLayout

The layout [`decide_bias_layout!`](@ref) last decided into `workspace`. Its vectors
**are** the workspace's, so they are overwritten by the next decision.
"""
bias_layout(workspace::BiasLayoutWorkspace) = BiasLayout((
    BiasColumns(workspace.clock_bias_indices, workspace.num_clock_biases,
        workspace.ifb_indices, length(workspace.extra_bands)),
    workspace.extra_bands,
    workspace.reference_bands,
    workspace.hub_system,
    workspace.hub_rows,
))

"""
    decide_bias_layout!(workspace::BiasLayoutWorkspace, measurements) -> Bool

[`decide_bias_layout`](@ref) into `workspace`: `true` and the layout stored there (read it
with [`bias_layout`](@ref)), or `false` when the constellation cannot be solved. The same
decision, allocating nothing once the workspace has grown to the epoch.
"""
function decide_bias_layout!(workspace::BiasLayoutWorkspace, measurements)
    num_sats = length(measurements)
    bands = resize!(workspace.bands, num_sats)
    for j in 1:num_sats
        bands[j] = measurements[j].band_id
    end
    # Distinct physical satellites, identified by `(time system, PRN)` — a PRN is only
    # unique within its GNSS. A satellite tracked on several bands appears once per band in
    # `measurements` but supplies one line of sight, so only distinct satellites constrain
    # the geometry and clock unknowns; its repeats constrain the inter-frequency biases.
    num_distinct_sats = count_distinct_satellites(measurements)
    # Both are necessary for a full-rank design (`H` has `3 + M + B` columns, and its rows
    # take only `num_distinct_sats` distinct values outside the IFB columns), neither is
    # sufficient: the geometry itself can still be degenerate, which `calc_pvt` screens
    # for once the design matrix exists.
    enough_satellites() =
        num_sats >= 3 + workspace.num_clock_biases + length(workspace.extra_bands) &&
        num_distinct_sats >= 3 + workspace.num_clock_biases

    # One row per collapsed system; a row carries its own `time_system`, so it is its
    # own key and no separate dictionary is needed.
    hub_rows = empty!(workspace.hub_rows)
    workspace.hub_system = nothing
    num_components = decide_bias_columns!(workspace, measurements, nothing)
    num_components == 1 && enough_satellites() && return true

    # Connected-but-scarce or disconnected: try collapsing every other system that
    # broadcasts an offset to a hub system onto that hub. The offset is one
    # constellation-wide value whichever of its satellites reports it, so the first
    # decoded copy per system converts all of that system's measurements. Hubs are
    # tried in a fixed order, so the choice is deterministic and GPS wins whenever it
    # can; BDT closes the list only for completeness — no signal broadcasts an offset
    # toward BDT today, so its loop finds nothing.
    for (hub_index, hub_system) in enumerate(CANDIDATE_HUB_SYSTEMS)
        any(measurement -> measurement.time_system === hub_system, measurements) || continue
        empty!(hub_rows)
        for measurement in measurements
            measurement.time_system === hub_system && continue
            is_collapsed(hub_rows, measurement.time_system) && continue
            measurement.time_offsets[hub_index].available || continue
            push!(hub_rows, measurement)
        end
        isempty(hub_rows) && continue
        decide_bias_columns!(workspace, measurements, hub_system)
        if enough_satellites()
            workspace.hub_system = hub_system
            return true
        end
    end

    # No collapse available. The independent layout is still observable (its IFBs are
    # component-restricted); use it if there are enough satellites, otherwise unsolvable.
    empty!(hub_rows)
    decide_bias_columns!(workspace, measurements, nothing)
    return enough_satellites()
end

# The clock and IFB columns of one candidate layout, into `workspace`: every system
# independent when `hub_system` is `nothing`, else the systems of `workspace.hub_rows`
# merged onto `hub_system`. Returns the number of coverage components.
function decide_bias_columns!(workspace::BiasLayoutWorkspace, measurements, hub_system)
    effective_systems = resize!(workspace.effective_systems, length(measurements))
    for (j, measurement) in enumerate(measurements)
        sys = measurement.time_system
        effective_systems[j] =
            !isnothing(hub_system) && is_collapsed(workspace.hub_rows, sys) ? hub_system :
            sys
    end
    unique_effective = unique_time_systems!(workspace.unique_effective_systems,
        effective_systems)
    clock_bias_indices = resize!(workspace.clock_bias_indices, length(measurements))
    # `something`: every system is in `unique_effective` by construction, and saying
    # so keeps the columns a `Vector{Int}` rather than widening them with `Nothing`.
    for (j, sys) in enumerate(effective_systems)
        clock_bias_indices[j] = something(time_system_index(unique_effective, sys))
    end
    workspace.num_clock_biases = length(unique_effective)
    # The clock column, not the time system itself, keys the coverage graph: it
    # separates the constellations identically (and in the same first-appearance
    # order) while being a concretely-typed `Int`, where the `Union`-typed
    # `time_system` field would branch on the system at every comparison inside
    # `band_ifb_layout!`.
    band_ifb_layout!(workspace.ifb_indices, workspace.extra_bands,
        workspace.reference_bands, workspace.band_scratch, clock_bias_indices,
        workspace.bands)
end

# Four identity-based helpers over the `time_system` field. Every
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

# Whether `system` is one of the systems collapsed onto the hub, i.e. whether
# `hub_rows` holds a representative row for it.
is_collapsed(hub_rows, system) = any(row -> row.time_system === system, hub_rows)

# Distinct physical satellites among `measurements`, identified by `(time system, PRN)`
# — a PRN is unique only within its GNSS, and a satellite tracked on several bands
# contributes several rows but one line of sight.
#
# Counted by scanning the rows already in hand, rather than by collecting the keys into
# a set: a `(SupportedTimeSystem, Int)` key is not a concrete type — the time system is
# a `Union` — so a vector of them heap-allocates one box per satellite, inside the solver,
# on every epoch. That is precisely the per-satellite allocation the flat row exists to
# remove. The scan is quadratic in the satellite count where the set would be linear,
# which for the dozens of rows an epoch holds is the cheaper of the two by a wide
# margin, and it allocates nothing at all.
function count_distinct_satellites(measurements)
    distinct = 0
    for j in eachindex(measurements)
        measurement = measurements[j]
        repeated = any(firstindex(measurements):(j-1)) do i
            measurements[i].prn == measurement.prn &&
                measurements[i].time_system === measurement.time_system
        end
        repeated || (distinct += 1)
    end
    distinct
end

"""
    time_scale_offset_to_gpst(time_system) -> Float64

Signed offset of a GNSS time system's *count* against GPS Time's for the same
instant — `get_tai_offset(GPST) − get_tai_offset(time_system)`, so **negative**
where the system's count reads lower (BDT: `19 − 33 = −14.0`), and `0.0` for a
system that counts alike. `calc_time_scale_offsets` adds the *difference of two
of these* to a transmit time, which is what puts the +14 s onto a BeiDou time
in a GPS-primary solve.

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
    calc_time_scale_offsets(measurements, primary_system) -> Vector{Float64}

Seconds to add to each satellite's transmit time to express it in
`primary_system`'s count, so that `calc_pseudo_ranges` may difference them.
`0.0` for every satellite of a system that counts alike — which is every GPS and
Galileo satellite, and all of them in a single-constellation epoch.

Relative to the *primary* system rather than to GPS Time, because the primary
system's count is also what dates the reported epoch: `reference_time` comes out
of the same differencing, and is combined with the primary system's week and
`system_start_epoch`. Anchoring on GPST instead leaves a BeiDou-primary solve
reporting a GPS-count time of week dated from the BDT epoch — a 14 s error in
`PVTSolution.time` traded for the one this removes from the pseudoranges.

This is applied to the times *handed to the differencing*, and nowhere else. The
satellite position must still be propagated from the untouched transmit time
against the message's own `t_0e`, which is on the broadcasting system's scale, so
shifting the time itself would move the offset into the ephemeris — about 55 km
of along-track error at BeiDou MEO velocities. `SatInfo.time` and
[`calc_steering_offset`](@ref) likewise keep the unconverted value.
"""
function calc_time_scale_offsets(measurements, primary_system)
    primary = time_scale_offset_to_gpst(primary_system)
    # Each row already carries its own anchor (`count_offset_to_gpst`), precomputed by
    # `collect_measurements`, so this reads a number rather than branching on the
    # `time_system` field per satellite.
    map(measurement -> primary - measurement.count_offset_to_gpst, measurements)
end

"""
    calc_hub_range_offsets(measurements, hub_rows, hub_system) -> Vector{Float64}

Per-satellite range offsets (metres) that carry a clock collapse into the measurements,
as decided by [`decide_bias_layout`](@ref): all-zero for the satellites of a system
that keeps its own clock unknown — all of them when `hub_rows` is empty —
and otherwise `−c · Δt_systems` for each satellite of a collapsed system, evaluated at
its own transmit time.

The broadcast offset is `Δt_systems = (that system's time) − (the hub system's)` (see
[`calc_steering_offset`](@ref): the GGTO for Galileo toward a GPS hub, the BGTO for
BeiDou toward either), so a transmit time becomes hub time by SUBTRACTING it; the
modeled range therefore carries `−c·Δt_systems`, and the solve yields
`inter_system_biases[sys] = −c·Δt_systems`. Which satellite of a system reported the
offset does not matter — it is one constellation-wide value — so
`decide_bias_layout` picks the first decoded copy per system (`hub_rows`) and it
converts all of that system's measurements.
"""
calc_hub_range_offsets(measurements, hub_rows, hub_system) = calc_hub_range_offsets!(
    Vector{Float64}(undef, length(measurements)), measurements, hub_rows, hub_system)

# `calc_hub_range_offsets` into `offsets`, resized to the satellite count.
function calc_hub_range_offsets!(offsets, measurements, hub_rows, hub_system)
    offsets = resize!(offsets, length(measurements))
    fill!(offsets, 0.0)
    isnothing(hub_system) && return offsets
    hub_index = hub_system_index(hub_system)
    for (j, measurement) in enumerate(measurements)
        # `hub_rows` holds at most one row per collapsed system (two in practice), so a
        # linear identity scan beats hashing a key.
        row_index =
            findfirst(row -> row.time_system === measurement.time_system, hub_rows)
        isnothing(row_index) && continue
        offsets[j] =
            -SPEED_OF_LIGHT * calc_steering_offset(
                hub_rows[row_index].time_offsets[hub_index], measurement.time)
    end
    offsets
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
        elevation, azimuth = _elevation_azimuth(enu_from_ecef, measurement.position)
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

The reusable scratch storage of [`calc_pvt!`](@ref): the flat measurement rows, the bias
layout, the pseudoranges and their corrections, the least-squares buffers, the design
and normal-equations matrices and the residuals. Every buffer grows to the largest epoch
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
    layout::BiasLayoutWorkspace
    unique_systems::Vector{SupportedTimeSystem}
    sat_positions::Vector{SVector{3,Float64}}
    pseudo_ranges::Vector{Float64}
    hub_offsets::Vector{Float64}
    atmospheric_delays::Vector{Float64}
    corrected_ranges::Vector{Float64}
    prev_ξ::Vector{Float64}
    residuals::Vector{Float64}
    rate_residuals::Vector{Float64}
    design_matrix::Matrix{Float64}
    normal_matrix::Matrix{Float64}
    least_squares::LMWorkspace
end

PVTWorkspace() = PVTWorkspace(
    SatelliteMeasurement[],
    BiasLayoutWorkspace(),
    SupportedTimeSystem[],
    SVector{3,Float64}[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Matrix{Float64}(undef, 0, 0),
    Matrix{Float64}(undef, 0, 0),
    LMWorkspace(),
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
condition fails, constellations whose messages carry a broadcast offset to another
tracked system — Galileo's GGTO (Galileo–GPS Time Offset), BeiDou's BGTO (BDT–GNSS
Time Offset, toward GPS or Galileo), or the GGTO GPS itself broadcasts on CNAV/CNAV-2
— have their clock bias collapsed onto that hub system using the broadcast offset,
which makes a 4-satellite mixed fix possible (see [`decide_bias_layout`](@ref) for the
hub order). Estimating an independent bias is preferred whenever the geometry allows
it, since it avoids the broadcast offset's own error.

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
  the satellites are tracked — with `Tracking` loaded,
  [`signal_groups`](@ref)`(track_state, decoders)` builds them from a whole `TrackState`.
  A pooled `Vector{SatelliteState}`, which is what this function took before 6.0, is
  refused with an error saying what to build instead.

  Each `(signal, PRN)` pair must appear at most once — a receiver produces one
  measurement per signal per satellite, and a duplicate would enter the least-squares
  solve twice. A `Dictionary`-backed group prevents this within itself; two groups that
  share a ranging signal can still collide, and are not checked for.

  **Order is significant**, as vector order was before it: the flat measurement order
  is group order × within-group order, and it fixes the primary system's tie-break
  ("most satellites, ties by first appearance"), [`band_ifb_layout`](@ref)'s
  reference-band tie-break, and the insertion order of `PVTSolution.sats`.
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
solve reusing the containers of a solution passed to it, without allocating. Returns `prev_pvt` if the epoch cannot be solved: too few healthy
satellites to solve the constellation (including the GGTO fallback and the
distinct-satellite condition — a satellite tracked on several bands supplies one line of
sight, so measurements alone are not enough), a geometry whose solved design matrix is
rank deficient (reported as a negative GDOP). None of these throw, so a receiver can pass
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
solution takes over its `sats`, `inter_system_biases` and `inter_frequency_biases`
containers, which are emptied and refilled with this epoch's values. `solution` itself
is therefore spent — its own fields still describe the epoch it was computed for, but
its containers now hold the new one — so use the returned solution from here on, and
copy anything of `solution` you want to keep before the call (a `Dictionary` made by
`map` over `solution.sats` shares its keys, and sees the change too). `prev_pvt` is only
read. Passing the same solution as both is allowed and is the usual receiver loop above:
every value the seed needs is read before anything is written.

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

# A copy of `source` in the containers of `solution`. `source` itself where the two
# share their containers (above all, where they are the same solution): there is
# nothing to copy, and copying a container onto itself would empty it.
function copy_into(solution::PVTSolution, source::PVTSolution)
    solution.sats === source.sats &&
        solution.inter_system_biases === source.inter_system_biases &&
        solution.inter_frequency_biases === source.inter_frequency_biases &&
        return source
    copy_container!(solution.sats, source.sats)
    copy_container!(solution.inter_system_biases, source.inter_system_biases)
    copy_container!(solution.inter_frequency_biases, source.inter_frequency_biases)
    PVTSolution(
        source.position,
        source.velocity,
        source.course_over_ground,
        source.time_correction,
        source.time,
        source.relative_clock_drift,
        source.dop,
        solution.sats,
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

# `empty!(::Dictionary)` hands the dictionary a fresh, minimal hash table, so refilling
# it reallocates that table (and regrows it) on every epoch. This empties it in place
# instead, keeping every buffer's capacity: an all-zero slot table with no hashes,
# keys, values or holes is exactly the empty state Dictionaries' own `empty!` builds,
# only with the slot table at its current size — which is a power of two, as the
# hashing requires, because only Dictionaries itself ever sized it.
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

    # Gauss-Newton converges within its seed's basin, so a spurious far-away root
    # would re-seed itself epoch after epoch. Such roots flag themselves with an
    # impossible geometry (negative GDOP, exploded PDOP): solve from a cold seed
    # instead. Only the seed is affected — an unsolvable epoch still returns
    # `prev_pvt` unchanged.
    distrusted =
        !isnothing(prev_pvt.dop) &&
        (prev_pvt.dop.GDOP < 0 || prev_pvt.dop.PDOP > MAX_TRUSTED_WARM_START_PDOP)

    num_sats = length(measurements)

    # Every satellite that was not decoded-complete-and-healthy has already been
    # dropped by the collection pass, and each surviving row carries the GNSSSignals
    # keys that drive the solution: `time_system` (`GPST()`/`GST()`/`BDT()`) groups the
    # receiver clock bias — one per time system, ordered by first appearance —
    # and `band_id` (`:L1`, `:L5`, …) groups the receiver inter-frequency bias — one
    # per band beyond a per-coverage-component reference. (`signal_id`, `:GPSL1CA` …, is
    # the per-signal `sats` identity used below, not a grouping key.)
    #
    # Solvability is decided here: `decide_bias_layout!` keeps only observable IFBs and
    # falls back to a hub collapse when the geometry is disconnected or
    # under-determined. A degenerate geometry, which no count can see, is caught after
    # the solve by the DOP.
    decide_bias_layout!(workspace.layout, measurements) || return (false, prev_pvt)
    (; bias_columns, extra_bands, reference_bands, hub_system, hub_rows) =
        bias_layout(workspace.layout)
    (; clock_bias_indices, num_clock_biases) = bias_columns
    num_params = num_lsq_params(bias_columns)

    # The propagated ephemerides, already evaluated per satellite by the collection pass.
    sat_positions = resize!(workspace.sat_positions, num_sats)
    for j in 1:num_sats
        sat_positions[j] = measurements[j].position
    end

    # Primary system — its clock bias, reference time, week and start epoch
    # define the reported time. A collapse is anchored on its hub system, so the
    # hub must be primary there; otherwise pick the system with the most
    # satellites (best-conditioned reported time), breaking ties by first
    # appearance.
    unique_systems = unique_time_systems!(workspace.unique_systems,
        measurement.time_system for measurement in measurements)
    primary_system = if isnothing(hub_system)
        most_satellites, most_index = 0, 0
        for (i, sys) in enumerate(unique_systems)
            n = count(measurement -> measurement.time_system === sys, measurements)
            n > most_satellites && ((most_satellites, most_index) = (n, i))
        end
        unique_systems[most_index]
    else
        hub_system
    end
    # `something`, not a bare `findfirst`: the primary system is either one of the
    # systems present or the hub, which `decide_bias_layout!` only picks when it is
    # present — so this never fails, and saying so keeps the index an `Int` rather than
    # a `Union{Int,Nothing}` that would widen the clock column and the DOP call below.
    primary_index = something(
        findfirst(measurement -> measurement.time_system === primary_system, measurements))
    primary_clock_index = clock_bias_indices[primary_index]

    # The common reference cancels out of the reported time (the primary clock
    # bias absorbs it), so any latest-transmit-time reference works — but the times
    # must first be put on one count (`calc_time_scale_offsets`, added here row by row).
    # A BeiDou second-of-week reads 14 s below a GPS time of week for the same instant,
    # so differencing them raw hands every BeiDou measurement 4.2e9 m of structural
    # offset: not merely a biased BeiDou clock column, but a parameter nine orders of
    # magnitude above the others in the normal equations. Zero for GPS and Galileo, and
    # for any single-constellation epoch.
    pseudo_ranges = resize!(workspace.pseudo_ranges, num_sats)
    primary_count_offset = time_scale_offset_to_gpst(primary_system)
    for (j, measurement) in enumerate(measurements)
        pseudo_ranges[j] =
            measurement.time + (primary_count_offset - measurement.count_offset_to_gpst)
    end
    _, reference_time = calc_pseudo_ranges!(pseudo_ranges)
    # The known per-satellite broadcast steering offset (zero unless that system was
    # collapsed onto the hub), as a measurement correction like the atmospheric delays
    # below. Kept for the inter-system-bias readout at the end, which reports this
    # term alone: the time-count difference already removed above is a convention, not
    # a bias, and belongs in neither the readout nor the solve.
    hub_offsets =
        calc_hub_range_offsets!(workspace.hub_offsets, measurements, hub_rows, hub_system)
    pseudo_ranges .-= hub_offsets

    # The primary system's week and start epoch date the epoch absolutely: the day of
    # year feeds the tropospheric mapping's seasonal term here, and week/start epoch
    # date the reported time after the solve.
    primary_measurement = measurements[primary_index]
    week = primary_measurement.week
    start_time = primary_measurement.system_start_epoch
    doy = day_of_year(primary_measurement.system_start_time, week, reference_time)

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias. A
    # distrusted previous solution seeds nothing: the state stays all zeros, the cold
    # start.
    prev_abs_bias(sys) =
        sys === prev_pvt.reference_system ? prev_pvt.time_correction :
        prev_pvt.time_correction + get(prev_pvt.inter_system_biases, sys, 0.0m)
    prev_ξ = fill!(resize!(workspace.prev_ξ, num_params), 0.0)
    if !distrusted
        prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
        for j in 1:num_sats
            prev_ξ[3+clock_bias_indices[j]] =
                ustrip(m, prev_abs_bias(measurements[j].time_system))
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
                haskey(prev_ifbs, band) && prev_ifbs[band].reference == reference_bands[i] ?
                ustrip(m, prev_ifbs[band].value) : 0.0
        end
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
        resize!(workspace.corrected_ranges, num_sats) .= pseudo_ranges .- delays
    end
    least_squares = workspace.least_squares

    # `ξ` is the least-squares workspace's own parameter vector, so it stays valid until
    # the next solve on that workspace — which is after everything below.
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, residuals_uncorrected = user_position!(least_squares,
            workspace.residuals, sat_positions, pseudo_ranges, bias_columns, prev_ξ)
        if correct_atmosphere
            user_position!(least_squares, workspace.residuals, sat_positions,
                corrected_ranges(ξ_uncorrected), bias_columns, ξ_uncorrected)
        else
            (ξ_uncorrected, residuals_uncorrected)
        end
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction. With no
        # active correction, solve directly on the raw pseudoranges.
        user_position!(least_squares, workspace.residuals, sat_positions,
            correct_atmosphere ? corrected_ranges(prev_ξ) : pseudo_ranges, bias_columns,
            prev_ξ)
    end
    workspace.design_matrix = grown(workspace.design_matrix, num_sats, num_params)
    H = calc_H!(view(workspace.design_matrix, 1:num_sats, 1:num_params), sat_positions, ξ,
        bias_columns)
    position = ECEF(ξ[1], ξ[2], ξ[3])

    # Check the geometry at the converged position — the DOP is reported to the caller, and
    # a rank deficiency here (negative GDOP) means `ξ` is meaningless, so nothing should be
    # derived from it. The satellite conditions in `decide_bias_layout!` already reject the
    # count-shaped degeneracies before the solve; this catches what a count cannot see.
    #
    # This check must stay ahead of the velocity solve: it is what makes that solve's
    # normal-equations matrix positive definite (see `calc_user_velocity_and_clock_drift`),
    # and with the order reversed a degenerate geometry throws `SingularException` there.
    workspace.normal_matrix = grown(workspace.normal_matrix, num_params, num_params)
    dop = calc_DOP!(view(workspace.normal_matrix, 1:num_params, 1:num_params), H, position,
        primary_clock_index)

    dop.GDOP < 0 && return (false, prev_pvt)

    user_velocity_and_clock_drift, rate_residuals =
        calc_user_velocity_and_clock_drift!(workspace.rate_residuals, measurements, H)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    time_correction = ξ[3+primary_clock_index]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEED_OF_LIGHT

    # Everything is read from `prev_pvt` by now, so the containers of `solution` —
    # possibly `prev_pvt`'s own — are written from here on.
    #
    # Per-satellite `sats` key: (signal id, PRN) — signal-level (not time system), so a
    # satellite tracked on two signals of one constellation stays distinct; the
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
    sats = empty_keeping_capacity!(solution.sats)
    for (j, measurement) in enumerate(measurements)
        key = (measurement.signal_id, measurement.prn)
        set!(sats, key, SatInfo(measurement.position, measurement.time, residuals[j] * m,
            rate_residuals[j] * (m / s)))
    end

    # Inter-system biases relative to the reference (primary) system's clock, in
    # meters. The reference is omitted (its bias is `time_correction`); for a
    # collapsed system this is the broadcast offset −c·Δt_systems, read from the
    # system's first satellite (the per-satellite offsets differ only by the offset
    # polynomial's drift term over the transmit-time spread — sub-millimetre).
    inter_system_biases = empty!(solution.inter_system_biases)
    for sys in unique_systems
        sys === primary_system && continue
        j = something(findfirst(measurement -> measurement.time_system === sys, measurements))
        inter_system_biases[sys] =
            (ξ[3+clock_bias_indices[j]] + hub_offsets[j] - time_correction) * m
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
        time_correction * m,
        # Assumes `start_time.fraction == 0` (true for GPS/Galileo: integer-second
        # origins).
        TAITime(
            week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) +
            start_time.second,
            corrected_reference_time - floor(Int, corrected_reference_time),
        ),
        user_velocity_and_clock_drift[4] / SPEED_OF_LIGHT,
        dop,
        sats,
        primary_system,
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
table is consulted. This used to go through `get_system_start_time` (a UTC
label) and AstroTime's leap-aware conversion; GNSSSignals 4.1 states the TAI
labels itself, which is also what makes the value safe to derive at
precompile time.
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
    get_LLA(pvt::PVTSolution) -> LLA

Convert the ECEF position in the PVT solution to geodetic coordinates
(latitude, longitude, altitude) using the WGS84 ellipsoid.
"""
function get_LLA(pvt::PVTSolution)
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
