module PositionVelocityTime
using CoordinateTransformations,
    DocStringExtensions,
    Geodesy,
    GNSSDecoder,
    GNSSSignals,
    LinearAlgebra,
    AstroTime,
    LsqFit,
    StaticArrays,
    Unitful,
    Dates

using Unitful: s, Hz, m, °, dBHz, ustrip, uconvert
using Dictionaries: Dictionary

export calc_pvt,
    PVTSolution,
    SatInfo,
    InterFrequencyBias,
    SatelliteState,
    get_LLA,
    get_sat_info,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    get_sat_enu

const SPEEDOFLIGHT = 299792458.0

# The GPS CNAV family (CNAV on L5/L2C, CNAV-2 on L1C): both broadcast the full week
# number (no 1024-week rollover) and a quasi-Keplerian ephemeris (A_REF + ΔA,
# Ω̇_REF + ΔΩ̇, …) rather than LNAV's directly-broadcast Keplerian elements — so
# `orbital_elements` and the full-week `get_week` dispatch on this union.
const GPSModernNavData = Union{GNSSDecoder.GPSCNAVData,GNSSDecoder.GPSL1C_DData}

"""
    BiasColumns

Per-satellite assignment of the estimated bias columns of the least-squares design
matrix, shared by [`calc_ρ_hat!`](@ref), [`calc_H!`](@ref) and [`user_position`](@ref).
The state vector is `[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifb]` with
[`num_lsq_params`](@ref)`(bias_columns)` entries. The two column kinds have different
physical sources: a clock column is the receiver clock for one GNSS time system (the
spacing *between* systems is a system/space-segment effect — the GNSS time offset /
GGTO), whereas an inter-frequency-bias column is the receiver's per-band RF-chain delay.
Known per-satellite corrections (atmosphere, the GGTO time-system offset) are applied to
the pseudoranges in [`calc_pvt`](@ref), not carried here.

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
    BiasLayout

How one epoch's estimated biases are laid out, as decided by
[`decide_bias_layout`](@ref): the design-matrix columns, the bands its
inter-frequency-bias columns belong to, and the source of the range offsets a
GGTO-collapsed layout needs.

# Fields
- `bias_columns::BiasColumns`: the per-satellite column assignment (see
  [`BiasColumns`](@ref)).
- `extra_bands::Vector{Symbol}`: band of each inter-frequency-bias column, so
  `extra_bands[i]` belongs to column `i` of the IFB block.
- `reference_bands::Vector{Symbol}`: per IFB column, the reference band of its coverage
  component — the anchor that column's bias is measured against (see
  [`band_ifb_layout`](@ref)).
- `ggto_decoder::Union{Nothing,GNSSDecoder.GNSSDecoderState}`: the Galileo decoder whose
  broadcast GGTO converts the collapsed measurements (see
  [`calc_ggto_range_offsets`](@ref)), or `nothing` for a layout that estimates every clock
  bias independently. Declared as a union so one concrete `BiasLayout` covers both.
"""
struct BiasLayout
    bias_columns::BiasColumns
    extra_bands::Vector{Symbol}
    reference_bands::Vector{Symbol}
    ggto_decoder::Union{Nothing,GNSSDecoder.GNSSDecoderState}
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
- `cn0::Union{Nothing,typeof(1.0dBHz)}`: Carrier-to-noise-density ratio of this
  measurement, e.g. `42.0dBHz`, matching `Tracking.estimate_cn0` (default: `nothing`,
  "not reported"). Supplying it lets [`calc_pvt`](@ref) weight this satellite by its
  measurement uncertainty instead of treating every satellite as equally precise —
  see [`pseudorange_variance`](@ref). A non-positive level is treated as not
  reported, as is a non-finite one.
- `pseudorange_variance::Union{Nothing,typeof(1.0m^2)}`: Optional explicit a-priori
  variance of this pseudorange, e.g. `(3.0m)^2`, for a caller that owns its own error
  model. When given it replaces the modeled variance entirely — C/N₀ and elevation
  are then not consulted for this satellite's position weight (default: `nothing`).

# Constructors
    SatelliteState(; decoder, system, code_phase, carrier_doppler, carrier_phase=0.0,
                   cn0=nothing, pseudorange_variance=nothing)
    SatelliteState(decoder, system, sat_state; cn0=…, pseudorange_variance=nothing)

The second constructor extracts code phase, carrier Doppler, carrier phase and C/N₀
from a `Tracking` satellite state (`Tracking.TrackedSat`). It is provided by a package
extension that is loaded automatically once `Tracking` is available, so `Tracking` is
only a weak dependency of this package.
"""
@kwdef struct SatelliteState{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    system::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
    cn0::Union{Nothing,typeof(1.0dBHz)} = nothing
    pseudorange_variance::Union{Nothing,typeof(1.0m^2)} = nothing
end

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
    FormalAccuracy

Formal (1σ) accuracy of a [`PVTSolution`](@ref) in metres: the square-rooted diagonal
of the least-squares parameter covariance `(HᵀWH)⁻¹`, with the position block rotated
into the local East-North-Up frame — see [`calc_formal_accuracy`](@ref).

This is the *weighted* counterpart of [`DOP`](@ref) and deliberately a separate
output. DOP is by definition a purely geometric quantity, `(HᵀH)⁻¹`, and stays that
way whether or not the epoch was weighted, so it keeps meaning what every receiver
and textbook means by it; the accuracy here folds in the a-priori measurement
uncertainty (see [`pseudorange_variance`](@ref)) and so carries units. When no
satellite reports any uncertainty the weights are uniform and this reduces to the
familiar DOP × nominal-UERE product.

It is *formal*: it propagates the assumed measurement variances through the geometry
and says nothing about errors the model does not describe (multipath beyond the model,
an undetected faulty satellite, broadcast-ephemeris error).

# Fields
- `horizontal::typeof(1.0m)`: Horizontal 1σ (East² + North²), i.e. the classic HDRMS.
- `vertical::typeof(1.0m)`: Vertical (Up) 1σ.
- `position::typeof(1.0m)`: 3D position 1σ (trace of the position block; frame
  independent, so `≈ hypot(horizontal, vertical)`).
- `time::typeof(1.0m)`: 1σ of the reported receiver clock bias, as a range (metres);
  divide by the speed of light for seconds. This is the clock of `reference_system`,
  matching `DOP.TDOP`.
"""
struct FormalAccuracy
    horizontal::typeof(1.0m)
    vertical::typeof(1.0m)
    position::typeof(1.0m)
    time::typeof(1.0m)
end

"""
    SatInfo

Per-satellite information attached to a [`PVTSolution`](@ref) (one entry per
satellite used in the fix).

# Fields
- `position::ECEF`: Satellite ECEF position at transmit time (metres).
- `time::Float64`: Satellite transmit time (system time of week, seconds).
- `residual::typeof(1.0m)`: Post-fit least-squares pseudorange residual (metres) — the
  modeled minus the (atmosphere-corrected) measured pseudorange. A per-satellite
  fit-quality / outlier indicator. This is the **raw** residual in metres, not the
  weight-normalised one, so it keeps the same meaning whether or not the epoch was
  weighted; divide it by `pseudorange_sigma` for the normalised residual that fault
  detection (RAIM) works with.
- `rate_residual::typeof(1.0m/s)`: Post-fit least-squares range-rate residual (metres per
  second) — the modeled minus the measured range rate of the carrier-Doppler velocity and
  clock-drift solve. The rate-domain counterpart of `residual`: it flags a satellite whose
  Doppler disagrees with the velocity fix (cycle slips, dynamics) independently of its
  pseudorange.
- `pseudorange_sigma::typeof(1.0m)`: A-priori standard deviation (metres) this
  satellite's pseudorange was weighted by — its [`pseudorange_variance`](@ref),
  square-rooted. When no satellite of the epoch reported any uncertainty the solve is
  unweighted and this is the nominal UERE assumed for the reported accuracy, the same
  value for every satellite.
"""
struct SatInfo
    position::ECEF
    time::Float64
    residual::typeof(1.0m)
    rate_residual::typeof(1.0m/s)
    pseudorange_sigma::typeof(1.0m)
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
- `position::ECEF`: User position in ECEF coordinates (meters)
- `velocity::ECEF`: User velocity in ECEF coordinates (m/s)
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
- `time::Union{TAIEpoch{Float64}, Nothing}`: Estimated time as a TAI epoch
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values — purely geometric, and
  unaffected by measurement weighting (see [`FormalAccuracy`](@ref))
- `accuracy::Union{FormalAccuracy, Nothing}`: Formal 1σ accuracy of this solution in
  metres (horizontal, vertical, 3D and clock), from the least-squares covariance with
  the a-priori measurement variances folded in. See [`FormalAccuracy`](@ref)
- `sats::Dictionary{Tuple{Symbol, Int}, SatInfo}`: Maps `(signal, PRN)` to satellite
  info — position, transmit time, the post-fit pseudorange and range-rate residuals, and
  the a-priori σ the pseudorange was weighted by (see [`SatInfo`](@ref)). The signal tag
  (e.g. `:GPSL1CA`,
  `:GalileoE1B`; see `get_signal_id`) keeps the
  same PRN apart both across constellations (GPS PRN 5 vs Galileo E05) and across
  signals of one constellation (a satellite tracked on GPS L1 C/A and L5 yields two
  entries sharing a PRN). Receiver-clock grouping is by time system, not signal —
  see `reference_system`.
- `reference_system::Union{GNSSSignals.TimeSystem, Nothing}`: GNSS time system (e.g.
  `GNSSSignals.GPST()`, `GST()`) that `time` and `time_correction` are referenced to.
- `inter_system_biases::Dict{GNSSSignals.TimeSystem, typeof(1.0m)}`: For each GNSS time system
  other than `reference_system`, the offset of that system's time scale relative to the reference
  system's (meters) — the inter-system bias. This is a **system / space-segment** effect,
  not a receiver one (this receiver has no inter-system hardware bias): it is the GNSS
  system-time offset, so the Galileo entry equals `−c · Δt_systems`, `Δt_systems = GST −
  GPST` (the GGTO). It is **estimated directly from the geometry whenever observable**
  (no broadcast error); the broadcast GGTO is used to derive it only as a fallback (the
  GGTO-aided collapse), as that broadcast value may be erroneous. Reference-independent
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
"""
@kwdef struct PVTSolution
    position::ECEF = ECEF(0, 0, 0)
    velocity::ECEF = ECEF(0, 0, 0)
    course_over_ground::typeof(1.0°) = 0.0°
    time_correction::typeof(1.0m) = 0.0m
    time::Union{TAIEpoch{Float64},Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    accuracy::Union{FormalAccuracy,Nothing} = nothing
    sats::Dictionary{Tuple{Symbol,Int},SatInfo} = Dictionary{Tuple{Symbol,Int},SatInfo}()
    reference_system::Union{GNSSSignals.TimeSystem,Nothing} = nothing
    inter_system_biases::Dict{GNSSSignals.TimeSystem,typeof(1.0m)} =
        Dict{GNSSSignals.TimeSystem,typeof(1.0m)}()
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
    band_ifb_layout(systems, bands) -> (ifb_indices, extra_bands, reference_bands, num_components)

Lay out the receiver inter-frequency biases from the (constellation × band) coverage
graph. Two bands share a *coverage component* iff some constellation is tracked on
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
function band_ifb_layout(systems, bands)
    unique_bands = unique(bands)
    # Union-find over bands: union the bands a single constellation is tracked on.
    parent = Dict(b => b for b in unique_bands)
    root(b) = parent[b] == b ? b : (parent[b] = root(parent[b]))
    function link_bands!(a, c)
        ra, rc = root(a), root(c)
        ra == rc || (parent[ra] = rc)
    end
    for sys in unique(systems)
        sys_bands = unique(bands[i] for i in eachindex(systems) if systems[i] == sys)
        for k in 2:length(sys_bands)
            link_bands!(sys_bands[1], sys_bands[k])
        end
    end
    band_count = Dict{eltype(unique_bands),Int}()
    for b in bands
        band_count[b] = get(band_count, b, 0) + 1
    end
    # Reference band per component = most-populated in the component (ties: first seen).
    reference_of = Dict{eltype(unique_bands),eltype(unique_bands)}()
    for b in unique_bands
        r = root(b)
        if !haskey(reference_of, r) || band_count[b] > band_count[reference_of[r]]
            reference_of[r] = b
        end
    end
    extra_bands = filter(b -> reference_of[root(b)] != b, unique_bands)
    reference_bands = [reference_of[root(b)] for b in extra_bands]
    band_column = Dict(b => i for (i, b) in enumerate(extra_bands))
    ifb_indices = [get(band_column, b, 0) for b in bands]
    return ifb_indices, extra_bands, reference_bands, length(reference_of)
end

"""
    decide_bias_layout(states, systems, bands) -> Union{BiasLayout,Nothing}

Decide the full least-squares bias layout — one clock column per GNSS time system plus
the per-band inter-frequency-bias columns from [`band_ifb_layout`](@ref) — and return it
as a [`BiasLayout`](@ref), whose `ggto_decoder` also carries how a collapsed layout's
measurements are converted. Returns `nothing` when the constellation cannot be solved.

Every input is a classification of the epoch, never a transmit time: the counts, the
coverage graph and `ggto_available` are all the decision needs. The offsets its
`ggto_decoder` enables are built afterwards, once [`calc_pvt`](@ref) has the transmit
times.

The decision is observability-driven, not merely count-driven:

- When the coverage graph is connected and the satellites suffice for the unknowns
  (`n ≥ 3 + num_systems + num_ifb` measurements, of which `3 + num_systems` must come
  from *distinct* satellites — extra bands of an already-tracked satellite add
  inter-frequency-bias information, not geometry), estimate everything independently: the
  inter-system offset and the receiver inter-frequency biases are observed directly from
  the geometry, so neither inherits the broadcast-GGTO error (the satellite group delays
  are already removed per satellite upstream, so the per-band column carries the receiver
  chain).
- Otherwise merge the Galileo clock onto GPS via the broadcast GGTO when a Galileo
  satellite carries it. This removes a clock unknown (the scarce-satellite case) and
  reconnects a disjoint GPS/Galileo band split (the disconnected case — where a band's
  IFB column would otherwise be collinear with the stranded constellation's clock),
  making the inter-frequency bias observable again (it then carries the broadcast-GGTO
  error, alongside the GGTO-based inter-system bias).
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
function decide_bias_layout(states, systems, bands)
    num_sats = length(states)
    # Distinct physical satellites, identified by `(time system, PRN)` — a PRN is only
    # unique within its GNSS. A satellite tracked on several bands appears once per band in
    # `states` but supplies one line of sight, so only distinct satellites constrain the
    # geometry and clock unknowns; its repeats constrain the inter-frequency biases.
    num_distinct_sats = length(unique(zip(systems, (state.decoder.prn for state in states))))
    # Both are necessary for a full-rank design (`H` has `3 + M + B` columns, and its rows
    # take only `num_distinct_sats` distinct values outside the IFB columns), neither is
    # sufficient: the geometry itself can still be degenerate, which `calc_pvt` screens
    # for once the design matrix exists.
    enough_satellites(layout) =
        num_sats >= 3 + layout.num_clock_biases + length(layout.extra_bands) &&
        num_distinct_sats >= 3 + layout.num_clock_biases

    function bias_layout_for(effective_systems)
        unique_effective = unique(effective_systems)
        index = Dict(sys => i for (i, sys) in enumerate(unique_effective))
        clock_bias_indices = [index[sys] for sys in effective_systems]
        ifb_indices, extra_bands, reference_bands, num_components =
            band_ifb_layout(effective_systems, bands)
        (; clock_bias_indices, num_clock_biases = length(unique_effective), ifb_indices,
            extra_bands, reference_bands, num_components)
    end

    as_bias_layout(layout, ggto_decoder) = BiasLayout(
        BiasColumns(layout.clock_bias_indices, layout.num_clock_biases,
            layout.ifb_indices, length(layout.extra_bands)),
        layout.extra_bands,
        layout.reference_bands,
        ggto_decoder,
    )

    independent_layout = bias_layout_for(systems)
    if independent_layout.num_components == 1 && enough_satellites(independent_layout)
        return as_bias_layout(independent_layout, nothing)
    end

    # Connected-but-scarce or disconnected: try the GGTO collapse (merge Galileo onto
    # GPS). The GGTO is the same constellation-wide offset whichever Galileo satellite
    # reports it, so one decoded copy converts every Galileo measurement.
    ggto_idx = findfirst(
        j -> systems[j] == GST() && ggto_available(states[j].decoder), 1:num_sats)
    if (GPST() in systems) && !isnothing(ggto_idx)
        merged_layout = bias_layout_for(map(sys -> sys == GST() ? GPST() : sys, systems))
        if enough_satellites(merged_layout)
            return as_bias_layout(merged_layout, states[ggto_idx].decoder)
        end
    end

    # No collapse available. The independent layout is still observable (its IFBs are
    # component-restricted); use it if there are enough satellites, otherwise unsolvable.
    enough_satellites(independent_layout) ?
    as_bias_layout(independent_layout, nothing) : nothing
end

"""
    calc_ggto_range_offsets(ggto_decoder, systems, times) -> Vector{Float64}

Per-satellite range offsets (metres) that carry a GGTO collapse into the measurements, as
decided by [`decide_bias_layout`](@ref): all-zero when `ggto_decoder === nothing` (every
time system keeps its own clock unknown, so no measurement needs converting), and
otherwise `−c · Δt_systems` for each Galileo satellite, evaluated at its own transmit
time.

Per the OS SIS ICD the GGTO is `Δt_systems = GST − GPST` (see
[`calc_ggto_offset`](@ref)), so a Galileo transmit time becomes GPS time by SUBTRACTING
it; the modeled range therefore carries `−c·GGTO`, and the solve yields
`inter_system_biases[GST()] = −c·(GST − GPST)`. Which Galileo satellite reports the GGTO
does not matter — it is one constellation-wide offset — so `decide_bias_layout` picks the
first decoded copy and it converts every Galileo measurement.
"""
function calc_ggto_range_offsets(ggto_decoder, systems, times)
    offsets = zeros(length(systems))
    isnothing(ggto_decoder) && return offsets
    for j in eachindex(systems)
        systems[j] == GST() || continue
        offsets[j] = -SPEEDOFLIGHT * calc_ggto_offset(ggto_decoder, times[j])
    end
    offsets
end

"""
    predict_atmospheric_delays(ξ, states, sat_positions, correction,
                               reference_time, enable_tropospheric_correction) -> Vector{Float64}

Per-satellite ionospheric + tropospheric delay (metres), to be subtracted from the
pseudoranges. The user position is the first three elements of the least-squares
state vector `ξ = [x, y, z, tc₁, …]` (ECEF, metres); the remaining clock-bias
components are ignored. `correction` is the constellation-wide ionospheric model
from [`select_ionospheric_correction`](@ref) (`nothing` skips the ionosphere); the
troposphere uses the blind Saastamoinen model (see [`tropospheric_delay`](@ref))
unless `enable_tropospheric_correction` is `false`.

A single corrected solve is enough: the delays depend on position only through the
satellite elevation/azimuth (and, for the troposphere, the user height), and
∂delay/∂position is negligible over the metre-level position uncertainty (a 15 m
shift moves the elevation by ~1e-5°), so delays predicted at a nearby position are
accurate to well under a millimetre — no iterate-to-convergence needed. The user
geodetic coordinates and the ENU transform depend only on `ξ`, so they are
computed once and reused across satellites.
"""
function predict_atmospheric_delays(
    ξ,
    states,
    sat_positions,
    correction,
    reference_time,
    enable_tropospheric_correction,
)
    user_pos = ECEF(ξ[1], ξ[2], ξ[3])
    user_lla = LLAfromECEF(wgs84)(user_pos)
    enu_from_ecef = ENUfromECEF(user_pos, wgs84)
    map(states, sat_positions) do state, sat_pos
        elevation, azimuth = _elevation_azimuth(enu_from_ecef, sat_pos)
        iono = ionospheric_delay(
            correction,
            state.system,
            elevation,
            azimuth,
            user_lla,
            reference_time,
        )
        tropo =
            enable_tropospheric_correction ? tropospheric_delay(elevation, user_lla) : 0.0
        iono + tropo
    end
end

"""
    calc_pvt(states::AbstractVector{<:SatelliteState},
             prev_pvt::PVTSolution = PVTSolution();
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
condition fails but both GPS and Galileo
are present and the Galileo message carries the GGTO (Galileo–GPS Time Offset),
the Galileo clock bias is collapsed onto GPS using the broadcast offset, which
makes a 4-satellite GPS+Galileo fix possible. Estimating an independent bias is
preferred whenever the geometry allows it, since it avoids the GGTO's own
broadcast error.

Unless disabled via `enable_ionospheric_correction`, the ionospheric delay is
corrected automatically using only the coefficients decoded from the navigation
messages. A single model is chosen for the whole solve and applied to every
satellite: NTCM-G if Galileo Effective Ionisation Level coefficients have been
decoded (the more accurate model), otherwise the Klobuchar model if GPS L1 α/β
have been decoded, otherwise no correction. See
[`select_ionospheric_correction`](@ref) and [`ionospheric_delay`](@ref).

Unless disabled via `enable_tropospheric_correction`, the tropospheric delay is
corrected with the blind Saastamoinen model (no broadcast coefficients needed).
See [`tropospheric_delay`](@ref).

Satellites are weighted by their measurement uncertainty rather than treated as equally
precise, as soon as any of them reports a `cn0` (or an explicit `pseudorange_variance`)
on its [`SatelliteState`](@ref): the position solve then minimises `Σ (ρ̂ⱼ − ρⱼ)²/σ²ⱼ`
with σ from the C/N₀- and elevation-dependent model in [`pseudorange_variance`](@ref),
and the velocity solve is weighted by the Doppler's own
[`range_rate_variance`](@ref). This matters because code noise is a strong function of
C/N₀ — an order of magnitude between a 30 dBHz and a 45 dBHz satellite — so a marginal
satellite can otherwise drag a fix that the strong ones would have pinned down. Reported
per satellite as `SatInfo.pseudorange_sigma`, and summarised as the solution's formal
[`FormalAccuracy`](@ref). When no satellite reports either, the solve is the ordinary
least-squares one it has always been.

# Arguments
- `states`: Vector of [`SatelliteState`](@ref) for observed satellites. Each
  (signal, PRN) pair must appear at most once — a receiver produces one
  measurement per signal per satellite, and a duplicate would enter the
  least-squares solve twice.
- `prev_pvt`: Previous PVT solution used as initial guess (default: origin)

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
A [`PVTSolution`](@ref) containing position, velocity, time, DOP values, formal accuracy,
and satellite information. Returns `prev_pvt` if the epoch cannot be solved: too few
healthy satellites to solve the constellation (including the GGTO fallback and the
distinct-satellite condition — a satellite tracked on several bands supplies one line of
sight, so measurements alone are not enough), a geometry whose solved design matrix is
rank deficient (reported as a negative GDOP). None of these throw, so a receiver can pass
whatever it currently tracks each epoch and carry `prev_pvt` forward.
"""
function calc_pvt(
    states::AbstractVector{<:SatelliteState},
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    # Keep a satellite only if its full nav-data set is decoded and it reports healthy.
    # Checking completeness first guarantees the health bit has been decoded.
    healthy_indices = findall(
        x -> is_decoding_completed_for_positioning(x.decoder) && is_sat_healthy(x.decoder),
        states,
    )
    healthy_states = view(states, healthy_indices)
    num_sats = length(healthy_states)

    # Classify each satellite by the GNSSSignals keys that drive the solution.
    # `get_time_system` (`GPST()`/`GST()`) groups the receiver clock bias — one per time
    # system, ordered by first appearance; `get_band_id` (`:L1`, `:L5`, …) groups the
    # receiver inter-frequency bias — one per band beyond a per-coverage-component
    # reference. (`get_signal_id`, `:GPSL1CA` …, is the per-signal `sats` identity used
    # below, not a grouping key.) `decide_bias_layout` then decides the full bias layout:
    # it keeps only observable IFBs and falls back to a GGTO collapse when the geometry
    # is disconnected or under-determined.
    systems = map(state -> get_time_system(state.system), healthy_states)
    bands = map(state -> get_band_id(state.system), healthy_states)

    # Solvability is decided here. A degenerate geometry, which no count can see, is caught
    # after the solve by the DOP.
    bias_layout = decide_bias_layout(healthy_states, systems, bands)
    isnothing(bias_layout) && return prev_pvt
    (; bias_columns, extra_bands, reference_bands, ggto_decoder) = bias_layout
    (; clock_bias_indices, num_clock_biases) = bias_columns

    times = map(calc_corrected_time, healthy_states)

    # Propagating the ephemerides.
    sat_positions_and_velocities = map(
        (state, time) -> calc_satellite_position_and_velocity(state.decoder, time),
        healthy_states,
        times,
    )
    sat_positions = map(get_sat_position, sat_positions_and_velocities)
    # Built once here and reused below. `stack` collects the SVector
    # columns into a single 3×N `Matrix{Float64}` in one pass.
    sat_positions_mat = stack(sat_positions)

    # Primary system — its clock bias, reference time, week and start epoch
    # define the reported time. A GGTO collapse (signalled by fewer clock biases
    # than systems) is anchored on GPS, so GPS must be primary there; otherwise
    # pick the system with the most satellites (best-conditioned reported time),
    # breaking ties by first appearance.
    unique_systems = unique(systems)
    primary_system =
        num_clock_biases < length(unique_systems) ? GPST() :
        unique_systems[argmax([count(==(sys), systems) for sys in unique_systems])]
    primary_clock_index = clock_bias_indices[findfirst(==(primary_system), systems)]

    # The common reference cancels out of the reported time (the primary clock
    # bias absorbs it), so any latest-transmit-time reference works.
    pseudo_ranges, reference_time = calc_pseudo_ranges(times)
    # Apply the known per-satellite GGTO time-system offset (zero unless Galileo was
    # collapsed onto GPS) as a measurement correction, like the atmospheric delays
    # below. Kept for the inter-system-bias readout at the end.
    ggto_offsets = calc_ggto_range_offsets(ggto_decoder, systems, times)
    pseudo_ranges = pseudo_ranges .- ggto_offsets

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias.
    prev_abs_bias(sys) =
        sys == prev_pvt.reference_system ? prev_pvt.time_correction :
        prev_pvt.time_correction + get(prev_pvt.inter_system_biases, sys, 0.0m)
    prev_ξ = zeros(num_lsq_params(bias_columns))
    prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
    for j in 1:num_sats
        prev_ξ[3+clock_bias_indices[j]] = ustrip(m, prev_abs_bias(systems[j]))
    end
    # Only warm-start an IFB column from the previous solution when that band's bias was
    # measured against the same reference band; a reference change (the anchor differs
    # across epochs) makes the stored value refer to a different quantity, so seed at 0
    # instead. The IFB enters the design linearly, so a stale seed only costs iterations,
    # but this keeps the starting point meaningful.
    for (i, band) in enumerate(extra_bands)
        prev_ifb = get(prev_pvt.inter_frequency_biases, band, nothing)
        prev_ξ[3+num_clock_biases+i] =
            !isnothing(prev_ifb) && prev_ifb.reference == reference_bands[i] ?
            ustrip(m, prev_ifb.value) : 0.0
    end

    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model is chosen for the whole solve. The
    # prediction is the top-level `predict_atmospheric_delays`, a function
    # barrier that specialises on the concrete type of the `Union`-typed
    # `ionospheric_correction`.
    ionospheric_correction =
        enable_ionospheric_correction ? select_ionospheric_correction(healthy_states) :
        nothing
    # Is any atmospheric correction active at all? When neither the ionosphere (no
    # model selected) nor the troposphere contributes, skip `predict_atmospheric_delays`
    # entirely and the solve runs on the raw pseudoranges.
    correct_atmosphere =
        !isnothing(ionospheric_correction) || enable_tropospheric_correction

    # Weight the solve by measurement uncertainty as soon as any satellite carries the
    # information for it (a C/N₀ or an explicit variance); with none, every weight would
    # be the same and the epoch is solved by ordinary least squares exactly as before.
    # See `has_measurement_uncertainty` and [`pseudorange_variance`](@ref).
    weighted = any(has_measurement_uncertainty, healthy_states)

    # One solve about the measurement model evaluated at `ξ_reference`: both the
    # atmospheric delays and the a-priori variances depend on position only through the
    # satellite elevations, so a metre-accurate reference position is enough for either
    # (see `predict_atmospheric_delays` / `predict_pseudorange_variances`) and neither
    # needs iterating to convergence. Returns the solved state, the raw (metre) residuals
    # and the variances the weights came from, which are reported per satellite and set
    # the scale of the formal accuracy below.
    function solve_about(ξ_reference)
        ranges =
            correct_atmosphere ?
            pseudo_ranges .- predict_atmospheric_delays(
                ξ_reference, healthy_states, sat_positions, ionospheric_correction,
                reference_time, enable_tropospheric_correction) : pseudo_ranges
        variances =
            weighted ?
            predict_pseudorange_variances(ξ_reference, healthy_states, sat_positions) :
            fill(_NOMINAL_PSEUDORANGE_SIGMA^2, num_sats)
        ξ_solved, residuals = user_position(
            sat_positions_mat, ranges, bias_columns, ξ_reference;
            weights = weighted ? inv.(variances) : nothing)
        (ξ_solved, residuals, variances)
    end

    ξ, residuals, pseudorange_variances = if iszero(prev_ξ)
        # Cold start: no prior position, and both the Klobuchar model and the elevation
        # part of the variance model need one, so first obtain an approximate fix from an
        # uncorrected, unweighted solve, then re-solve once about it (only if there is
        # anything to correct or to weight, so the plain case stays a single solve).
        ξ_uncorrected, resid_uncorrected =
            user_position(sat_positions_mat, pseudo_ranges, bias_columns, prev_ξ)
        if correct_atmosphere || weighted
            solve_about(ξ_uncorrected)
        else
            (ξ_uncorrected, resid_uncorrected, fill(_NOMINAL_PSEUDORANGE_SIGMA^2, num_sats))
        end
    else
        # Warm start: the previous (already metre-accurate) position is the reference, so
        # ξ never needs a post-solve correction and one solve suffices.
        solve_about(prev_ξ)
    end
    H = calc_H(sat_positions_mat, ξ, bias_columns)
    position = ECEF(ξ[1], ξ[2], ξ[3])

    # Check the geometry at the converged position — the DOP is reported to the caller, and
    # a rank deficiency here (negative GDOP) means `ξ` is meaningless, so nothing should be
    # derived from it. The satellite conditions in `decide_bias_layout` already reject the
    # count-shaped degeneracies before the solve; this catches what a count cannot see.
    #
    # This check must stay ahead of the velocity solve: it is what makes that solve's
    # normal-equations matrix positive definite (see `calc_user_velocity_and_clock_drift`),
    # and with the order reversed a degenerate geometry throws `SingularException` there.
    dop = calc_DOP(H, position, primary_clock_index)
    dop.GDOP < 0 && return prev_pvt

    # The DOP stays the purely geometric `(HᵀH)⁻¹`; the formal accuracy is the separate,
    # metre-valued `(HᵀWH)⁻¹` of the solve that actually ran (see `FormalAccuracy`). Its
    # positive definiteness follows from the DOP check just made, weights being positive.
    accuracy = calc_formal_accuracy(H, pseudorange_variances, position, primary_clock_index)

    # The velocity solve is weighted from the Doppler's own variance model — carrier-loop
    # noise, not code-loop noise — and left unweighted when no C/N₀ was reported.
    range_rate_weights =
        weighted ? inv.(map(range_rate_variance, healthy_states)) : nothing
    user_velocity_and_clock_drift, rate_residuals = calc_user_velocity_and_clock_drift(
        sat_positions_and_velocities, healthy_states, times, H;
        weights = range_rate_weights)
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    relative_clock_drift = user_velocity_and_clock_drift[4] / SPEEDOFLIGHT
    course_over_ground = calc_course_over_ground(position, velocity)
    time_correction = ξ[3+primary_clock_index]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEEDOFLIGHT

    primary_state = healthy_states[findfirst(==(primary_system), systems)]
    week = get_week(primary_state.decoder; approximate_year)
    start_time = system_start_epoch(primary_state.system)
    # Assumes `start_time.fraction == 0` (true for GPS/Galileo: integer-second origins).
    time = TAIEpoch(
        week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second,
        corrected_reference_time - floor(Int, corrected_reference_time),
    )

    sat_infos = SatInfo.(sat_positions, times, residuals .* m, rate_residuals .* (m/s),
        sqrt.(pseudorange_variances) .* m)

    # Inter-system biases relative to the reference (primary) system's clock, in
    # meters. The reference is omitted (its bias is `time_correction`); for a
    # GGTO-collapsed system this is the broadcast offset −c·Δt_systems, read from the
    # system's first satellite (the per-satellite offsets differ only by the GGTO
    # drift A_1G over the transmit-time spread — sub-millimetre).
    inter_system_biases = Dict{GNSSSignals.TimeSystem,typeof(1.0m)}()
    for sys in unique_systems
        sys == primary_system && continue
        j = findfirst(==(sys), systems)
        inter_system_biases[sys] =
            (ξ[3+clock_bias_indices[j]] + ggto_offsets[j] - time_correction) * m
    end

    # Receiver inter-frequency biases relative to the reference band, in meters
    # (the reference band is omitted; its bias is folded into the clock biases).
    inter_frequency_biases = Dict{Symbol,InterFrequencyBias}()
    for (i, band) in enumerate(extra_bands)
        inter_frequency_biases[band] =
            InterFrequencyBias(ξ[3+num_clock_biases+i] * m, reference_bands[i])
    end

    # Per-satellite `sats` key: (signal id, PRN) — signal-level (not time system), so a
    # satellite tracked on two signals of one constellation stays distinct; the
    # receiver-clock grouping is separate, by time system.
    healthy_sat_keys =
        map(state -> (get_signal_id(state.system), state.decoder.prn), healthy_states)

    PVTSolution(
        position,
        velocity,
        course_over_ground,
        time_correction * m,
        time,
        relative_clock_drift,
        dop,
        accuracy,
        Dictionary(healthy_sat_keys, sat_infos),
        primary_system,
        inter_system_biases,
        inter_frequency_biases,
    )
end

"""
    system_start_epoch(system) -> TAIEpoch

Absolute TAI epoch of a ranging signal's GNSS time-scale origin (week 0, time of
week 0). `get_system_start_time` gives the origin as a UTC `DateTime`;
AstroTime's leap-aware UTC→TAI conversion places it on the atomic scale (GPS
`1980-01-06T00:00:19` TAI, Galileo `1999-08-22T00:00:19` TAI). The leap-second count
lives in AstroTime, not in this package.
"""
system_start_epoch(system) =
    from_utc(get_system_start_time(system); scale = TAI)

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
    # that `cycle_base + trans_week` is closest to that anchor.
    days_at_anchor = Date(approximate_year, 6, 30) - Date(1980, 1, 6)
    weeks_at_anchor = Dates.value(days_at_anchor) ÷ 7
    n_cycles = round(Int, (weeks_at_anchor - decoder.data.trans_week) / 1024)
    return n_cycles * 1024 + decoder.data.trans_week
end

# Galileo I/NAV and F/NAV (12-bit WN) and GPS CNAV (L5, L2C) / CNAV-2 (L1C) (full
# 13-bit WN) all broadcast the absolute week number, so — unlike GPS L1 C/A LNAV —
# there is no 1024-week rollover to resolve and `approximate_year` is unused.
function get_week(
    decoder::GNSSDecoder.GNSSDecoderState{<:Union{GNSSDecoder.AbstractGalileoData,GPSModernNavData}};
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
include("measurement_variance.jl")
end
