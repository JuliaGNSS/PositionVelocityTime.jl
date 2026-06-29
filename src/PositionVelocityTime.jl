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

using Unitful: s, Hz
using Dictionaries: Dictionary

export calc_pvt,
    PVTSolution,
    SatInfo,
    SatelliteState,
    get_LLA,
    get_num_used_sats,
    get_sat_info,
    calc_satellite_position,
    calc_satellite_position_and_velocity,
    get_sat_enu,
    get_frequency_offset

const SPEEDOFLIGHT = 299792458.0

# Galileo navigation messages that share the same Keplerian ephemeris, clock,
# NTCM-G ionospheric coefficients and GGTO layout: I/NAV (decoded from E1B) and
# F/NAV (decoded from E5a-I). They are referenced to the same Galileo System Time,
# so the receiver-clock and correction handling is identical; only the broadcast
# group delay differs per signal (see `correct_by_group_delay`).
const GalileoNavData = Union{GNSSDecoder.GalileoE1BData,GNSSDecoder.GalileoE5aData}

# All GPS civil navigation messages: they share the GPS time system and the
# single-frequency Klobuchar ionospheric model (LNAV on L1 C/A, CNAV on L5 and
# L2C, CNAV-2 on L1C). `GNSSDecoder.GPSCNAVData` is the shared CNAV container for
# both L5-I and L2C-M — the signal is fixed by the decoder's constants type, not
# the data type — so it covers both.
const GPSNavData =
    Union{GNSSDecoder.GPSL1CAData,GNSSDecoder.GPSCNAVData,GNSSDecoder.GPSL1C_DData}

# The GPS CNAV family (CNAV on L5 and L2C, CNAV-2 on L1C): all broadcast the full
# week number (no 1024-week rollover) and the quasi-Keplerian ephemeris (A_REF + ΔA,
# Ω̇_REF + ΔΩ̇, …) rather than the directly-broadcast Keplerian elements of LNAV.
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
- `code_phase::CP`: Code phase measurement
- `carrier_doppler`: Carrier Doppler frequency in Hz
- `carrier_phase::CP`: Carrier phase measurement (default: `0.0`)

# Constructors
    SatelliteState(; decoder, system, code_phase, carrier_doppler, carrier_phase=0.0)
    SatelliteState(decoder, system, sat_state)

The second constructor extracts code phase, carrier Doppler, and carrier phase from a
`Tracking` satellite state (`Tracking.SatState` for Tracking ≤ 1, `Tracking.TrackedSat` for
Tracking ≥ 2). It is provided by a package extension that is loaded automatically once
`Tracking` is available, so `Tracking` is only a weak dependency of this package.
"""
@kwdef struct SatelliteState{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    system::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
    carrier_phase::CP = 0.0
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
    SatInfo

Per-satellite information attached to a [`PVTSolution`](@ref) (one entry per
satellite used in the fix).

# Fields
- `position::ECEF`: Satellite ECEF position at transmit time (metres).
- `time::Float64`: Satellite transmit time (system time of week, seconds).
- `residual::Float64`: Post-fit least-squares pseudorange residual (metres) — the
  modeled minus the (atmosphere-corrected) measured pseudorange. A per-satellite
  fit-quality / outlier indicator.
"""
struct SatInfo
    position::ECEF
    time::Float64
    residual::Float64
end

"""
    PVTSolution

Complete Position, Velocity, and Time solution from GNSS measurements.

# Fields
- `position::ECEF`: User position in ECEF coordinates (meters)
- `velocity::ECEF`: User velocity in ECEF coordinates (m/s)
- `time_correction::Float64`: Estimated receiver clock bias of the reference GNSS
  (meters). For a multi-GNSS solution this is the bias of `reference_system`;
  other systems' biases are `time_correction + inter_system_biases[system]`.
- `time::Union{TAIEpoch{Float64}, Nothing}`: Estimated time as a TAI epoch
- `relative_clock_drift::Float64`: Relative receiver clock drift (dimensionless)
- `dop::Union{DOP, Nothing}`: Dilution of precision values
- `sats::Dictionary{Tuple{Symbol, Int}, SatInfo}`: Maps `(signal, PRN)` to satellite
  info — position, transmit time, and post-fit residual (see [`SatInfo`](@ref)). The
  signal tag (e.g. `:GPSL1CA`, `:GalileoE1B`; see [`get_signal_id`](@ref)) keeps the
  same PRN apart both across constellations (GPS PRN 5 vs Galileo E05) and across
  signals of one constellation (a satellite tracked on GPS L1 C/A and L5 yields two
  entries sharing a PRN). Receiver-clock grouping is by time system, not signal —
  see `reference_system`.
- `reference_system::Union{Symbol, Nothing}`: GNSS time system (e.g. `:GPS`,
  `:Galileo`) that `time` and `time_correction` are referenced to.
- `inter_system_biases::Dict{Symbol, Float64}`: For each GNSS time system other than
  `reference_system`, the offset of that system's time scale relative to the reference
  system's (meters) — the inter-system bias. This is a **system / space-segment** effect,
  not a receiver one (this receiver has no inter-system hardware bias): it is the GNSS
  system-time offset, so the Galileo entry equals `−c · Δt_systems`, `Δt_systems = GST −
  GPST` (the GGTO). It is **estimated directly from the geometry whenever observable**
  (no broadcast error); the broadcast GGTO is used to derive it only as a fallback (the
  GGTO-aided collapse), as that broadcast value may be erroneous. Reference-independent
  (the difference of two entries is the offset between those two systems); empty for a
  single-system solution.
- `inter_frequency_biases::Dict{Symbol, Float64}`: For each frequency band other than
  the reference band, the receiver inter-frequency bias relative to it (meters) — the
  differential hardware delay of that band's RF chain, estimated as an extra unknown
  when satellites are processed on more than one band. The key is the band
  (e.g. `:L1`, `:L5`; see [`get_frequency_band`](@ref)) and the value is shared across
  all constellations on that band. Empty for a single-band solution. A solution then
  needs `n ≥ 3 + M + B` satellites for `M` time systems and `B` extra bands.
"""
@kwdef struct PVTSolution
    position::ECEF = ECEF(0, 0, 0)
    velocity::ECEF = ECEF(0, 0, 0)
    time_correction::Float64 = 0
    time::Union{TAIEpoch{Float64},Nothing} = nothing
    relative_clock_drift::Float64 = 0
    dop::Union{DOP,Nothing} = nothing
    sats::Dictionary{Tuple{Symbol,Int},SatInfo} = Dictionary{Tuple{Symbol,Int},SatInfo}()
    reference_system::Union{Symbol,Nothing} = nothing
    inter_system_biases::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
    inter_frequency_biases::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
end

"""
    get_num_used_sats(pvt_solution::PVTSolution) -> Int

Return the number of satellites used in the PVT solution.
"""
function get_num_used_sats(pvt_solution::PVTSolution)
    length(pvt_solution.sats)
end

"""
    get_sat_info(pvt_solution::PVTSolution, signal::Symbol, prn::Integer) -> Union{SatInfo,Nothing}

Return the [`SatInfo`](@ref) (position, transmit time and post-fit residual) of the
satellite with the given `prn` on GNSS `signal` (e.g. `:GPSL1CA`, `:GalileoE1B`; see
[`get_signal_id`](@ref)), or `nothing` if that satellite was not used in the fix. The
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
    get_time_system(decoder) -> Symbol

Return the GNSS time system a decoder's measurements are referenced to
(`:GPS` or `:Galileo`). Satellites sharing a time system share one receiver
clock bias in the PVT estimation.
"""
get_time_system(::GNSSDecoder.GNSSDecoderState{<:GPSNavData}) = :GPS
get_time_system(::GNSSDecoder.GNSSDecoderState{<:GalileoNavData}) = :Galileo

"""
    get_signal_id(state::SatelliteState) -> Symbol

Signal-level identifier of a satellite measurement — the name of its GNSS signal
type (e.g. `:GPSL1CA`, `:GalileoE1B`, `:GalileoE5aI`). This is the per-signal part of
the `sats` key in [`PVTSolution`](@ref), so the same PRN measured on two different
bands of one constellation (e.g. GPS L1 and L5, or Galileo E1 and E5a) yields two
distinct entries.

Distinct from [`get_time_system`](@ref): those two signals share a time system —
and therefore one receiver clock bias — but have different signal ids. The signal
type carries its (large) code table, so its name is used as the key rather than the
instance itself.
"""
get_signal_id(state::SatelliteState) = nameof(typeof(state.system))

"""
    get_frequency_band(state::SatelliteState) -> Symbol

Frequency band of a satellite measurement (e.g. `:L1`, `:L5`), from the GNSSSignals
band of `state.system`. Satellites sharing a band share one receiver inter-frequency
bias in the PVT estimation, regardless of constellation: GPS L1 C/A, GPS L1C and
Galileo E1 are all `:L1`; GPS L5 and Galileo E5a are all `:L5`. Distinct from both
[`get_time_system`](@ref) (which groups the receiver *clock* bias) and
[`get_signal_id`](@ref) (the per-signal `sats` key).
"""
get_frequency_band(state::SatelliteState) = nameof(typeof(GNSSSignals.get_band(state.system)))

"""
    band_ifb_layout(systems, bands) -> (ifb_indices, extra_bands, num_components)

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
deterministically); `num_components` is the number of coverage components (`1` ⇔ the
graph is connected).
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
    band_column = Dict(b => i for (i, b) in enumerate(extra_bands))
    ifb_indices = [get(band_column, b, 0) for b in bands]
    return ifb_indices, extra_bands, length(reference_of)
end

"""
    decide_bias_layout(states, systems, bands, times)
        -> (clock_bias_indices, num_clock_biases, ifb_indices, extra_bands, ggto_offsets)

Decide the full least-squares bias layout — one clock column per GNSS time system plus
the per-band inter-frequency-bias columns from [`band_ifb_layout`](@ref) — and the
known per-satellite range offsets. Returns `nothing` when the constellation cannot be
solved. The decision is observability-driven, not merely count-driven:

- When the coverage graph is connected and `n ≥ 3 + num_systems + num_ifb`, estimate
  everything independently: the inter-system offset and the receiver inter-frequency
  biases are observed directly from the geometry, so neither inherits the broadcast-GGTO
  error (the satellite group delays are already removed per satellite upstream, so the
  per-band column carries the receiver chain).
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
layout yields a full-rank design matrix — the degenerate disjoint-band case is removed
structurally, not caught after the fact.
"""
function decide_bias_layout(states, systems, bands, times)
    num_sats = length(states)
    enough_satellites(num_clock, num_ifb) = num_sats >= 3 + num_clock + num_ifb

    function bias_layout_for(effective_systems)
        unique_effective = unique(effective_systems)
        index = Dict(sys => i for (i, sys) in enumerate(unique_effective))
        clock_bias_indices = [index[sys] for sys in effective_systems]
        ifb_indices, extra_bands, num_components = band_ifb_layout(effective_systems, bands)
        (; clock_bias_indices, num_clock_biases = length(unique_effective), ifb_indices,
            extra_bands, num_components)
    end

    independent_layout = bias_layout_for(systems)
    if independent_layout.num_components == 1 &&
       enough_satellites(independent_layout.num_clock_biases, length(independent_layout.extra_bands))
        return independent_layout.clock_bias_indices, independent_layout.num_clock_biases,
            independent_layout.ifb_indices, independent_layout.extra_bands, zeros(num_sats)
    end

    # Connected-but-scarce or disconnected: try the GGTO collapse (merge Galileo onto
    # GPS). The GGTO is the same constellation-wide offset whichever Galileo satellite
    # reports it, so one decoded copy converts every Galileo measurement.
    ggto_idx = findfirst(
        j -> systems[j] == :Galileo && ggto_available(states[j].decoder), 1:num_sats)
    if (:GPS in systems) && !isnothing(ggto_idx)
        merged_layout = bias_layout_for(map(sys -> sys == :Galileo ? :GPS : sys, systems))
        if enough_satellites(merged_layout.num_clock_biases, length(merged_layout.extra_bands))
            # Per the OS SIS ICD the GGTO is GST − GPST, so a Galileo transmit time
            # becomes GPS time by SUBTRACTING it; the modeled range carries −c·GGTO and
            # the solve yields inter_system_biases[:Galileo] = −c·(GST − GPST).
            ggto_decoder = states[ggto_idx].decoder
            ggto_offsets = zeros(num_sats)
            for j in 1:num_sats
                systems[j] == :Galileo || continue
                ggto_offsets[j] = -SPEEDOFLIGHT * calc_ggto_offset(ggto_decoder, times[j])
            end
            return merged_layout.clock_bias_indices, merged_layout.num_clock_biases, merged_layout.ifb_indices,
                merged_layout.extra_bands, ggto_offsets
        end
    end

    # No collapse available. The independent layout is still observable (its IFBs are
    # component-restricted); use it if there are enough satellites, otherwise unsolvable.
    enough_satellites(independent_layout.num_clock_biases, length(independent_layout.extra_bands)) ?
    (independent_layout.clock_bias_indices, independent_layout.num_clock_biases, independent_layout.ifb_indices,
        independent_layout.extra_bands, zeros(num_sats)) : nothing
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
[`get_frequency_band`](@ref)). The state vector is therefore
`[x, y, z, tc₁, …, tc_M, ifb₁, …, ifb_B]` for `M` distinct time systems and `B`
extra bands. Position and time are found by least squares; velocity and clock drift
are solved from carrier Doppler.

A solution requires `n ≥ 3 + M + B` healthy satellites (each system needs at least
one satellite, and a system contributing a single satellite spends it entirely
on that system's clock bias). When that condition fails but both GPS and Galileo
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

# Arguments
- `states`: Vector of [`SatelliteState`](@ref) for observed satellites
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
A [`PVTSolution`](@ref) containing position, velocity, time, DOP values, and
satellite information. Returns `prev_pvt` if too few healthy satellites are
available to solve the constellation (including the GGTO fallback) or if the
GDOP is negative.

# Throws
- `ArgumentError`: If fewer than 4 satellite states are provided
"""
function calc_pvt(
    states::AbstractVector{<:SatelliteState},
    prev_pvt::PVTSolution = PVTSolution();
    approximate_year::Integer = year(now(UTC)),
    enable_ionospheric_correction::Bool = true,
    enable_tropospheric_correction::Bool = true,
)
    length(states) < 4 &&
        throw(ArgumentError("You'll need at least 4 satellites to calculate PVT"))
    healthy_indices = findall(x -> is_sat_healthy(x.decoder), states)
    length(healthy_indices) < 4 && return prev_pvt
    healthy_states = view(states, healthy_indices)
    num_sats = length(healthy_states)

    times = map(calc_corrected_time, healthy_states)

    # Group satellites by GNSS time system (ordered by first appearance) and by frequency
    # band, then decide the full bias layout: one clock column per time system, plus a
    # receiver inter-frequency-bias column per
    # band beyond a per-coverage-component reference (`decide_bias_layout` keeps only
    # observable IFBs and falls back to a GGTO collapse when the geometry is
    # disconnected or under-determined).
    systems = map(state -> get_time_system(state.decoder), healthy_states)
    bands = map(get_frequency_band, healthy_states)

    bias_layout = decide_bias_layout(healthy_states, systems, bands, times)
    isnothing(bias_layout) && return prev_pvt
    clock_bias_indices, num_clock_biases, ifb_indices, extra_bands, ggto_offsets = bias_layout
    num_ifb = length(extra_bands)
    bias_columns = BiasColumns(clock_bias_indices, num_clock_biases, ifb_indices, num_ifb)

    # Propagating the ephemerides.
    sat_positions_and_velocities = map(
        (state, time) -> calc_satellite_position_and_velocity(state.decoder, time),
        healthy_states,
        times,
    )
    sat_positions = map(get_sat_position, sat_positions_and_velocities)
    # Built once here and reused by every least-squares solve below and by calc_H,
    # instead of being rebuilt inside each user_position call.
    sat_positions_mat = reduce(hcat, sat_positions)

    # Primary system — its clock bias, reference time, week and start epoch
    # define the reported time. A GGTO collapse (signalled by fewer clock biases
    # than systems) is anchored on GPS, so GPS must be primary there; otherwise
    # pick the system with the most satellites (best-conditioned reported time),
    # breaking ties by first appearance.
    unique_systems = unique(systems)
    primary_system =
        num_clock_biases < length(unique_systems) ? :GPS :
        unique_systems[argmax([count(==(sys), systems) for sys in unique_systems])]
    primary_clock_index = clock_bias_indices[findfirst(==(primary_system), systems)]

    # The common reference cancels out of the reported time (the primary clock
    # bias absorbs it), so any latest-transmit-time reference works.
    pseudo_ranges, reference_time = calc_pseudo_ranges(times)
    # Apply the known per-satellite GGTO time-system offset (zero unless Galileo was
    # collapsed onto GPS) as a measurement correction, like the atmospheric delays
    # below.
    pseudo_ranges = pseudo_ranges .- ggto_offsets

    # Seed each clock bias from the previous solution, reconstructing a system's
    # absolute bias from the reference bias plus its stored inter-system bias.
    prev_abs_bias(sys) =
        sys == prev_pvt.reference_system ? prev_pvt.time_correction :
        prev_pvt.time_correction + get(prev_pvt.inter_system_biases, sys, 0.0)
    prev_ξ = zeros(num_lsq_params(bias_columns))
    prev_ξ[1], prev_ξ[2], prev_ξ[3] = prev_pvt.position
    for j in 1:num_sats
        prev_ξ[3+clock_bias_indices[j]] = prev_abs_bias(systems[j])
    end
    for (i, band) in enumerate(extra_bands)
        prev_ξ[3+num_clock_biases+i] = get(prev_pvt.inter_frequency_biases, band, 0.0)
    end

    # Atmospheric corrections, summed per satellite and subtracted from the
    # pseudoranges. The ionospheric model is chosen for the whole solve from the
    # coefficients decoded across the constellation (NTCM-G if Galileo coefficients
    # are available, else Klobuchar if GPS α/β are available, else none; see
    # `select_ionospheric_correction`). The prediction is the top-level
    # `predict_atmospheric_delays`, a function barrier that specialises on the
    # concrete type of the `Union`-typed `ionospheric_correction`.
    ionospheric_correction =
        enable_ionospheric_correction ? select_ionospheric_correction(healthy_states) :
        nothing
    ξ, residuals = if iszero(prev_ξ)
        # Cold start: no prior position, and the Klobuchar model is undefined near
        # the geocenter, so first obtain an approximate fix from an uncorrected
        # solve, then re-solve once with the delay-corrected pseudoranges (only if
        # there is anything to correct, so the uncorrected case stays a single solve).
        ξ_uncorrected, resid_uncorrected =
            user_position(sat_positions_mat, pseudo_ranges, bias_columns, prev_ξ)
        atmospheric_delays = predict_atmospheric_delays(
            ξ_uncorrected, healthy_states, sat_positions, ionospheric_correction,
            reference_time, enable_tropospheric_correction)
        any(!iszero, atmospheric_delays) ?
        user_position(
            sat_positions_mat, pseudo_ranges .- atmospheric_delays, bias_columns, ξ_uncorrected) :
        (ξ_uncorrected, resid_uncorrected)
    else
        # Warm start: predict the delays from the previous (already metre-accurate)
        # position before solving, so ξ never needs a post-solve correction.
        user_position(
            sat_positions_mat,
            pseudo_ranges .- predict_atmospheric_delays(
                prev_ξ, healthy_states, sat_positions, ionospheric_correction,
                reference_time, enable_tropospheric_correction),
            bias_columns,
            prev_ξ,
        )
    end
    H = calc_H(sat_positions_mat, ξ, bias_columns)
    # H feeds both the per-system DOP and the velocity/clock-drift solve, which
    # reuses its line-of-sight block and collapses its clock columns to a single
    # common receiver clock drift.
    user_velocity_and_clock_drift = calc_user_velocity_and_clock_drift(
        sat_positions_and_velocities, healthy_states, times, H)
    position = ECEF(ξ[1], ξ[2], ξ[3])
    velocity = ECEF(
        user_velocity_and_clock_drift[1],
        user_velocity_and_clock_drift[2],
        user_velocity_and_clock_drift[3],
    )
    relative_clock_drift = user_velocity_and_clock_drift[4] / SPEEDOFLIGHT
    time_correction = ξ[3+primary_clock_index]
    # The estimated time correction is negative
    # See https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/8
    corrected_reference_time = reference_time - time_correction / SPEEDOFLIGHT

    primary_decoder = healthy_states[findfirst(==(primary_system), systems)].decoder
    week = get_week(primary_decoder; approximate_year)
    start_time = get_system_start_time(primary_decoder)
    time = TAIEpoch(
        week * 7 * 24 * 60 * 60 + floor(Int, corrected_reference_time) + start_time.second,
        corrected_reference_time - floor(Int, corrected_reference_time),
    )

    sat_infos = SatInfo.(sat_positions, times, residuals)

    dop = calc_DOP(H, position, primary_clock_index)
    if dop.GDOP < 0
        return prev_pvt
    end

    # Inter-system biases relative to the reference (primary) system's clock, in
    # meters. The reference is omitted (its bias is `time_correction`); for a
    # GGTO-collapsed system this is the broadcast offset c·Δt_systems.
    inter_system_biases = Dict{Symbol,Float64}()
    for j in 1:num_sats
        systems[j] == primary_system && continue
        inter_system_biases[systems[j]] =
            ξ[3+clock_bias_indices[j]] + ggto_offsets[j] - time_correction
    end

    # Receiver inter-frequency biases relative to the reference band, in meters
    # (the reference band is omitted; its bias is folded into the clock biases).
    inter_frequency_biases = Dict{Symbol,Float64}()
    for (i, band) in enumerate(extra_bands)
        inter_frequency_biases[band] = ξ[3+num_clock_biases+i]
    end

    # Per-satellite `sats` key: (signal id, PRN) — signal-level (not time system), so a
    # satellite tracked on two signals of one constellation stays distinct; the
    # receiver-clock grouping is separate, by time system.
    healthy_sat_keys =
        map(state -> (get_signal_id(state), state.decoder.prn), healthy_states)

    PVTSolution(
        position,
        velocity,
        time_correction,
        time,
        relative_clock_drift,
        dop,
        Dictionary(healthy_sat_keys, sat_infos),
        primary_system,
        inter_system_biases,
        inter_frequency_biases,
    )
end

"""
    get_frequency_offset(pvt::PVTSolution, base_frequency) -> typeof(base_frequency)

Calculate the receiver frequency offset from the relative clock drift and a base frequency.

# Arguments
- `pvt::PVTSolution`: PVT solution containing the relative clock drift
- `base_frequency`: Reference frequency (e.g., the carrier frequency of the GNSS signal)

# Returns
The frequency offset as `relative_clock_drift * base_frequency`.
"""
function get_frequency_offset(pvt::PVTSolution, base_frequency)
    pvt.relative_clock_drift * base_frequency
end

function get_system_start_time(decoder::GNSSDecoder.GNSSDecoderState{<:GPSNavData})
    TAIEpoch(1980, 1, 6, 0, 0, 19.0) # There were 19 leap seconds at 01/06/1999 compared to UTC
end

function get_system_start_time(
    decoder::GNSSDecoder.GNSSDecoderState{<:GalileoNavData},
)
    TAIEpoch(1999, 8, 22, 0, 0, (32 - 13.0)) # There were 32 leap seconds at 08/22/1999 compared to UTC
end

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
    decoder::GNSSDecoder.GNSSDecoderState{<:Union{GalileoNavData,GPSModernNavData}};
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
end
