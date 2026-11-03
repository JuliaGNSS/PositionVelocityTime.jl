# The receiver bias columns of the least-squares state: one clock per time system, the
# inter-frequency biases, the hub that systems tied by a broadcast offset collapse onto,
# and the reference system.

"""
    BiasColumns

Per-satellite assignment of the estimated bias columns of the least-squares design
matrix, shared by [`predict_pseudoranges!`](@ref), [`calc_design_matrix!`](@ref) and
[`estimate_position`](@ref). The state vector is
`[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifbs]` with
[`num_parameters`](@ref)`(bias_columns)` entries. The two column kinds have different
physical sources: a clock column is the receiver clock for one GNSS time system (the
spacing *between* systems is a system/space-segment effect — the broadcast GNSS time
offset, GGTO or BGTO), whereas an inter-frequency-bias column is the receiver's per-band
RF-chain delay. Known per-satellite corrections (atmosphere, the broadcast time-system
offset) are applied to the pseudoranges in [`calc_pvt`](@ref), not carried here.

# Fields
- `clock_columns::Vector{Int}`: per satellite, the clock column (1…`num_clock_biases`)
  of its GNSS time system; the design-matrix `1.0` lands at `3 + clock_columns[j]`.
- `num_clock_biases::Int`: number of clock biases (also the offset of the IFB block).
- `ifb_columns::Vector{Int}`: per satellite, the inter-frequency-bias column
  (1…`num_ifbs`) of its frequency band, or `0` for the reference band; the `1.0` lands
  at `3 + num_clock_biases + ifb_columns[j]`.
- `num_ifbs::Int`: number of inter-frequency biases (frequency bands beyond the reference).
"""
struct BiasColumns
    clock_columns::Vector{Int}
    num_clock_biases::Int
    ifb_columns::Vector{Int}
    num_ifbs::Int
end

"""
    num_parameters(bias_columns::BiasColumns) -> Int

Length of the least-squares state vector for `bias_columns`:
`3 + num_clock_biases + num_ifbs`.
"""
num_parameters(bias_columns::BiasColumns) =
    3 + bias_columns.num_clock_biases + bias_columns.num_ifbs

"""
    ifb_layout(constellation_keys, bands)
        -> (ifb_indices, extra_bands, reference_bands, num_components)

Lay out the receiver inter-frequency biases from the (constellation × band) coverage
graph. `constellation_keys[j]` and `bands[j]` are satellite `j`'s constellation and
frequency band; the constellation enters only through equality and first appearance, so any
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
function ifb_layout(constellation_keys, bands)
    ifb_indices = Int[]
    extra_bands = Vector{eltype(bands)}()
    reference_bands = Vector{eltype(bands)}()
    num_components = ifb_layout!(ifb_indices, extra_bands,
        reference_bands,
        IFBLayoutScratch{eltype(bands)}(), constellation_keys, bands)
    return ifb_indices, extra_bands, reference_bands, num_components
end

# The scratch vectors of `ifb_layout!`. Bands are identified by their position in
# `unique_bands` (first appearance), so the union-find, the counts and the per-component
# references are all plain `Int` vectors — the same bookkeeping a `Dict` keyed by band
# would hold, without a hash table to allocate per epoch.
struct IFBLayoutScratch{B}
    unique_bands::Vector{B}
    band_index::Vector{Int}   # per satellite: its band's position in `unique_bands`
    parent::Vector{Int}       # per band: union-find parent
    band_count::Vector{Int}   # per band: satellites on it
    reference::Vector{Int}    # per union-find root: the component's reference band, or 0
    column::Vector{Int}       # per band: its IFB column, or 0 for a reference band
end

IFBLayoutScratch{B}() where {B} =
    IFBLayoutScratch{B}(B[], Int[], Int[], Int[], Int[], Int[])

"""
    ifb_layout!(ifb_indices, extra_bands, reference_bands,
        scratch, constellation_keys, bands)
        -> num_components

[`ifb_layout`](@ref) writing its three vectors into the given ones (emptied first)
and working in `scratch`, so that it allocates nothing once they have grown to the epoch.
The result is identical: the same components, the same references and tie-breaks, and
the same column order.
"""
function ifb_layout!(
    ifb_indices,
    extra_bands,
    reference_bands,
    scratch::IFBLayoutScratch,
    constellation_keys,
    bands,
)
    (; unique_bands, band_index, parent, band_count, reference, column) = scratch
    num_observations = length(bands)
    empty!(unique_bands)
    resize!(band_index, num_observations)
    for j in 1:num_observations
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
    for j in 1:num_observations
        first_of_constellation = something(
            findfirst(i -> isequal(constellation_keys[i], constellation_keys[j]), 1:j))
        a = find_root!(parent, band_index[first_of_constellation])
        c = find_root!(parent, band_index[j])
        a == c || (parent[a] = c)
    end
    resize!(band_count, num_bands)
    fill!(band_count, 0)
    for j in 1:num_observations
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
    resize!(ifb_indices, num_observations)
    for j in 1:num_observations
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
    clock_systems::Vector{SupportedTimeSystem},
    extra_bands::Vector{Symbol},
    reference_bands::Vector{Symbol},
    hub_system::Union{Nothing,SupportedTimeSystem},
    hub_rows::Vector{PreparedObservation},
}

"""
    decide_bias_layout(prepared, time_offsets::BroadcastTimeOffsets)
        -> Union{BiasLayout,Nothing}

Decide the full least-squares bias layout — one clock column per GNSS time system plus
the per-band inter-frequency-bias columns from [`ifb_layout`](@ref) — for the flat
[`PreparedObservation`](@ref) rows of one epoch, and return it as
a [`BiasLayout`](@ref),
or `nothing` when the constellation cannot be solved. Every classification the decision
needs — each satellite's `time_system`, `band_id` and `prn` — is a field of the row, and
the broadcast offsets between the systems are the epoch's `time_offsets`:

- `bias_columns::BiasColumns`: the per-satellite column assignment (see
  [`BiasColumns`](@ref)).
- `clock_systems::Vector{SupportedTimeSystem}`: the time system of each clock column, so
  `clock_systems[k]` belongs to clock column `k`, in order of first appearance. A
  collapsed system has no column of its own; its satellites share the hub's.
- `extra_bands::Vector{Symbol}`: band of each inter-frequency-bias column, so
  `extra_bands[i]` belongs to column `i` of the IFB block.
- `reference_bands::Vector{Symbol}`: per IFB column, the reference band of its coverage
  component — the anchor that column's bias is measured against (see
  [`ifb_layout`](@ref)).
- `hub_system::Union{Nothing,SupportedTimeSystem}`: the system the collapsed
  clocks were merged onto, or `nothing` for a layout that estimates every clock bias
  independently, which is the common case.
- `hub_rows::Vector{PreparedObservation}`: one row per collapsed time system, the
  first of its satellites; each row carries its own `time_system`, so it is its own key.
  A system is collapsed where `time_offsets` holds its offset toward `hub_system` (see
  [`broadcast_time_offset`](@ref)), which [`broadcast_inter_system_bias`](@ref)
  evaluates (see [`calc_hub_range_offsets`](@ref)). Empty exactly when `hub_system` is
  `nothing`.

Only the *classification* fields of each row are read, never its transmit time: the
counts, the coverage graph and which offsets are broadcast are all the decision needs.

The decision is observability-driven, not merely count-driven:

- When the coverage graph is connected and the satellites suffice for the unknowns
  (`n ≥ 3 + num_systems + num_ifbs` observations, of which `3 + num_systems` must come
  from *distinct* satellites — extra bands of an already-tracked satellite add
  inter-frequency-bias information, not geometry), estimate everything independently: the
  inter-system offset and the receiver inter-frequency biases are observed directly from
  the geometry, so neither inherits the broadcast-GGTO error (the satellite group delays
  are already removed per satellite upstream, so the per-band column carries the receiver
  chain).
- Otherwise merge every clock that can be merged onto a *hub* system, using a
  broadcast offset to that hub (see [`broadcast_time_offset`](@ref)). Hubs are tried
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

Because `ifb_layout` never creates an unobservable IFB column, every returned
layout is structurally sound — the degenerate disjoint-band case is removed by
construction, not caught after the fact. The satellite conditions above are structural
too, and thus necessary but not sufficient: a returned layout can still have degenerate
*geometry* (lines of sight that span too little), which no satellite count can see. That
is left to the checks `calc_pvt` makes on the solved geometry — the DOP's positive-definite
test and the velocity solve's own — rather than pre-screened.
"""
function decide_bias_layout(prepared,
    time_offsets::BroadcastTimeOffsets)::Union{Nothing,BiasLayout}
    workspace = BiasLayoutWorkspace()
    decide_bias_layout!(workspace, prepared, time_offsets) ? bias_layout(workspace) :
    nothing
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
    clock_columns::Vector{Int}
    num_clock_biases::Int
    ifb_columns::Vector{Int}
    extra_bands::Vector{Symbol}
    reference_bands::Vector{Symbol}
    hub_system::Union{Nothing,SupportedTimeSystem}
    hub_rows::Vector{PreparedObservation}
    clock_systems::Vector{SupportedTimeSystem}
    # Scratch.
    bands::Vector{Symbol}
    effective_systems::Vector{SupportedTimeSystem}
    ifb_scratch::IFBLayoutScratch{Symbol}
end

BiasLayoutWorkspace() = BiasLayoutWorkspace(Int[], 0, Int[], Symbol[], Symbol[], nothing,
    PreparedObservation[], SupportedTimeSystem[], Symbol[], SupportedTimeSystem[],
    IFBLayoutScratch{Symbol}())

"""
    bias_layout(workspace::BiasLayoutWorkspace) -> BiasLayout

The layout [`decide_bias_layout!`](@ref) last decided into `workspace`. Its vectors
**are** the workspace's, so they are overwritten by the next decision.
"""
bias_layout(workspace::BiasLayoutWorkspace) = BiasLayout((
    BiasColumns(workspace.clock_columns, workspace.num_clock_biases,
        workspace.ifb_columns, length(workspace.extra_bands)),
    workspace.clock_systems,
    workspace.extra_bands,
    workspace.reference_bands,
    workspace.hub_system,
    workspace.hub_rows,
))

"""
    decide_bias_layout!(workspace::BiasLayoutWorkspace, prepared,
                        time_offsets::BroadcastTimeOffsets) -> Bool

[`decide_bias_layout`](@ref) into `workspace`: `true` and the layout stored there (read it
with [`bias_layout`](@ref)), or `false` when the constellation cannot be solved. The same
decision, allocating nothing once the workspace has grown to the epoch.
"""
function decide_bias_layout!(workspace::BiasLayoutWorkspace, prepared,
    time_offsets::BroadcastTimeOffsets)
    num_observations = length(prepared)
    bands = resize!(workspace.bands, num_observations)
    for j in 1:num_observations
        bands[j] = prepared[j].band_id
    end
    # Distinct physical satellites, identified by `(time system, PRN)` — a PRN is only
    # unique within its GNSS. A satellite tracked on several bands appears once per band in
    # `prepared` but supplies one line of sight, so only distinct satellites constrain
    # the geometry and clock unknowns; its repeats constrain the inter-frequency biases.
    num_distinct_sats = count_distinct_satellites(prepared)
    # Both are necessary for a full-rank design (`H` has `3 + M + B` columns, and its rows
    # take only `num_distinct_sats` distinct values outside the IFB columns), neither is
    # sufficient: the geometry itself can still be degenerate, which `calc_pvt` screens
    # for once the design matrix exists.
    enough_satellites() =
        num_observations >=
            3 + workspace.num_clock_biases + length(workspace.extra_bands) &&
        num_distinct_sats >= 3 + workspace.num_clock_biases

    # One row per collapsed system; a row carries its own `time_system`, so it is its
    # own key and no separate dictionary is needed.
    hub_rows = empty!(workspace.hub_rows)
    workspace.hub_system = nothing
    num_components = decide_bias_columns!(workspace, prepared, nothing)
    num_components == 1 && enough_satellites() && return true

    # Connected-but-scarce or disconnected: try collapsing every other system whose
    # offset to a hub system is broadcast onto that hub. Hubs are tried in a fixed
    # order, so the choice is deterministic and GPS wins whenever it can; BDT closes the
    # list only for completeness — no signal broadcasts an offset toward BDT today, so
    # its loop finds nothing.
    for hub_system in CANDIDATE_HUB_SYSTEMS
        any(observation -> observation.time_system === hub_system, prepared) || continue
        empty!(hub_rows)
        for observation in prepared
            observation.time_system === hub_system && continue
            is_collapsed(hub_rows, observation.time_system) && continue
            # By index: the two `SupportedTimeSystem`s would not union-split on Julia 1.10.
            isnothing(broadcast_time_offset_at(time_offsets,
                time_system_index(observation.time_system),
                time_system_index(hub_system))) && continue
            push!(hub_rows, observation)
        end
        isempty(hub_rows) && continue
        decide_bias_columns!(workspace, prepared, hub_system)
        if enough_satellites()
            workspace.hub_system = hub_system
            return true
        end
    end

    # No collapse available. The independent layout is still observable (its IFBs are
    # component-restricted); use it if there are enough satellites, otherwise unsolvable.
    empty!(hub_rows)
    decide_bias_columns!(workspace, prepared, nothing)
    return enough_satellites()
end

# The clock and IFB columns of one candidate layout, into `workspace`: every system
# independent when `hub_system` is `nothing`, else the systems of `workspace.hub_rows`
# merged onto `hub_system`. Returns the number of coverage components.
function decide_bias_columns!(workspace::BiasLayoutWorkspace, prepared, hub_system)
    effective_systems = resize!(workspace.effective_systems, length(prepared))
    for (j, observation) in enumerate(prepared)
        sys = observation.time_system
        effective_systems[j] =
            !isnothing(hub_system) && is_collapsed(workspace.hub_rows, sys) ? hub_system :
            sys
    end
    unique_effective = unique_time_systems!(workspace.clock_systems, effective_systems)
    clock_columns = resize!(workspace.clock_columns, length(prepared))
    # `something`: every system is in `unique_effective` by construction, and saying
    # so keeps the columns a `Vector{Int}` rather than widening them with `Nothing`.
    for (j, sys) in enumerate(effective_systems)
        clock_columns[j] = something(time_system_index(unique_effective, sys))
    end
    workspace.num_clock_biases = length(unique_effective)
    # The clock column, not the time system itself, keys the coverage graph: it
    # separates the constellations identically (and in the same first-appearance
    # order) while being a concretely-typed `Int`, where the `Union`-typed
    # `time_system` field would branch on the system at every comparison inside
    # `ifb_layout!`.
    ifb_layout!(workspace.ifb_columns, workspace.extra_bands,
        workspace.reference_bands, workspace.ifb_scratch, clock_columns,
        workspace.bands)
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

# Whether `system` is one of the systems collapsed onto the hub, i.e. whether
# `hub_rows` holds a representative row for it.
is_collapsed(hub_rows, system) = any(row -> row.time_system === system, hub_rows)

# Distinct physical satellites among `prepared`, identified by `(time system, PRN)`
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
function count_distinct_satellites(prepared)
    distinct = 0
    for j in eachindex(prepared)
        observation = prepared[j]
        repeated = any(firstindex(prepared):(j-1)) do i
            prepared[i].prn == observation.prn &&
                prepared[i].time_system === observation.time_system
        end
        repeated || (distinct += 1)
    end
    distinct
end

"""
    select_reference_system(prepared, hub_system) -> system

The reference time system of an epoch: the one whose clock bias, reference time, week and
start epoch define the reported time ([`PVTSolution`](@ref)'s `reference_system`).

- `prepared`: the epoch's rows.
- `hub_system`: the hub of a collapse ([`BiasLayout`](@ref PositionVelocityTime.BiasLayout)),
  or `nothing` if there is none.

A collapse is anchored on its hub, so the hub is the reference whenever there is one.
Otherwise the reference is the system with the most rows — the best-conditioned reported
time — with ties going to the one that appears first. Throws an `ArgumentError` if there
is neither a hub nor a row. Allocates nothing.
"""
function select_reference_system(prepared, hub_system)
    isnothing(hub_system) || return hub_system
    isempty(prepared) &&
        throw(ArgumentError("no hub and no rows to select a reference system from"))
    # Plain loops rather than `any`/`count`: closures capturing `system`, whose type is
    # the `SupportedTimeSystem` `Union`, allocated on every call.
    reference, most_rows = prepared[firstindex(prepared)].time_system, 0
    for j in eachindex(prepared)
        system = prepared[j].time_system
        # Each system is counted once, at its first appearance; a strict `>` then
        # leaves a tie with the one that appeared first.
        seen = false
        for i in firstindex(prepared):(j-1)
            if prepared[i].time_system === system
                seen = true
                break
            end
        end
        seen && continue
        n = 0
        for observation in prepared
            n += observation.time_system === system
        end
        n > most_rows && ((reference, most_rows) = (system, n))
    end
    reference
end

"""
    broadcast_inter_system_bias(time_offsets::BroadcastTimeOffsets, system, target)
        -> Union{Nothing,Float64}

The inter-system bias (metres) of time system `system` against `target` that the
broadcast offset predicts: `−c · Δt`, with `Δt` the offset from `system` to `target`
([`broadcast_time_offset`](@ref)). Added to a modeled range it puts an observation of
`system` on `target`'s clock, and a collapse onto the hub `target` reports it as
`system`'s entry of `inter_system_biases`. Zero when `system` is `target`, `nothing`
where no offset from `system` toward `target` is broadcast; the reverse direction is not
substituted.
"""
function broadcast_inter_system_bias(
    time_offsets::BroadcastTimeOffsets,
    system::SupportedTimeSystem,
    target::SupportedTimeSystem,
)
    broadcast_inter_system_bias_at(time_offsets, time_system_index(system),
        time_system_index(target))
end

# `broadcast_inter_system_bias` by `time_system_index` (see `broadcast_time_offset_at`).
function broadcast_inter_system_bias_at(time_offsets::BroadcastTimeOffsets, system::Int,
    target::Int)
    system == target && return 0.0
    offset = broadcast_time_offset_at(time_offsets, system, target)
    isnothing(offset) ? nothing : -SPEED_OF_LIGHT * offset
end

"""
    calc_hub_range_offsets(prepared, hub_rows, hub_system, time_offsets) -> Vector{Float64}

Per-satellite range offsets (metres) that carry a clock collapse into the observations,
as decided by [`decide_bias_layout`](@ref): all-zero for the satellites of a system
that keeps its own clock unknown — all of them when `hub_rows` is empty —
and otherwise `−c · Δt_systems` for each satellite of a collapsed system.

The broadcast offset is `Δt_systems = (that system's time) − (the hub system's)` (see
[`broadcast_time_offset`](@ref): the GGTO for Galileo toward a GPS hub, the BGTO for
BeiDou toward either), so a transmit time becomes hub time by SUBTRACTING it; the
modeled range therefore carries `−c·Δt_systems`, and the solve yields
`inter_system_biases[sys] = −c·Δt_systems`. It is one constellation-wide value, taken
from `time_offsets` (see [`BroadcastTimeOffsets`](@ref)), and it converts all of that
system's observations: the offset drifts by ~1e-15 s/s, far below a millimetre across
the spread of the transmit times.
"""
calc_hub_range_offsets(prepared, hub_rows, hub_system, time_offsets) =
    calc_hub_range_offsets!(Vector{Float64}(undef, length(prepared)), prepared, hub_rows,
        hub_system, time_offsets)

# `calc_hub_range_offsets` into `offsets`, resized to the satellite count.
function calc_hub_range_offsets!(offsets, prepared, hub_rows, hub_system, time_offsets)
    offsets = resize!(offsets, length(prepared))
    fill!(offsets, 0.0)
    isnothing(hub_system) && return offsets
    for (j, observation) in enumerate(prepared)
        # `hub_rows` holds at most one row per collapsed system (two in practice), so a
        # linear identity scan beats hashing a key.
        row_index =
            findfirst(row -> row.time_system === observation.time_system, hub_rows)
        isnothing(row_index) && continue
        # `something`: a hub row's system has its offset toward the hub by construction.
        # By index: the two `SupportedTimeSystem`s would not union-split on Julia 1.10.
        offsets[j] = something(broadcast_inter_system_bias_at(time_offsets,
            time_system_index(observation.time_system), time_system_index(hub_system)))
    end
    offsets
end
