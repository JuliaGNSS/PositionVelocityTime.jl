# The receiver bias columns of the least-squares state: one clock per time system, the
# inter-frequency biases, the broadcast-offset constraints between clocks, and the
# reference system.

"""
    ClockConstraint

One known difference between two receiver clock biases, fitted as an extra observation
row: `ξ[3 + column] − ξ[3 + target_column] = value`. It is how a broadcast GNSS time
offset (GGTO, BGTO) ties the clocks of two time systems together; see
[`decide_bias_layout`](@ref).

# Fields
- `column::Int`: the clock column of the system whose broadcast offset is used.
- `target_column::Int`: the clock column of the system the offset is toward.
- `value::Float64`: the known difference in metres
  ([`broadcast_inter_system_bias`](@ref)), which the fit meets.
"""
struct ClockConstraint
    column::Int
    target_column::Int
    value::Float64
end

"""
    BiasColumns

Per-satellite assignment of the bias columns of the least-squares state
`[x, y, z, tc₁, …, tc_num_clock_biases, ifb₁, …, ifb_num_ifbs]`
([`num_parameters`](@ref) entries), shared by [`predict_pseudoranges!`](@ref),
[`calc_design_matrix!`](@ref) and [`estimate_position`](@ref). A clock column is the
receiver clock of one GNSS time system; an inter-frequency-bias column is the receiver's
per-band RF-chain delay.

# Fields
- `clock_columns::Vector{Int}`: per satellite, the clock column (1…`num_clock_biases`)
  of its GNSS time system; the design-matrix `1.0` lands at `3 + clock_columns[j]`.
- `num_clock_biases::Int`: number of clock biases (also the offset of the IFB block).
- `ifb_columns::Vector{Int}`: per satellite, the inter-frequency-bias
  column (1…`num_ifbs`) of its frequency band, or `0` for the reference
  band; the `1.0` lands at `3 + num_clock_biases + ifb_columns[j]`.
- `num_ifbs::Int`: number of inter-frequency biases (frequency bands
  beyond the reference).
- `constraints::Vector{ClockConstraint}`: the clock-difference rows, appended after the
  satellite rows. Without the argument, a shared empty vector that must not be modified.
"""
struct BiasColumns
    clock_columns::Vector{Int}
    num_clock_biases::Int
    ifb_columns::Vector{Int}
    num_ifbs::Int
    constraints::Vector{ClockConstraint}
end

BiasColumns(clock_columns, num_clock_biases, ifb_columns, num_ifbs) = BiasColumns(
    clock_columns, num_clock_biases, ifb_columns, num_ifbs, NO_CLOCK_CONSTRAINTS)

# Shared by every `BiasColumns` without constraints, so building one allocates nothing.
# Never pushed to.
const NO_CLOCK_CONSTRAINTS = ClockConstraint[]

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
graph. `constellation_keys[j]` and `bands[j]` are satellite `j`'s constellation and band;
any key that separates the constellations will do — [`decide_bias_layout`](@ref) passes
each satellite's clock group.

Two bands share a *coverage component* iff some constellation is tracked on both
(directly or through a chain of shared constellations). Within a component the clock/IFB
split has one gauge freedom, so its most-populated band (ties by first appearance) is the
reference and only its other bands get an IFB column. A band alone in its component gets
none: its delay folds into its constellation's clock. So exactly the observable IFBs are
estimated.

`ifb_indices[j]` is satellite `j`'s IFB column, or `0` on a reference
band; `extra_bands[i]` is the band of IFB column `i` and `reference_bands[i]` the
reference band it is measured against; `num_components` is the number of coverage
components (`1` ⇔ connected).
"""
function ifb_layout(constellation_keys, bands)
    ifb_indices = Int[]
    extra_bands = Vector{eltype(bands)}()
    reference_bands = Vector{eltype(bands)}()
    num_components = ifb_layout!(ifb_indices,
        extra_bands, reference_bands, IFBLayoutScratch{eltype(bands)}(),
        constellation_keys, bands)
    return ifb_indices, extra_bands, reference_bands, num_components
end

# The scratch of `ifb_layout!`. Bands are identified by their position in
# `unique_bands`, so the bookkeeping is plain `Int` vectors rather than a `Dict` per epoch.
struct IFBLayoutScratch{B}
    unique_bands::Vector{B}
    band_index::Vector{Int}   # per satellite: its band's position in `unique_bands`
    parent::Vector{Int}       # per band: union-find parent
    band_count::Vector{Int}   # per band: satellites on it
    reference::Vector{Int}    # per union-find root: the component's reference band, or 0
    column::Vector{Int}       # per band: its IFB column, or 0 for a reference band
end

IFBLayoutScratch{B}() where {B} = IFBLayoutScratch{B}(B[], Int[], Int[], Int[], Int[], Int[])

# `ifb_layout` into the given vectors (emptied first), working in
# `scratch`; returns `num_components`. Allocates nothing once its buffers have grown.
# Internal: `decide_bias_layout!` is the public allocation-free way in.
function ifb_layout!(
    ifb_indices,
    extra_bands,
    reference_bands,
    scratch::IFBLayoutScratch,
    constellation_keys,
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
    # Union the bands a constellation is tracked on: link each satellite's band to that
    # of its constellation's first satellite.
    resize!(parent, num_bands)
    parent .= 1:num_bands
    for j in 1:num_sats
        first_of_constellation = something(
            findfirst(i -> isequal(constellation_keys[i], constellation_keys[j]), 1:j))
        a = find_root!(parent, band_index[first_of_constellation])
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

# The union-find root of `b` in `parent`, compressing the path. A loop, as a recursive
# closure would be boxed and dispatched dynamically.
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

The least-squares bias layout of one epoch, as [`decide_bias_layout!`](@ref) decides it,
refilled in place on reuse.

# Fields
- `clock_columns::Vector{Int}`: per satellite, the clock column of its time system.
- `clock_systems::Vector{SupportedTimeSystem}`: the time system of each clock column, in
  order of first appearance.
- `constraints::Vector{ClockConstraint}`: the broadcast-offset rows; usually empty.
- `clock_group_columns::Vector{Int}`: per satellite, its *clock group*, in which the systems
  `constraints` tie together count as one.
- `num_clock_groups::Int`: the number of clock groups.
- `ifb_columns::Vector{Int}`, `extra_bands::Vector{Symbol}`,
  `reference_bands::Vector{Symbol}`: as [`ifb_layout`](@ref) lays them out
  over the clock groups.

[`fit_bias_columns`](@ref) is the design the position solve fits, [`dop_bias_columns`](@ref)
the one its DOP is taken from.
"""
mutable struct BiasLayout
    clock_columns::Vector{Int}
    clock_systems::Vector{SupportedTimeSystem}
    constraints::Vector{ClockConstraint}
    clock_group_columns::Vector{Int}
    num_clock_groups::Int
    ifb_columns::Vector{Int}
    extra_bands::Vector{Symbol}
    reference_bands::Vector{Symbol}
    # Scratch.
    bands::Vector{Symbol}
    sat_counts::Vector{Int}                  # per clock column: distinct satellites
    candidates::Vector{ClockConstraint}      # one per pair of systems with an offset
    system_parent::Vector{Int}               # per clock column: union-find parent
    clock_group_of_root::Vector{Int}         # per union-find root: its group, or 0
    ifb_scratch::IFBLayoutScratch{Symbol}
end

BiasLayout() = BiasLayout(Int[], SupportedTimeSystem[], ClockConstraint[], Int[], 0, Int[],
    Symbol[], Symbol[], Symbol[], Int[], ClockConstraint[], Int[], Int[],
    IFBLayoutScratch{Symbol}())

"""
    fit_bias_columns(layout::BiasLayout) -> BiasColumns

The [`BiasColumns`](@ref) the position solve fits: one clock column per time system, and
the clock constraints as extra rows. Its vectors **are** the layout's.
"""
fit_bias_columns(layout::BiasLayout) = BiasColumns(layout.clock_columns,
    length(layout.clock_systems), layout.ifb_columns,
    length(layout.extra_bands), layout.constraints)

"""
    dop_bias_columns(layout::BiasLayout) -> BiasColumns

The [`BiasColumns`](@ref) the DOP is taken from: one clock column per clock group and no
constraint rows, i.e. the constraints taken as exact. Its vectors **are** the layout's;
its `constraints` is a shared empty vector that must not be modified.
"""
dop_bias_columns(layout::BiasLayout) = BiasColumns(layout.clock_group_columns,
    layout.num_clock_groups, layout.ifb_columns,
    length(layout.extra_bands), NO_CLOCK_CONSTRAINTS)

"""
    decide_bias_layout(measurements) -> Union{BiasLayout,Nothing}

Decide the least-squares bias layout (see [`BiasLayout`](@ref)) for the
[`SatelliteMeasurement`](@ref)s of one epoch, or `nothing` when the epoch cannot be
solved.

Every time system present gets its own clock column. A broadcast time offset between two
of them — Galileo's GGTO toward GPS Time, GPS's CNAV/CNAV-2 GGTO toward Galileo System
Time, BeiDou's BGTO toward either (see [`has_broadcast_time_offset`](@ref)) — can tie them as
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
   used (none, if the satellites already suffice). A band stranded on its own
   constellation then gets no IFB column, so the reported inter-system bias reads
   `GGTO + (IFB difference)`; position, residuals and DOP are unaffected.

Failing all three, the epoch is unsolvable. Each constraint pins an otherwise
unobservable direction, so the fit meets it exactly and the constrained system's
inter-system bias reads as the broadcast `−c·Δt`. Galileo's GGTO (specified to 20 ns,
95 %) is preferred over GPS's. A constraint is evaluated at the transmit time of the first
satellite carrying it ([`broadcast_inter_system_bias`](@ref)).

The conditions are necessary but not sufficient: a degenerate geometry is caught by
`calc_pvt` from the DOP of the solved design.
"""
function decide_bias_layout(measurements)::Union{Nothing,BiasLayout}
    layout = BiasLayout()
    decide_bias_layout!(layout, measurements) ? layout : nothing
end

"""
    decide_bias_layout!(layout::BiasLayout, measurements) -> Bool

[`decide_bias_layout`](@ref) into `layout`: `true` and the layout stored there, or `false`
when the epoch cannot be solved. Allocates nothing once `layout` has grown to the epoch.

After `false`, `layout` is still a consistent, unsolvable layout: it ties every pair of
systems a broadcast offset can tie, through the most preferred offsets (a spanning forest
over the systems), with the groups and IFB columns that go with them — for a navigation
filter that keeps its clocks tied through an epoch it cannot solve.
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
    clock_columns = resize!(layout.clock_columns, num_sats)
    sat_counts = fill!(resize!(layout.sat_counts, num_systems), 0)
    for j in 1:num_sats
        # `something` keeps the columns a `Vector{Int}`; the system is always found.
        system = measurements[j].time_system
        column = something(findfirst(other -> other === system, clock_systems))
        clock_columns[j] = column
        is_new_sat(measurements, j) && (sat_counts[column] += 1)
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
        min(sat_counts[candidate.column], sat_counts[candidate.target_column]),
        -max(sat_counts[candidate.column], sat_counts[candidate.target_column]),
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
    # Unsolvable: the spanning offsets, as `decide_bias_layout!` documents.
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
    clock_group_columns = resize!(layout.clock_group_columns, length(measurements))
    num_clock_groups = 0
    for j in eachindex(clock_group_columns)
        root = find_root!(parent, layout.clock_columns[j])
        if clock_group_of_root[root] == 0
            num_clock_groups += 1
            clock_group_of_root[root] = num_clock_groups
        end
        clock_group_columns[j] = clock_group_of_root[root]
    end
    layout.num_clock_groups = num_clock_groups
    # Keyed by clock group, so the tied systems count as one constellation.
    num_components = ifb_layout!(layout.ifb_columns,
        layout.extra_bands, layout.reference_bands, layout.ifb_scratch,
        clock_group_columns, layout.bands)
    # Necessary, not sufficient (see `decide_bias_layout`).
    solvable =
        length(measurements) >= 3 + num_clock_groups + length(layout.extra_bands) &&
        sum(layout.sat_counts) >= 3 + num_clock_groups
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
# measurement of the former's system that carries that offset; `false` if none does.
function push_offset_candidate!(candidates, measurements, clock_systems, column,
    target_column)
    system, target = clock_systems[column], clock_systems[target_column]
    j = first_offset_carrier(measurements, system, target)
    j == 0 && return false
    value = something(broadcast_inter_system_bias(measurements[j], target))
    push!(candidates, ClockConstraint(column, target_column, value))
    true
end

# The index of the first measurement of `system` that carries a broadcast offset toward
# `target`, or 0.
function first_offset_carrier(measurements, system, target)
    for j in eachindex(measurements)
        measurements[j].time_system === system &&
            has_broadcast_time_offset(measurements[j], target) && return j
    end
    0
end

# Identity-based helpers over the `time_system` field: every `TimeSystem` is a singleton,
# so `===` is exact and compiles to a pointer test, where `==` or a `Dict` lookup would
# first branch on the field's `Union`.

# The distinct time systems of `systems` into `unique_systems` (emptied first), in order
# of first appearance — the order that fixes the clock columns of `BiasColumns`.
function unique_time_systems!(unique_systems, systems)
    empty!(unique_systems)
    for system in systems
        any(other -> other === system, unique_systems) || push!(unique_systems, system)
    end
    unique_systems
end

# Whether measurement `j` is the first of its physical satellite, `(time system, PRN)`: a
# satellite tracked on several bands is several measurements but one line of sight. A scan
# of the earlier measurements, as a set of `(SupportedTimeSystem, Int)` keys would box each
# one; quadratic, which for an epoch's dozens of measurements is the cheaper.
function is_new_sat(measurements, j)
    measurement = measurements[j]
    for i in firstindex(measurements):(j-1)
        measurements[i].prn == measurement.prn &&
            measurements[i].time_system === measurement.time_system && return false
    end
    true
end

"""
    select_reference_system(layout::BiasLayout) -> system

The reference time system of the epoch `layout` was decided for, solvable or not: the one
that dates the reported time ([`PVTSolution`](@ref)'s `reference_system`). It is the
system with the most distinct satellites, ties going to the first. Throws an
`ArgumentError` for a layout without a time system. Allocates nothing.

Distinct satellites, not measurements: a satellite's extra bands mostly determine
inter-frequency biases. It is also the count `decide_bias_layout!` constrains by, so the
reference system is never one whose clock only a broadcast offset determines.
"""
function select_reference_system(layout::BiasLayout)
    isempty(layout.clock_systems) &&
        throw(ArgumentError("no time system to select a reference system from"))
    layout.clock_systems[select_reference_clock_column(layout)]
end

# The clock column of `select_reference_system`; `argmax` takes the first maximum.
select_reference_clock_column(layout::BiasLayout) = argmax(layout.sat_counts)

"""
    broadcast_inter_system_bias(measurement, target) -> Union{Nothing,Float64}

The inter-system bias (metres) of `measurement`'s time system against `target` that the
broadcast offset predicts at its transmit time: `−c · Δt_systems`, with `Δt_systems` the
broadcast offset ([`broadcast_time_offset`](@ref)). It is the clock difference
`tc_own − tc_target` a [`ClockConstraint`](@ref) fixes, signed like [`PVTSolution`](@ref)'s
`inter_system_biases`. Zero when `measurement` is of `target`, `nothing` where it carries
no usable offset toward it ([`has_broadcast_time_offset`](@ref)).
"""
function broadcast_inter_system_bias(
    measurement::SatelliteMeasurement,
    target::SupportedTimeSystem,
)
    measurement.time_system === target && return 0.0
    offset = broadcast_time_offset(measurement, target)
    isnothing(offset) ? nothing : -SPEED_OF_LIGHT * offset
end
