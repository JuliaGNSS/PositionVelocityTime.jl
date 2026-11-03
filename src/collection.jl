# The preparation pass: signal observations in, flat prepared observations out.
#
# Everything that depends on a satellite's navigation-data type — decoder, ranging signal,
# propagator, clock polynomial, broadcast inter-system time offsets — is evaluated here,
# once per satellite, into the one concrete `PreparedObservation` the solver works on. This
# is `calc_pvt`'s function barrier: upstream specialises on the collection shape (and is
# small), downstream compiles once for every constellation mix.

"""
    SignalObservation{CP<:Real,D<:GNSSDecoderState,S<:AbstractGNSSSignal}

One satellite's decoder state together with its code phase and carrier Doppler.

# Fields
- `decoder::GNSSDecoderState`: the decoded navigation data.
- `signal::AbstractGNSSSignal`: the ranging signal the code phase and Doppler were
  measured on (e.g. `GPSL1CA()`, `GalileoE1C()`); may be a dataless pilot whose data
  `decoder` decoded from the band's data component.
- `code_phase::CP`: code phase in chips.
- `carrier_doppler`: carrier Doppler in Hz.

# Constructors
    SignalObservation(; decoder, signal, code_phase, carrier_doppler)

To read a `Tracking.TrackedSat`, pass its `get_code_phase` and `get_carrier_doppler`
results to the keyword constructor. A carrier phase does not enter: the code alone
times the transmission.
"""
@kwdef struct SignalObservation{CP<:Real,D<:GNSSDecoder.GNSSDecoderState,S<:AbstractGNSSSignal}
    decoder::D
    signal::S
    code_phase::CP
    carrier_doppler::typeof(1.0Hz)
end

# A collection retains one concrete observation type, even when it is empty. A tuple
# keeps different decoder/signal types separate until the concrete preparation barrier.
const ObservationCollection =
    Union{AbstractVector{<:SignalObservation},Dictionary{<:Any,<:SignalObservation}}

function check_observation_collection(observations::ObservationCollection)
    isconcretetype(eltype(observations)) || throw_observation_collections()
    observations
end

observation_collections(observations::ObservationCollection) =
    (check_observation_collection(observations),)
observation_collections(observations::Tuple{Vararg{ObservationCollection}}) =
    map(check_observation_collection, observations)
observation_collections(
    observations::NamedTuple{<:Any,<:Tuple{Vararg{ObservationCollection}}},
) =
    observation_collections(values(observations))
observation_collections(observations) = throw_observation_collections()

function throw_observation_collections()
    throw(
        ArgumentError(
            "calc_pvt takes a concrete-typed collection of SignalObservation, or a tuple " *
            "or NamedTuple of such collections: calc_pvt((gps = gps_observations, " *
            "galileo = galileo_observations)). Keep different observation types in separate " *
            "collections; an empty collection must retain its concrete element type.",
        ),
    )
end

"""
    PositionVelocityTime.SUPPORTED_TIME_SYSTEMS

The GNSS time systems this package solves with, `(GPST(), GST(), BDT())`: those a
satellite can be on, and those [`BroadcastInterSystemTimeOffsets`](@ref) holds the
broadcast offsets between.
"""
const SUPPORTED_TIME_SYSTEMS = (GPST(), GST(), BDT())

"""
    PositionVelocityTime.SupportedTimeSystem

`Union{GPST,GST,BDT}`: the GNSS time systems a satellite can be on, which are exactly
[`SUPPORTED_TIME_SYSTEMS`](@ref). The type of [`PreparedObservation`](@ref)'s
`time_system` field, and of the keys and the reference system of a [`PVTSolution`](@ref).

A closed `Union` rather than the abstract `GNSSSignals.TimeSystem`, so every call on a
time system splits into static branches instead of a dynamic dispatch, which a
`juliac --trim` build rejects, and a `PreparedObservation` stays one concrete type. A time
system added later joins `SUPPORTED_TIME_SYSTEMS`, and this follows;
`BroadcastInterSystemTimeOffsets` then grows a row and a column for it.
"""
const SupportedTimeSystem = Union{map(typeof, SUPPORTED_TIME_SYSTEMS)...}

"""
    PreparedObservation

One signal observation's contribution to one epoch: the concrete, parameter-free record the solver
works on. [`prepare_observations`](@ref) has already evaluated everything that depends on
the message family, the decoder or the signal, so a mixed-constellation epoch is one
`Vector{PreparedObservation}`.

# Fields
- `prn::Int`: the satellite's PRN, unique only within its constellation.
- `signal_id::Symbol`: `get_signal_id` of the ranging signal (`:GPSL1CA`, …), half of the
  `(signal, PRN)` key of `PVTSolution.observations`.
- `band_id::Symbol`: `get_band_id` of the ranging signal (`:L1`, `:L5`, …), which groups
  the inter-frequency biases.
- `time_system::SupportedTimeSystem`: `get_time_system` of the ranging signal, which
  groups the receiver clock biases.
- `transmit_time::Float64`: the transmit time ([`calc_transmit_time`](@ref)) as a
  seconds-of-week count **on this satellite's own system's scale**.
- `position::SVector{3,Float64}`: satellite ECEF position at `transmit_time` (m).
- `velocity::SVector{3,Float64}`: satellite ECEF velocity at `transmit_time` (m/s).
- `carrier_doppler::Float64`: the measured carrier Doppler (Hz, unitless).
- `center_frequency::Float64`: the ranging signal's carrier (Hz, unitless) — the
  Doppler wavelength, and the `1/f²` of every ionospheric model.
- `satellite_clock_drift::Float64`: [`calc_satellite_clock_drift`](@ref) at
  `transmit_time` (s/s).
- `week::Int`: the absolute week number of `time_system`, with the GPS L1 C/A
  1024-week rollover already resolved against `approximate_year`; its start is
  GNSSSignals' `get_week_start_epoch(observation.time_system, observation.week)`.

The broadcast inter-system time offsets a satellite carries belong to a pair of time
systems, not to it; [`prepare_observations`](@ref) gathers them into one
[`BroadcastInterSystemTimeOffsets`](@ref).
"""
struct PreparedObservation
    prn::Int
    signal_id::Symbol
    band_id::Symbol
    time_system::SupportedTimeSystem
    transmit_time::Float64
    position::SVector{3,Float64}
    velocity::SVector{3,Float64}
    carrier_doppler::Float64
    center_frequency::Float64
    satellite_clock_drift::Float64
    week::Int
end

# One signal observation → one prepared observation: every decoder- and signal-typed call
# of the solve.
function prepare_observation(observation::SignalObservation, approximate_year::Integer)
    decoder = observation.decoder
    signal = observation.signal
    ephemeris = something(get_broadcast_ephemeris(decoder))
    transmit_time = correct_clock(ephemeris, group_delay(decoder, signal),
        calc_uncorrected_transmit_time(observation))
    pos_and_vel = calc_satellite_position_and_velocity(ephemeris, transmit_time)
    PreparedObservation(
        decoder.prn,
        get_signal_id(signal),
        get_band_id(signal),
        get_time_system(signal),
        transmit_time,
        pos_and_vel.position,
        pos_and_vel.velocity,
        ustrip(Hz, observation.carrier_doppler),
        ustrip(Hz, get_center_frequency(signal)),
        calc_satellite_clock_drift(ephemeris, transmit_time),
        something(get_broadcast_week(decoder; approximate_year)),
    )
end

"""
    prepare_observations(observations; approximate_year = year(now(UTC)))
        -> (; prepared::Vector{PreparedObservation},
              ionospheric_model::IonosphericModel,
              inter_system_time_offsets::BroadcastInterSystemTimeOffsets)

Prepare the observation collections of one epoch (see [`calc_pvt`](@ref)) into the solver's
input: the [`PreparedObservation`](@ref)s, and, gathered from the same observations, the
[`IonosphericModel`](@ref) and the [`BroadcastInterSystemTimeOffsets`](@ref) between
their time systems. A `NamedTuple`, to destructure by name:

```julia
(; prepared, ionospheric_model, inter_system_time_offsets) =
    prepare_observations(observations)
```

A satellite is kept only if its decoder is ready for positioning (GNSSDecoder's
`has_positioning_data`: the message set is complete and its time is decoded; checked
first, which guarantees the health bit is decoded) and it reports healthy. The
prepared observations are in collection order × within-collection order (see
[`calc_pvt`](@ref), also for `approximate_year`).
A `(signal_id, PRN)` pair kept twice, such as `GalileoE1B()` and `GalileoE1B(CBOC)` of
one satellite, throws an `ArgumentError`.
"""
function prepare_observations(observations; approximate_year::Integer = year(now(UTC)))
    collections = observation_collections(observations)
    prepared = PreparedObservation[]
    # Reserve the unfiltered count (a prepared observation is 120 bytes) rather than regrow.
    sizehint!(prepared, sum(length, collections; init = 0))
    (; ionospheric_model, inter_system_time_offsets) =
        _prepare_observations!(prepared, collections, approximate_year)
    (; prepared, ionospheric_model, inter_system_time_offsets)
end

"""
    prepare_observations!(prepared::Vector{PreparedObservation}, observations;
                          approximate_year = year(now(UTC)))
        -> (; ionospheric_model::IonosphericModel,
              inter_system_time_offsets::BroadcastInterSystemTimeOffsets)

[`prepare_observations`](@ref) into `prepared` (emptied first), returning the
concrete [`IonosphericModel`](@ref) and the [`BroadcastInterSystemTimeOffsets`](@ref) for
the measurement-model functions and [`decide_bias_layout!`](@ref). Allocates nothing once
its buffers have grown. It is the pass [`calc_pvt!`](@ref) runs.
"""
prepare_observations!(
    prepared::Vector{PreparedObservation},
    observations;
    approximate_year::Integer = year(now(UTC)),
) = _prepare_observations!(prepared, observation_collections(observations),
    approximate_year)

# The preparation pass proper, on validated collections. Positional for the solver.
function _prepare_observations!(prepared, collections::Tuple, approximate_year::Integer)
    empty!(prepared)
    # In collection order. Over the `Tuple`, not the `NamedTuple`, whose `map` would not
    # specialise on the closure and dispatch each collection dynamically.
    per_collection = map(collections) do collection
        prepare_collection!(prepared, collection, approximate_year)
    end
    parameters = select_from_ionospheric_candidates(
        merge_all_ionospheric_candidates(map(first, per_collection)))
    # Collection by collection, so a tie keeps the first observation's offset.
    inter_system_time_offsets = foldl(merge_inter_system_time_offsets,
        map(last, per_collection); init = BroadcastInterSystemTimeOffsets())
    check_observation_keys(prepared)
    (; ionospheric_model = IonosphericModel(parameters), inter_system_time_offsets)
end

# A broadcast component and PRN identify one selected ranging observation, regardless
# of replica variant. Checked on the kept observations, before any estimator sees them,
# so that none double-counts one.
function check_observation_keys(observations)
    for j = 2:length(observations), k = 1:j-1
        observations[j].signal_id === observations[k].signal_id &&
            observations[j].prn == observations[k].prn &&
            throw(
                ArgumentError(
                    "a (signal, PRN) pair appears twice in one epoch; select one observation per broadcast component",
                ),
            )
    end
    nothing
end

# One collection's observations, appended to `prepared`, its ionospheric candidates
# and its broadcast inter-system time offsets. Specialises on the collection's
# `SignalObservation` type, so every call in the loop is static.
function prepare_collection!(prepared, collection::ObservationCollection,
    approximate_year::Integer)
    candidates = NO_IONOSPHERIC_CANDIDATES
    inter_system_time_offsets = BroadcastInterSystemTimeOffsets()
    for observation in collection
        has_positioning_data(observation.decoder) || continue
        is_sat_healthy(observation.decoder) || continue
        prepared_observation = prepare_observation(observation, approximate_year)
        push!(prepared, prepared_observation)
        candidates = update_ionospheric_candidates(candidates, observation.decoder)
        inter_system_time_offsets = update_inter_system_time_offsets(
            inter_system_time_offsets, observation.decoder,
            get_time_system(observation.signal), prepared_observation.transmit_time)
    end
    candidates, inter_system_time_offsets
end
