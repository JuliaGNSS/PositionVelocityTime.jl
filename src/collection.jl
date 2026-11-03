# The collection pass: signal observations in, flat prepared observations out.
#
# Everything that depends on a satellite's navigation-data type — decoder, ranging signal,
# propagator, clock polynomial, broadcast time offsets — is evaluated here, once per
# satellite, into the one concrete `PreparedObservation` the solver works on. This
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
satellite can be on, and those [`BroadcastTimeOffsets`](@ref) holds the broadcast offsets
between.
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
system added later joins `SUPPORTED_TIME_SYSTEMS`, and this follows; `BroadcastTimeOffsets`
then grows a row and a column for it.
"""
const SupportedTimeSystem = Union{map(typeof, SUPPORTED_TIME_SYSTEMS)...}

"""
    PreparedObservation

One signal observation's contribution to one epoch: the concrete, parameter-free record the solver
works on. [`collect_measurements`](@ref) has already evaluated everything that depends on
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

The broadcast time offsets a satellite carries belong to a pair of time systems, not to
it; [`collect_measurements`](@ref) gathers them into one [`BroadcastTimeOffsets`](@ref).
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
function satellite_measurement(observation::SignalObservation, approximate_year::Integer)
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
    collect_measurements(observations; approximate_year = year(now(UTC)))
        -> (Vector{PreparedObservation}, ionospheric_parameters, time_offsets)

Prepare the observation collections of one epoch (see [`calc_pvt`](@ref)) into the
solver's input: the [`PreparedObservation`](@ref)s, the ionospheric parameters
[`select_ionospheric_correction`](@ref) would pick from the same satellites, and the
[`BroadcastTimeOffsets`](@ref) between their time systems.

A satellite is kept only if its navigation data is complete (checked first, which
guarantees its health bit is decoded) and it reports healthy. The prepared observations
are in collection order × within-collection order (see [`calc_pvt`](@ref), also for
`approximate_year`).
"""
function collect_measurements(observations; approximate_year::Integer = year(now(UTC)))
    collections = observation_collections(observations)
    prepared = PreparedObservation[]
    # Reserve the unfiltered count (a prepared observation is 120 bytes) rather than regrow.
    sizehint!(prepared, sum(length, collections; init = 0))
    (; ionospheric_model, time_offsets) =
        _collect_measurements!(prepared, collections, approximate_year)
    prepared, ionospheric_model.parameters, time_offsets
end

"""
    collect_measurements!(prepared::Vector{PreparedObservation}, observations;
                          approximate_year = year(now(UTC)))
        -> (prepared, ionospheric_parameters, time_offsets)

[`collect_measurements`](@ref) into an existing vector: `prepared` is emptied and
refilled with this epoch's observations. What it held is overwritten.

The returned tuple pairs the (non-`isbits`) vector with a `Union` of coefficient sets,
which Julia returns boxed: one small allocation per call. A consumer that must not
allocate calls [`collect_measurement_rows!`](@ref).
"""
function collect_measurements!(
    prepared::Vector{PreparedObservation},
    observations;
    approximate_year::Integer = year(now(UTC)),
)
    (; ionospheric_model, time_offsets) =
        collect_measurement_rows!(prepared, observations; approximate_year)
    prepared, ionospheric_model.parameters, time_offsets
end

"""
    collect_measurement_rows!(prepared::Vector{PreparedObservation}, observations;
                              approximate_year = year(now(UTC)))
        -> (; ionospheric_model::IonosphericModel, time_offsets::BroadcastTimeOffsets)

[`collect_measurements!`](@ref) without the vector in the result: `prepared` is
emptied and refilled with this epoch's observations, and the ionospheric parameters,
wrapped in a concrete [`IonosphericModel`](@ref), and the [`BroadcastTimeOffsets`](@ref)
are returned. Both are concrete, so Julia returns them unboxed, and once `prepared` has
held as many observations before, this allocates nothing. It is the collection pass
[`calc_pvt!`](@ref) runs, for a consumer that runs its own estimator over the same
observations in a loop that must not allocate.
"""
collect_measurement_rows!(
    prepared::Vector{PreparedObservation},
    observations;
    approximate_year::Integer = year(now(UTC)),
) = _collect_measurements!(prepared, observation_collections(observations),
    approximate_year)

# The collection pass proper, on validated collections. Positional for the solver.
function _collect_measurements!(prepared, collections::Tuple, approximate_year::Integer)
    empty!(prepared)
    # In collection order. Over the `Tuple`, not the `NamedTuple`, whose `map` would not
    # specialise on the closure and dispatch each collection dynamically.
    per_collection = map(collections) do collection
        collect_group!(prepared, collection, approximate_year)
    end
    parameters = select_from_ionospheric_candidates(
        merge_all_ionospheric_candidates(map(first, per_collection)))
    # Collection by collection, so a tie keeps the first observation's offset.
    time_offsets =
        foldl(merge_time_offsets, map(last, per_collection); init = BroadcastTimeOffsets())
    (; ionospheric_model = IonosphericModel(parameters), time_offsets)
end

# One collection's observations, appended to `prepared`, its ionospheric candidates
# and its broadcast time offsets. Specialises on the collection's `SignalObservation`
# type, so every call in the loop is static.
function collect_group!(prepared, collection::ObservationCollection,
    approximate_year::Integer)
    candidates = NO_IONOSPHERIC_CANDIDATES
    time_offsets = BroadcastTimeOffsets()
    for observation in collection
        has_positioning_data(observation.decoder) || continue
        is_sat_healthy(observation.decoder) || continue
        prepared_observation = satellite_measurement(observation, approximate_year)
        push!(prepared, prepared_observation)
        candidates = update_ionospheric_candidates(candidates, observation.decoder)
        time_offsets = update_time_offsets(time_offsets, observation.decoder,
            get_time_system(observation.signal), prepared_observation.transmit_time)
    end
    candidates, time_offsets
end
