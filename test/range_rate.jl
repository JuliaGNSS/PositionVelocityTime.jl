# The modelled range rate must be the exact time derivative of the pseudorange the
# position solve models — Earth rotation over the travel time included — or the velocity
# solve inherits millimetres per second of bias from the model alone. The test below
# derives its Dopplers from that pseudorange directly, by finite difference, so it shares
# none of the solver's algebra.

using PositionVelocityTime: SVector

const ONE_CLOCK_COLUMN = PositionVelocityTime.BiasColumns([1], 1, [0], 0)

# The pseudorange at receiver time `t`, as `calc_ρ_hat!` models it: the satellite taken
# at its transmit time `t_t = t − ρ/c` (found by fixed-point iteration), the receiver at
# `user_at(t)` with clock bias `bias_at(t)` (m), and the satellite clock offset
# `sat_clock_at(t_t)` (m) subtracted.
function modelled_pseudorange(sat_at, sat_clock_at, user_at, bias_at, t)
    c = PositionVelocityTime.SPEED_OF_LIGHT
    ρ = 0.0
    t_t = t
    for _ in 1:10
        t_t = t - ρ / c
        ρ = PositionVelocityTime.calc_ρ_hat!(
            zeros(1), [sat_at(t_t)], [user_at(t); bias_at(t)], ONE_CLOCK_COLUMN)[1] -
            sat_clock_at(t_t)
    end
    return ρ, t_t
end

central_difference(f, t, h) = (f(t + h) - f(t - h)) / 2h

# `state` re-timed so that its corrected transmit time is `transmit_time`, carrying
# `doppler` (Hz). Only the code phase and bit count move, as in `_precompile_states`;
# the carrier phase, which also counts into the transmit time, is kept.
function retimed_state(state, transmit_time, doppler)
    uncorrected = PositionVelocityTime.calc_uncorrected_time(state)
    # `correct_clock` is a slowly varying offset, so a few fixed-point steps invert it.
    for _ in 1:5
        uncorrected += transmit_time -
                       PositionVelocityTime.correct_clock(state.decoder, state.system, uncorrected)
    end
    rate = ustrip(Hz, get_data_frequency(state.decoder))
    chips = ustrip(Hz, get_code_frequency(state.system))
    offset = state.decoder.num_bits_after_valid_syncro_sequence / rate +
             mod(state.code_phase, chips / rate) / chips +
             uncorrected - PositionVelocityTime.calc_uncorrected_time(state)
    num_bits = floor(Int, offset * rate)
    SatelliteState(;
        decoder = GNSSDecoder.GNSSDecoderState(state.decoder;
            num_bits_after_valid_syncro_sequence = num_bits),
        system = state.system,
        code_phase = (offset - num_bits / rate) * chips,
        carrier_doppler = doppler * Hz,
        carrier_phase = state.carrier_phase,
    )
end

@testset "calc_pvt recovers the true velocity from exact observables: $name" for (name, states) in (
    ("GPS L1 C/A", gps_l1_states(0.0Hz)),
    ("Galileo E1B", galileo_e1b_states(0.0Hz)),
)
    c = PositionVelocityTime.SPEED_OF_LIGHT
    kwargs = (; approximate_year = 2021, enable_ionospheric_correction = false,
        enable_tropospheric_correction = false)
    # The fixture's own fix sets a realistic geometry, and its latest transmit time plus
    # a typical travel time the receiver clock reading at reception.
    fixture = calc_pvt(signal_group(states); kwargs...)
    user_pos = SVector{3,Float64}(fixture.position)
    # An aircraft's speed: the travel-time rate's user-velocity share of the Earth-rotation
    # term is ~2e-6·|v| per satellite, and leaving it out of the design costs ~0.1 mm/s of
    # velocity here — ten times the bounds below, which the solve meets to ~1e-6 m/s.
    user_vel = SVector(150.0, -200.0, 50.0)
    clock_bias, clock_drift = 0.0, 30.0
    reception = maximum(PositionVelocityTime.calc_corrected_time, states) + 0.075

    exact_states = map(states) do state
        orbit_at(t_t) = calc_satellite_position_and_velocity(state.decoder, t_t).position
        # The transmit time that puts this satellite on the true position…
        _, transmit_time = modelled_pseudorange(orbit_at, t_t -> 0.0, t -> user_pos,
            t -> clock_bias, reception)
        # …and the exact Doppler: the negated rate of the pseudorange there, with the
        # receiver moving and drifting and the satellite clock running at its drift.
        sat_clock_drift = PositionVelocityTime.calc_satellite_clock_drift(state.decoder,
            transmit_time)
        pseudorange(t) = first(modelled_pseudorange(orbit_at,
            t_t -> c * sat_clock_drift * (t_t - transmit_time),
            t -> user_pos + user_vel * (t - reception),
            t -> clock_bias + clock_drift * (t - reception), t))
        λ = c / ustrip(Hz, get_center_frequency(state.system))
        doppler = -central_difference(pseudorange, reception, 0.1) / λ
        retimed_state(state, transmit_time, doppler)
    end

    pvt = calc_pvt(signal_group(exact_states); kwargs...)
    # The transmit times are seconds-of-week counts, resolved in Float64 to ~6e-11 s —
    # ~2 cm of range — so the position is exact only to millimetres.
    @test norm(pvt.position - user_pos) < 1e-2
    @test norm(pvt.velocity - user_vel) < 1e-5
    @test pvt.relative_clock_drift * c ≈ clock_drift atol = 1e-5
    @test maximum(info -> abs(ustrip(m / s, info.rate_residual)), pvt.sats) < 1e-5
end
