function calc_eccentric_anomaly(mean_anomaly, eccentricity)
    Ek = mean_anomaly
    for k = 1:30
        Et = Ek
        Ek = mean_anomaly + eccentricity * sin(Ek)
        if abs(Ek - Et) <= 1e-12
            break
        end
    end
    return Ek
end

function calc_eccentric_anomaly(eph::Ephemeris, t)
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - eph.t_0e)
    n = eph.n_0 + eph.n_dot_half * time_from_ephemeris_reference_epoch
    mean_anomaly = eph.M_0 + n * time_from_ephemeris_reference_epoch
    calc_eccentric_anomaly(mean_anomaly, eph.e)
end

calc_eccentric_anomaly(decoder::GNSSDecoder.GNSSDecoderState, t) =
    calc_eccentric_anomaly(Ephemeris(decoder), t)

"""
    calc_satellite_position(decoder::GNSSDecoder.GNSSDecoderState, t)
    calc_satellite_position(state::SatelliteState)

Calculate the satellite ECEF position from orbital parameters at time `t`.

The first method takes a decoder state and explicit time. The second method
computes the corrected transmission time from the [`SatelliteState`](@ref) automatically.

# Arguments
- `decoder`: GNSS decoder state containing ephemeris data
- `t`: Transmission time in system time (seconds)
- `state`: A [`SatelliteState`](@ref) combining decoder, system, and phase measurements

# Returns
An `SVector{3, Float64}` with the satellite position in ECEF coordinates (meters).
"""
function calc_satellite_position(decoder::GNSSDecoder.GNSSDecoderState, t)
    pos_and_vel = calc_satellite_position_and_velocity(decoder, t)
    pos_and_vel.position
end

"""
    calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t)
    calc_satellite_position_and_velocity(eph::Ephemeris, t)
    calc_satellite_position_and_velocity(state::SatelliteState)

Calculate the satellite ECEF position and velocity from orbital parameters at time `t`.

Uses Keplerian orbital mechanics with perturbation corrections (harmonic corrections
for argument of latitude, radius, and inclination) to propagate the satellite ephemeris.

# Arguments
- `decoder`: GNSS decoder state containing ephemeris data
- `eph`: An already-extracted [`Ephemeris`](@ref) — the form `calc_pvt` uses to reuse
  one extraction per satellite across the clock and position kernels
- `t`: Transmission time in system time (seconds)
- `state`: A [`SatelliteState`](@ref) combining decoder, system, and phase measurements

# Returns
A named tuple `(position, velocity)` where each is an `SVector{3, Float64}` in ECEF
coordinates (meters and m/s respectively).
"""
function calc_satellite_position_and_velocity(eph::Ephemeris, t)
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - eph.t_0e)
    # Semi-major axis at t_k: constant for LNAV/Galileo (A_dot = 0), `A_0 + Ȧ·t_k`
    # for CNAV/CNAV-2; same for the mean motion and its rate `n_dot_half`.
    semi_major_axis = eph.A + eph.A_dot * time_from_ephemeris_reference_epoch
    corrected_mean_motion = eph.n_0 + eph.n_dot_half * time_from_ephemeris_reference_epoch
    mean_anomaly = eph.M_0 + corrected_mean_motion * time_from_ephemeris_reference_epoch
    eccentric_anomaly = calc_eccentric_anomaly(mean_anomaly, eph.e)
    eccentric_anomaly_dot = corrected_mean_motion / (1.0 - eph.e * cos(eccentric_anomaly))
    β = eph.e / (1 + sqrt(1 - eph.e^2))
    true_anomaly =
        eccentric_anomaly +
        2 * atan(β * sin(eccentric_anomaly) / (1 - β * cos(eccentric_anomaly)))
    true_anomaly_dot =
        sin(eccentric_anomaly) *
        eccentric_anomaly_dot *
        (1.0 + eph.e * cos(true_anomaly)) /
        (sin(true_anomaly) * (1.0 - eph.e * cos(eccentric_anomaly)))
    argument_of_latitude = true_anomaly + eph.ω
    argrument_of_latitude_correction =
        eph.C_us * sin(2 * argument_of_latitude) +
        eph.C_uc * cos(2 * argument_of_latitude)
    radius_correction =
        eph.C_rs * sin(2 * argument_of_latitude) +
        eph.C_rc * cos(2 * argument_of_latitude)
    inclination_correction =
        eph.C_is * sin(2 * argument_of_latitude) +
        eph.C_ic * cos(2 * argument_of_latitude)
    corrected_argument_of_latitude = argument_of_latitude + argrument_of_latitude_correction
    corrected_radius =
        semi_major_axis * (1 - eph.e * cos(eccentric_anomaly)) + radius_correction
    corrected_inclination =
        eph.i_0 + inclination_correction + eph.i_dot * time_from_ephemeris_reference_epoch

    corrected_argument_of_latitude_dot =
        true_anomaly_dot +
        2 *
        (
            eph.C_us * cos(2 * corrected_argument_of_latitude) -
            eph.C_uc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_radius_dot =
        eph.A_dot * (1.0 - eph.e * cos(eccentric_anomaly)) +
        semi_major_axis * eph.e * sin(eccentric_anomaly) * corrected_mean_motion /
        (1.0 - eph.e * cos(eccentric_anomaly)) +
        2 *
        (
            eph.C_rs * cos(2 * corrected_argument_of_latitude) -
            eph.C_rc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_inclination_dot =
        eph.i_dot +
        (
            eph.C_is * cos(2 * corrected_argument_of_latitude) -
            eph.C_ic * sin(2 * corrected_argument_of_latitude)
        ) *
        2 *
        true_anomaly_dot

    x_position_in_orbital_plane = corrected_radius * cos(corrected_argument_of_latitude)
    y_position_in_orbital_plane = corrected_radius * sin(corrected_argument_of_latitude)

    x_position_in_orbital_plane_dot =
        corrected_radius_dot * cos(corrected_argument_of_latitude) -
        y_position_in_orbital_plane * corrected_argument_of_latitude_dot
    y_position_in_orbital_plane_dot =
        corrected_radius_dot * sin(corrected_argument_of_latitude) +
        x_position_in_orbital_plane * corrected_argument_of_latitude_dot

    corrected_longitude_of_ascending_node =
        eph.Ω_0 + (eph.Ω_dot - eph.Ω_dot_e) * time_from_ephemeris_reference_epoch -
        eph.Ω_dot_e * eph.t_0e

    corrected_longitude_of_ascending_node_dot = eph.Ω_dot - eph.Ω_dot_e

    position = SVector(
        x_position_in_orbital_plane * cos(corrected_longitude_of_ascending_node) -
        y_position_in_orbital_plane *
        cos(corrected_inclination) *
        sin(corrected_longitude_of_ascending_node),
        x_position_in_orbital_plane * sin(corrected_longitude_of_ascending_node) +
        y_position_in_orbital_plane *
        cos(corrected_inclination) *
        cos(corrected_longitude_of_ascending_node),
        y_position_in_orbital_plane * sin(corrected_inclination),
    )

    velocity = SVector(
        (
            x_position_in_orbital_plane_dot -
            y_position_in_orbital_plane *
            cos(corrected_inclination) *
            corrected_longitude_of_ascending_node_dot
        ) * cos(corrected_longitude_of_ascending_node) -
        (
            x_position_in_orbital_plane * corrected_longitude_of_ascending_node_dot +
            y_position_in_orbital_plane_dot * cos(corrected_inclination) -
            y_position_in_orbital_plane *
            sin(corrected_inclination) *
            corrected_inclination_dot
        ) * sin(corrected_longitude_of_ascending_node),
        (
            x_position_in_orbital_plane_dot -
            y_position_in_orbital_plane *
            cos(corrected_inclination) *
            corrected_longitude_of_ascending_node_dot
        ) * sin(corrected_longitude_of_ascending_node) +
        (
            x_position_in_orbital_plane * corrected_longitude_of_ascending_node_dot +
            y_position_in_orbital_plane_dot * cos(corrected_inclination) -
            y_position_in_orbital_plane *
            sin(corrected_inclination) *
            corrected_inclination_dot
        ) * cos(corrected_longitude_of_ascending_node),
        y_position_in_orbital_plane_dot * sin(corrected_inclination) +
        y_position_in_orbital_plane *
        cos(corrected_inclination) *
        corrected_inclination_dot,
    )
    (position = position, velocity = velocity)
end

calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t) =
    calc_satellite_position_and_velocity(Ephemeris(decoder), t)

function calc_satellite_position(state::SatelliteState)
    pos_and_vel = calc_satellite_position_and_velocity(state)
    pos_and_vel.position
end
function calc_satellite_position_and_velocity(state::SatelliteState)
    clock = ClockModel(state.decoder, state.system)
    t = calc_corrected_time(state, clock)
    calc_satellite_position_and_velocity(clock.ephemeris, t)
end

function calc_pseudo_ranges(times)
    t_ref = maximum(times)
    reference_times = map(time -> t_ref - time, times)
    pseudoranges = reference_times .* SPEEDOFLIGHT
    return pseudoranges, t_ref
end