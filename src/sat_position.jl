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

function calc_eccentric_anomaly(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    computed_mean_motion = sqrt(decoder.constants.μ) / data.sqrt_A^3
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - data.t_0e)
    corrected_mean_motion = computed_mean_motion + data.Δn
    mean_anomaly = data.M_0 + corrected_mean_motion * time_from_ephemeris_reference_epoch
    calc_eccentric_anomaly(mean_anomaly, data.e)
end

"""
Calculates satellite position at time instance t

$SIGNATURES
`decoder`: GNSS decoder state
`t`: time at satellite in system time
"""
function calc_satellite_position(decoder::GNSSDecoder.GNSSDecoderState, t)
    pos_and_vel = calc_satellite_position_and_velocity(decoder, t)
    pos_and_vel.position
end

function calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    constants = decoder.constants
    semi_major_axis = data.sqrt_A^2
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - data.t_0e)
    computed_mean_motion = sqrt(decoder.constants.μ) / data.sqrt_A^3
    corrected_mean_motion = computed_mean_motion + data.Δn
    eccentric_anomaly = calc_eccentric_anomaly(decoder, t)
    eccentric_anomaly_dot = corrected_mean_motion / (1.0 - data.e * cos(eccentric_anomaly))
    β = data.e / (1 + sqrt(1 - data.e^2))
    true_anomaly =
        eccentric_anomaly +
        2 * atan(β * sin(eccentric_anomaly) / (1 - β * cos(eccentric_anomaly)))
    true_anomaly_dot =
        sin(eccentric_anomaly) *
        eccentric_anomaly_dot *
        (1.0 + data.e * cos(true_anomaly)) /
        (sin(true_anomaly) * (1.0 - data.e * cos(eccentric_anomaly)))
    argument_of_latitude = true_anomaly + data.ω
    argrument_of_latitude_correction =
        data.C_us * sin(2 * argument_of_latitude) +
        data.C_uc * cos(2 * argument_of_latitude)
    radius_correction =
        data.C_rs * sin(2 * argument_of_latitude) +
        data.C_rc * cos(2 * argument_of_latitude)
    inclination_correction =
        data.C_is * sin(2 * argument_of_latitude) +
        data.C_ic * cos(2 * argument_of_latitude)
    corrected_argument_of_latitude = argument_of_latitude + argrument_of_latitude_correction
    corrected_radius =
        semi_major_axis * (1 - data.e * cos(eccentric_anomaly)) + radius_correction
    corrected_inclination =
        data.i_0 + inclination_correction + data.i_dot * time_from_ephemeris_reference_epoch

    corrected_argument_of_latitude_dot =
        true_anomaly_dot +
        2 *
        (
            data.C_us * cos(2 * corrected_argument_of_latitude) -
            data.C_uc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_radius_dot =
        semi_major_axis * data.e * sin(eccentric_anomaly) * corrected_mean_motion /
        (1.0 - data.e * cos(eccentric_anomaly)) +
        2 *
        (
            data.C_rs * cos(2 * corrected_argument_of_latitude) -
            data.C_rc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_inclination_dot =
        data.i_dot +
        (
            data.C_is * cos(2 * corrected_argument_of_latitude) -
            data.C_ic * sin(2 * corrected_argument_of_latitude)
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
        data.Ω_0 + (data.Ω_dot - constants.Ω_dot_e) * time_from_ephemeris_reference_epoch -
        constants.Ω_dot_e * data.t_0e

    corrected_longitude_of_ascending_node_dot = data.Ω_dot - constants.Ω_dot_e

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

"""
Calculates satellite position at transmission time based on code phase and number of bits since TOW.

$SIGNATURES
`system`: GNSS system
`state`: Satellite state (SatelliteState)
"""
function calc_satellite_position(state::SatelliteState)
    pos_and_vel = calc_satellite_position_and_velocity(state)
    pos_and_vel.position
end
function calc_satellite_position_and_velocity(state::SatelliteState)
    t = calc_corrected_time(state)
    calc_satellite_position_and_velocity(state.decoder, t)
end

"""
Computes pseudo ranges 

$SIGNATURES
`sat_state`: satellite state, combining decoded data, code- and carrierphase 

Computes relative pseudo ranges of satellite vehicles.
The algorithm is based on the common reception method. 
"""
function calc_pseudo_ranges(times)
    t_ref = maximum(times)
    reference_times = map(time -> t_ref - time, times)
    pseudoranges = reference_times .* SPEEDOFLIGHT
    return pseudoranges, t_ref
end