"""
    _rotate_beidou_geo(position, velocity, Ω_dot_e, t_k) -> (; position, velocity)

Rotate a BeiDou GEO satellite from the auxiliary frame the user algorithm builds
it in into BDCS (BDS-SIS-ICD-B1I-3.0 §5.2.4.12):

    [X; Y; Z] = R_z(ω̇_e · t_k) · R_x(−5°) · [X_GK; Y_GK; Z_GK]

The Keplerian projection is ill-conditioned for a GEO orbit (`i ≈ 0`, and `Ω_k` and the
Earth rotation cancel), so the ICD propagates it in a frame inertial over the fit
interval (`Ω_k` without `−ω̇_e·t_k`; see the `is_geostationary` branch), tilted 5° off the
equator, then turns it into BDCS. The velocity's `Ṙ_z(ω̇_e·t_k)` term, easily dropped, is
~3 km/s at GEO radius: the whole velocity.
"""
function _rotate_beidou_geo(position, velocity, Ω_dot_e, t_k)
    φ = Ω_dot_e * t_k
    sφ, cφ = sincos(φ)
    # Inclination of the auxiliary frame relative to BDCS: −5° (§5.2.4.12).
    st, ct = sincos(-5.0 * π / 180)
    R_z = @SMatrix [cφ sφ 0.0; -sφ cφ 0.0; 0.0 0.0 1.0]
    R_z_dot = Ω_dot_e * @SMatrix [-sφ cφ 0.0; -cφ -sφ 0.0; 0.0 0.0 0.0]
    R_x = @SMatrix [1.0 0.0 0.0; 0.0 ct st; 0.0 -st ct]
    tilted_position = R_x * position
    (
        position = R_z * tilted_position,
        velocity = R_z_dot * tilted_position + R_z * (R_x * velocity),
    )
end

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

# The corrected mean motion `t_k` seconds from `t_0e` (rad/s), the quasi-Keplerian form
# every message is normalised to (the Keplerian ones have zero rates).
mean_motion(ephemeris::BroadcastEphemeris, t_k) =
    sqrt(ephemeris.μ / ephemeris.A^3) + ephemeris.Δn_0 + 0.5 * ephemeris.Δn_0_dot * t_k

function calc_eccentric_anomaly(ephemeris::BroadcastEphemeris, t)
    t_k = fold_week_crossover(t - ephemeris.t_0e)
    calc_eccentric_anomaly(ephemeris.M_0 + mean_motion(ephemeris, t_k) * t_k, ephemeris.e)
end

"""
    calc_satellite_position_and_velocity(ephemeris::BroadcastEphemeris, t)
    calc_satellite_position_and_velocity(state::SatelliteState)

The satellite's ECEF `(position, velocity)` (`SVector{3}`s, m and m/s) from GNSSDecoder's
`BroadcastEphemeris` (`get_broadcast_ephemeris`), at the transmit time `t` (system time of
week, seconds) or at the one [`calc_transmit_time`](@ref) computes for `state`.
"""
function calc_satellite_position_and_velocity(ephemeris::BroadcastEphemeris, t)
    t_0e = ephemeris.t_0e
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - t_0e)
    # Semi-major axis at t_k (`A_dot = 0` where broadcast directly).
    semi_major_axis = ephemeris.A + ephemeris.A_dot * time_from_ephemeris_reference_epoch
    corrected_mean_motion = mean_motion(ephemeris, time_from_ephemeris_reference_epoch)
    eccentric_anomaly = calc_eccentric_anomaly(ephemeris, t)
    eccentric_anomaly_dot = corrected_mean_motion / (1.0 - ephemeris.e * cos(eccentric_anomaly))
    β = ephemeris.e / (1 + sqrt(1 - ephemeris.e^2))
    true_anomaly =
        eccentric_anomaly +
        2 * atan(β * sin(eccentric_anomaly) / (1 - β * cos(eccentric_anomaly)))
    # ν̇ without the textbook `sin E / sin ν`, which is 0/0 at perigee and apogee:
    # `sin ν = √(1−e²)·sin E / (1 − e·cos E)` turns it into `Ė·(1 + e·cos ν)/√(1−e²)`.
    true_anomaly_dot =
        eccentric_anomaly_dot * (1.0 + ephemeris.e * cos(true_anomaly)) / sqrt(1.0 - ephemeris.e^2)
    argument_of_latitude = true_anomaly + ephemeris.ω
    argument_of_latitude_correction =
        ephemeris.C_us * sin(2 * argument_of_latitude) +
        ephemeris.C_uc * cos(2 * argument_of_latitude)
    radius_correction =
        ephemeris.C_rs * sin(2 * argument_of_latitude) +
        ephemeris.C_rc * cos(2 * argument_of_latitude)
    inclination_correction =
        ephemeris.C_is * sin(2 * argument_of_latitude) +
        ephemeris.C_ic * cos(2 * argument_of_latitude)
    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
    corrected_radius =
        semi_major_axis * (1 - ephemeris.e * cos(eccentric_anomaly)) + radius_correction
    corrected_inclination =
        ephemeris.i_0 + inclination_correction + ephemeris.i_dot * time_from_ephemeris_reference_epoch

    corrected_argument_of_latitude_dot =
        true_anomaly_dot +
        2 *
        (
            ephemeris.C_us * cos(2 * corrected_argument_of_latitude) -
            ephemeris.C_uc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_radius_dot =
        ephemeris.A_dot * (1.0 - ephemeris.e * cos(eccentric_anomaly)) +
        semi_major_axis * ephemeris.e * sin(eccentric_anomaly) * corrected_mean_motion /
        (1.0 - ephemeris.e * cos(eccentric_anomaly)) +
        2 *
        (
            ephemeris.C_rs * cos(2 * corrected_argument_of_latitude) -
            ephemeris.C_rc * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_inclination_dot =
        ephemeris.i_dot +
        (
            ephemeris.C_is * cos(2 * corrected_argument_of_latitude) -
            ephemeris.C_ic * sin(2 * corrected_argument_of_latitude)
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

    # A BeiDou GEO satellite's Ω_k leaves out `−ω_e·t_k`: it is projected into an
    # inertial frame and turned into ECEF by `_rotate_beidou_geo` (see there).
    geo = ephemeris.is_geostationary
    corrected_longitude_of_ascending_node =
        ephemeris.Ω_0 + (ephemeris.Ω_dot - (geo ? 0.0 : ephemeris.Ω_dot_e)) *
                   time_from_ephemeris_reference_epoch - ephemeris.Ω_dot_e * t_0e

    corrected_longitude_of_ascending_node_dot =
        ephemeris.Ω_dot - (geo ? 0.0 : ephemeris.Ω_dot_e)

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
    geo ? _rotate_beidou_geo(position, velocity, ephemeris.Ω_dot_e,
        time_from_ephemeris_reference_epoch) : (position = position, velocity = velocity)
end

calc_satellite_position_and_velocity(state::SatelliteState) =
    calc_satellite_position_and_velocity(something(get_broadcast_ephemeris(state.decoder)),
        calc_transmit_time(state))
