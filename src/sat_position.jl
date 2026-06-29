# GPS CNAV / CNAV-2 quasi-Keplerian reference values (IS-GPS-200N CNAV user
# algorithm; identical in IS-GPS-705J and IS-GPS-800J).
const GPS_SEMI_MAJOR_AXIS_REFERENCE = 26_559_710.0          # A_REF (m)
const GPS_RATE_OF_RIGHT_ASCENSION_REFERENCE = -2.6e-9 * π   # Ω̇_REF (rad/s; -2.6e-9 semicircles/s)

"""
    orbital_elements(data, μ, t_k) -> (; A, sqrt_A, A_dot, n, Ω_dot)

The Keplerian elements that differ between the directly-broadcast ephemerides
(GPS LNAV `GPSL1CAData`, Galileo I/NAV `GalileoE1BData`) and the quasi-Keplerian
GPS CNAV/CNAV-2 ephemerides (`GPSCNAVData` for L5 and L2C, `GPSL1C_DData` for
L1C), evaluated at the time-from-ephemeris `t_k` (seconds):

- `A`: semi-major axis at `t_0e` (m)
- `sqrt_A`: its square root (m^½); the broadcast `sqrt_A` for the directly-broadcast
  Keplerian case (no round-trip), `√A` for CNAV/CNAV-2 (which carry no `sqrt_A` field)
- `A_dot`: its rate (m/s; `0` for the directly-broadcast Keplerian case)
- `n`: corrected mean motion at `t_k` (rad/s)
- `Ω_dot`: rate of right ascension (rad/s)

Everything else the propagator needs (`e`, `ω`, `i_0`, `i_dot`, `M_0`, the `C_*`
harmonic coefficients, `t_0e`) is named identically across all four nav messages
and read from `data` directly. CNAV recovers `A` from `A_REF + ΔA` (with `Ȧ·t_k`
added for the radius), the mean motion from `Δn_0 (+ ½ Δṅ_0 t_k)`, and `Ω̇` from
`Ω̇_REF + ΔΩ̇`.
"""
function orbital_elements(data::GNSSDecoder.AbstractGNSSData, μ, t_k)
    (A = data.sqrt_A^2, sqrt_A = data.sqrt_A, A_dot = 0.0, n = sqrt(μ) / data.sqrt_A^3 + data.Δn, Ω_dot = data.Ω_dot)
end
function orbital_elements(data::GPSModernNavData, μ, t_k)
    A = GPS_SEMI_MAJOR_AXIS_REFERENCE + data.ΔA
    n = sqrt(μ / A^3) + data.Δn_0 + 0.5 * data.Δn_0_dot * t_k
    Ω_dot = GPS_RATE_OF_RIGHT_ASCENSION_REFERENCE + data.ΔΩ_dot
    (A = A, sqrt_A = sqrt(A), A_dot = data.A_dot, n = n, Ω_dot = Ω_dot)
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

function calc_eccentric_anomaly(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - data.t_0e)
    el = orbital_elements(data, decoder.constants.μ, time_from_ephemeris_reference_epoch)
    mean_anomaly = data.M_0 + el.n * time_from_ephemeris_reference_epoch
    calc_eccentric_anomaly(mean_anomaly, data.e)
end

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
    calc_satellite_position_and_velocity(state::SatelliteState)

Calculate the satellite ECEF position and velocity from orbital parameters at time `t`.

Uses Keplerian orbital mechanics with perturbation corrections (harmonic corrections
for argument of latitude, radius, and inclination) to propagate the satellite ephemeris.

# Arguments
- `decoder`: GNSS decoder state containing ephemeris data
- `t`: Transmission time in system time (seconds)
- `state`: A [`SatelliteState`](@ref) combining decoder, system, and phase measurements

# Returns
A named tuple `(position, velocity)` where each is an `SVector{3, Float64}` in ECEF
coordinates (meters and m/s respectively).
"""
function calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    constants = decoder.constants
    time_from_ephemeris_reference_epoch = correct_week_crossovers(t - data.t_0e)
    el = orbital_elements(data, constants.μ, time_from_ephemeris_reference_epoch)
    # Semi-major axis at t_k: constant for LNAV/Galileo (A_dot = 0), `A_0 + Ȧ·t_k`
    # for CNAV/CNAV-2.
    semi_major_axis = el.A + el.A_dot * time_from_ephemeris_reference_epoch
    corrected_mean_motion = el.n
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
        el.A_dot * (1.0 - data.e * cos(eccentric_anomaly)) +
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
        data.Ω_0 + (el.Ω_dot - constants.Ω_dot_e) * time_from_ephemeris_reference_epoch -
        constants.Ω_dot_e * data.t_0e

    corrected_longitude_of_ascending_node_dot = el.Ω_dot - constants.Ω_dot_e

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

function calc_satellite_position(state::SatelliteState)
    pos_and_vel = calc_satellite_position_and_velocity(state)
    pos_and_vel.position
end
function calc_satellite_position_and_velocity(state::SatelliteState)
    t = calc_corrected_time(state)
    calc_satellite_position_and_velocity(state.decoder, t)
end

function calc_pseudo_ranges(times)
    t_ref = maximum(times)
    reference_times = map(time -> t_ref - time, times)
    pseudoranges = reference_times .* SPEEDOFLIGHT
    return pseudoranges, t_ref
end