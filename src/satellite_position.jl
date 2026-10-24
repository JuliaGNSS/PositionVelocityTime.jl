"""
    is_geo_orbit(decoder) -> Bool

Whether the satellite is geostationary, which selects the BeiDou GEO branch of
[`calc_satellite_position_and_velocity`](@ref). GNSSDecoder answers the orbit class. The
test must be for GEO: on B1I/B3I the PRN partition separates GEO from non-GEO only, and
`get_orbit_class` returns `nothing` rather than guess between IGSO and MEO.
"""
is_geo_orbit(decoder::GNSSDecoder.GNSSDecoderState) =
    get_orbit_class(decoder) === GNSSDecoder.geostationary_orbit

# A navigation-data field (`Union{Nothing,T}` until decoded) narrowed to its value, which
# complete data always has. Unnarrowed, Julia 1.10 boxes the intermediates of an
# expression over several such fields: an allocation per field read.
decoded(x) = something(x)

"""
    orbital_elements(data, μ, t_k) -> (; A, sqrt_A, A_dot, n, Ω_dot)

The Keplerian elements that differ between the directly broadcast ephemerides (GPS LNAV,
Galileo I/NAV and F/NAV, BeiDou D1/D2) and the quasi-Keplerian ones (GPS CNAV and CNAV-2,
BeiDou B-CNAV), at the time from ephemeris `t_k` (seconds):

- `A`: semi-major axis at `t_0e` (m)
- `sqrt_A`: its square root (m^½), the broadcast one where there is one
- `A_dot`: its rate (m/s; `0` when broadcast directly)
- `n`: corrected mean motion at `t_k` (rad/s)
- `Ω_dot`: rate of right ascension (rad/s)

The rest of the propagator's inputs are named alike in every message and read from `data`.
"""
function orbital_elements(data::GNSSDecoder.AbstractGNSSData, μ, t_k)
    (A = decoded(data.sqrt_A)^2, sqrt_A = decoded(data.sqrt_A), A_dot = 0.0, n = sqrt(μ) / decoded(data.sqrt_A)^3 + decoded(data.Δn), Ω_dot = decoded(data.Ω_dot))
end
function orbital_elements(data::AbstractGPSCNAVData, μ, t_k)
    # The CNAV reference values the broadcast deltas are off (IS-GPS-200N, as
    # IS-GPS-705J and IS-GPS-800J).
    A_REF = 26_559_710.0        # m
    Ω_dot_REF = -2.6e-9 * π     # rad/s (-2.6e-9 semicircles/s)
    A = A_REF + decoded(data.ΔA)
    n = sqrt(μ / A^3) + decoded(data.Δn_0) + 0.5 * decoded(data.Δn_0_dot) * t_k
    Ω_dot = Ω_dot_REF + decoded(data.ΔΩ_dot)
    (A = A, sqrt_A = sqrt(A), A_dot = decoded(data.A_dot), n = n, Ω_dot = Ω_dot)
end

"""
Reference semi-major axis for a BeiDou B-CNAV `sat_type` (BDS-SIS-ICD-B1C-1.0 Table 7-3
and the B2a/B2b counterparts): 27 906 100 m for MEO (3), 42 162 200 m for IGSO and GEO (2
and 1). The reserved 0 throws, as it leaves the orbit undefined; GNSSDecoder's
`is_decoding_completed_for_positioning` keeps it from reaching [`calc_pvt`](@ref).
"""
function beidou_reference_semi_major_axis(sat_type)
    sat_type == 3 && return 27_906_100.0            # MEO
    sat_type in (1, 2) && return 42_162_200.0       # GEO, IGSO
    throw(
        ArgumentError(
            "BeiDou sat_type $sat_type is reserved; the semi-major-axis reference " *
            "the broadcast ΔA is relative to is unknown. GNSSDecoder gates " *
            "is_decoding_completed_for_positioning on this, so a satellite " *
            "reaching the propagator cannot carry it.",
        ),
    )
end

# The BDS-3 B-CNAV user algorithm. Unlike GPS CNAV, `A_REF` depends on `sat_type`, and
# `Ω̇` is broadcast outright (BDS-SIS-ICD-B2a-1.0 Table 7-8: 19 bits at 2⁻⁴⁴ π/s), not as a
# delta: adding `Ω̇_REF` would cost ~50 m along-track after an hour.
function orbital_elements(data::AbstractBeiDouCNAVData, μ, t_k)
    A = beidou_reference_semi_major_axis(decoded(data.sat_type)) + decoded(data.ΔA)
    n = sqrt(μ / A^3) + decoded(data.Δn_0) + 0.5 * decoded(data.Δn_0_dot) * t_k
    (A = A, sqrt_A = sqrt(A), A_dot = decoded(data.A_dot), n = n, Ω_dot = decoded(data.Ω_dot))
end

"""
    rotate_beidou_geo(position, velocity, Ω_dot_e, t_k) -> (; position, velocity)

Rotate a BeiDou GEO satellite from the auxiliary frame the user algorithm builds
it in into BDCS (BDS-SIS-ICD-B1I-3.0 §5.2.4.12):

    [X; Y; Z] = R_z(ω̇_e · t_k) · R_x(−5°) · [X_GK; Y_GK; Z_GK]

The Keplerian projection is ill-conditioned for a GEO orbit (`i ≈ 0`, and `Ω_k` and the
Earth rotation cancel), so the ICD propagates it in a frame inertial over the fit
interval (`Ω_k` without `−ω̇_e·t_k`; see the `is_geo_orbit` branch), tilted 5° off the
equator, then turns it into BDCS. The velocity's `Ṙ_z(ω̇_e·t_k)` term, easily dropped, is
~3 km/s at GEO radius: the whole velocity.
"""
function rotate_beidou_geo(position, velocity, Ω_dot_e, t_k)
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

function calc_eccentric_anomaly(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - decoded(data.t_0e))
    el = orbital_elements(data, decoder.constants.μ, time_from_ephemeris_reference_epoch)
    mean_anomaly = decoded(data.M_0) + el.n * time_from_ephemeris_reference_epoch
    calc_eccentric_anomaly(mean_anomaly, decoded(data.e))
end

"""
    calc_satellite_position(decoder::GNSSDecoder.GNSSDecoderState, t)
    calc_satellite_position(state::SatelliteState)

The satellite's ECEF position (an `SVector{3}`, metres) at the transmit time `t` (system
time of week, seconds), or at the one [`calc_transmit_time`](@ref) computes for `state`.
"""
function calc_satellite_position(decoder::GNSSDecoder.GNSSDecoderState, t)
    pos_and_vel = calc_satellite_position_and_velocity(decoder, t)
    pos_and_vel.position
end

"""
    calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t)
    calc_satellite_position_and_velocity(state::SatelliteState)

The satellite's ECEF `(position, velocity)` (`SVector{3}`s, m and m/s) from the broadcast
ephemeris, at the transmit time `t` (system time of week, seconds) or at the one
[`calc_transmit_time`](@ref) computes for `state`.
"""
function calc_satellite_position_and_velocity(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    constants = decoder.constants
    t_0e = decoded(data.t_0e)
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - t_0e)
    el = orbital_elements(data, constants.μ, time_from_ephemeris_reference_epoch)
    # Semi-major axis at t_k (`A_dot = 0` where broadcast directly).
    semi_major_axis = el.A + el.A_dot * time_from_ephemeris_reference_epoch
    corrected_mean_motion = el.n
    eccentric_anomaly = calc_eccentric_anomaly(decoder, t)
    eccentric_anomaly_dot = corrected_mean_motion / (1.0 - decoded(data.e) * cos(eccentric_anomaly))
    β = decoded(data.e) / (1 + sqrt(1 - decoded(data.e)^2))
    true_anomaly =
        eccentric_anomaly +
        2 * atan(β * sin(eccentric_anomaly) / (1 - β * cos(eccentric_anomaly)))
    # ν̇ without the textbook `sin E / sin ν`, which is 0/0 at perigee and apogee:
    # `sin ν = √(1−e²)·sin E / (1 − e·cos E)` turns it into `Ė·(1 + e·cos ν)/√(1−e²)`.
    true_anomaly_dot =
        eccentric_anomaly_dot * (1.0 + decoded(data.e) * cos(true_anomaly)) / sqrt(1.0 - decoded(data.e)^2)
    argument_of_latitude = true_anomaly + decoded(data.ω)
    argument_of_latitude_correction =
        decoded(data.C_us) * sin(2 * argument_of_latitude) +
        decoded(data.C_uc) * cos(2 * argument_of_latitude)
    radius_correction =
        decoded(data.C_rs) * sin(2 * argument_of_latitude) +
        decoded(data.C_rc) * cos(2 * argument_of_latitude)
    inclination_correction =
        decoded(data.C_is) * sin(2 * argument_of_latitude) +
        decoded(data.C_ic) * cos(2 * argument_of_latitude)
    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
    corrected_radius =
        semi_major_axis * (1 - decoded(data.e) * cos(eccentric_anomaly)) + radius_correction
    corrected_inclination =
        decoded(data.i_0) + inclination_correction + decoded(data.i_dot) * time_from_ephemeris_reference_epoch

    corrected_argument_of_latitude_dot =
        true_anomaly_dot +
        2 *
        (
            decoded(data.C_us) * cos(2 * corrected_argument_of_latitude) -
            decoded(data.C_uc) * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_radius_dot =
        el.A_dot * (1.0 - decoded(data.e) * cos(eccentric_anomaly)) +
        semi_major_axis * decoded(data.e) * sin(eccentric_anomaly) * corrected_mean_motion /
        (1.0 - decoded(data.e) * cos(eccentric_anomaly)) +
        2 *
        (
            decoded(data.C_rs) * cos(2 * corrected_argument_of_latitude) -
            decoded(data.C_rc) * sin(2 * corrected_argument_of_latitude)
        ) *
        true_anomaly_dot
    corrected_inclination_dot =
        decoded(data.i_dot) +
        (
            decoded(data.C_is) * cos(2 * corrected_argument_of_latitude) -
            decoded(data.C_ic) * sin(2 * corrected_argument_of_latitude)
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
    # inertial frame and turned into ECEF by `rotate_beidou_geo` (see there).
    geo = is_geo_orbit(decoder)
    corrected_longitude_of_ascending_node =
        decoded(data.Ω_0) + (el.Ω_dot - (geo ? 0.0 : constants.Ω_dot_e)) *
                   time_from_ephemeris_reference_epoch - constants.Ω_dot_e * t_0e

    corrected_longitude_of_ascending_node_dot =
        el.Ω_dot - (geo ? 0.0 : constants.Ω_dot_e)

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
    geo ? rotate_beidou_geo(position, velocity, constants.Ω_dot_e,
        time_from_ephemeris_reference_epoch) : (position = position, velocity = velocity)
end

function calc_satellite_position(state::SatelliteState)
    pos_and_vel = calc_satellite_position_and_velocity(state)
    pos_and_vel.position
end
function calc_satellite_position_and_velocity(state::SatelliteState)
    t = calc_transmit_time(state)
    calc_satellite_position_and_velocity(state.decoder, t)
end
