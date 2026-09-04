"""
    is_geo_orbit(decoder) -> Bool

Whether the satellite is in a geostationary orbit, which selects the BeiDou GEO
branch of [`calc_satellite_position_and_velocity`](@ref) below.

One method for every signal, because GNSSDecoder answers the orbit class itself:
from the constellation supertype for GPS and Galileo, from the broadcast
`sat_type` on BeiDou B1C/B2a/B2b, and from the PRN partition on B1I/B3I.

`=== geostationary_orbit` is the whole test, and it has to be written that way
round. B1I/B3I can only answer the question partly — D1/D2 NAV broadcasts no
orbit-type field and the ICD's PRN partition separates GEO from non-GEO and no
further — so `get_orbit_class` returns `nothing` there rather than guess between
IGSO and MEO. Testing for GEO is exact on every signal; testing for MEO is not.
"""
is_geo_orbit(decoder::GNSSDecoder.GNSSDecoderState) =
    get_orbit_class(decoder) === GNSSDecoder.geostationary_orbit

"""
    orbital_elements(data, μ, t_k) -> (; A, sqrt_A, A_dot, n, Ω_dot)

The Keplerian elements that differ between the directly-broadcast ephemerides
(GPS LNAV `GPSL1CAData`, Galileo I/NAV `GalileoINAVData`, BeiDou D1/D2
`BeiDouDNAVData`) and the quasi-Keplerian ones (`GPSCNAVData` for L5 and L2C,
`GPSL1C_DData` for L1C, and the BeiDou B-CNAV family on B1C/B2a/B2b), evaluated
at the time-from-ephemeris `t_k` (seconds):

- `A`: semi-major axis at `t_0e` (m)
- `sqrt_A`: its square root (m^½); the broadcast `sqrt_A` for the directly-broadcast
  Keplerian case (no round-trip), `√A` for the quasi-Keplerian ones (which carry no
  `sqrt_A` field)
- `A_dot`: its rate (m/s; `0` for the directly-broadcast Keplerian case)
- `n`: corrected mean motion at `t_k` (rad/s)
- `Ω_dot`: rate of right ascension (rad/s)

Everything else the propagator needs (`t_0e`, `e`, `ω`, `i_0`, `i_dot`, `M_0`, the
`C_*` harmonic coefficients) is named identically across every nav message and read
from `data` directly. The quasi-Keplerian messages recover `A` from
`A_REF + ΔA` (with `Ȧ·t_k` added for the radius) and the mean motion from
`Δn_0 (+ ½ Δṅ_0 t_k)`; they differ in `Ω̇`, which GPS broadcasts as a delta off
`Ω̇_REF` and BeiDou broadcasts outright (see the `AbstractBeiDouCNAVData` method below).
"""
function orbital_elements(data::GNSSDecoder.AbstractGNSSData, μ, t_k)
    (A = data.sqrt_A^2, sqrt_A = data.sqrt_A, A_dot = 0.0, n = sqrt(μ) / data.sqrt_A^3 + data.Δn, Ω_dot = data.Ω_dot)
end
function orbital_elements(data::AbstractGPSCNAVData, μ, t_k)
    # Quasi-Keplerian reference values from the CNAV user algorithm (IS-GPS-200N;
    # identical in IS-GPS-705J and IS-GPS-800J); the broadcast fields are deltas off
    # these.
    A_REF = 26_559_710.0        # m
    Ω_dot_REF = -2.6e-9 * π     # rad/s (-2.6e-9 semicircles/s)
    A = A_REF + data.ΔA
    n = sqrt(μ / A^3) + data.Δn_0 + 0.5 * data.Δn_0_dot * t_k
    Ω_dot = Ω_dot_REF + data.ΔΩ_dot
    (A = A, sqrt_A = sqrt(A), A_dot = data.A_dot, n = n, Ω_dot = Ω_dot)
end

"""
Reference semi-major axis for a BeiDou B-CNAV `sat_type` (BDS-SIS-ICD-B1C-1.0
Table 7-3 and its B2a/B2b counterparts): 27 906 100 m for MEO (`sat_type` 3),
42 162 200 m for IGSO and GEO (2 and 1). The reserved `sat_type` 0 throws
rather than picking one: the two references differ by 14 256 100 m, so a
satellite that has not said which one its `ΔA` is relative to has not given its
orbit at all. It cannot arrive through [`calc_pvt`](@ref) — all three B-CNAV
messages gate `is_decoding_completed_for_positioning` on the orbit type being
one of the three defined ones.
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

# The BDS-3 B-CNAV user algorithm, the BeiDou counterpart of the GPS CNAV branch
# above it. Two differences from GPS, both easy to miss:
#
#   - `A_REF` is not one constant but two, selected by the broadcast `sat_type`
#     (MEO orbits sit ~14 000 km below the IGSO/GEO ones).
#   - `Ω̇` is broadcast *outright*, not as a delta off a reference rate the way
#     GPS CNAV broadcasts `ΔΩ̇` (BDS-SIS-ICD-B2a-1.0 Table 7-8: 19 bits at
#     2⁻⁴⁴ π/s, wide enough for the absolute value; GPS's 17-bit `ΔΩ̇` is not).
#     Adding a `Ω̇_REF` here would be a ~2.6e-9 semicircle/s error — about 50 m of
#     along-track position after an hour of propagation.
function orbital_elements(data::AbstractBeiDouCNAVData, μ, t_k)
    A = beidou_reference_semi_major_axis(data.sat_type) + data.ΔA
    n = sqrt(μ / A^3) + data.Δn_0 + 0.5 * data.Δn_0_dot * t_k
    (A = A, sqrt_A = sqrt(A), A_dot = data.A_dot, n = n, Ω_dot = data.Ω_dot)
end

"""
    _rotate_beidou_geo(position, velocity, Ω_dot_e, t_k) -> (; position, velocity)

Rotate a BeiDou GEO satellite from the auxiliary frame the user algorithm builds
it in into BDCS (BDS-SIS-ICD-B1I-3.0 §5.2.4.12):

    [X; Y; Z] = R_z(ω̇_e · t_k) · R_x(−5°) · [X_GK; Y_GK; Z_GK]

The reason this exists at all is that the standard Keplerian projection is
ill-conditioned for a GEO orbit: with `i ≈ 0` the ascending node is barely
defined, and a satellite that is nearly stationary in the earth-fixed frame has
`Ω_k` and the earth rotation cancelling to within the noise of the broadcast
elements. The ICD therefore propagates a GEO satellite in a frame that is
*inertial* over the fit interval — `Ω_k` carries no `−ω̇_e·t_k` term, which is
what the `is_geo_orbit` branch in `calc_satellite_position_and_velocity`
arranges — tilts it 5° off the equator so the orbit is not degenerate there
either, and only then spins it into BDCS by the earth rotation accumulated since
`t_0e`.

The velocity follows by the product rule; `Ṙ_z(ω̇_e·t_k) · ω̇_e` is the term a
naive implementation drops, and it is worth ~3 km/s at GEO radius — the whole
velocity, not a correction to it.
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

function calc_eccentric_anomaly(decoder::GNSSDecoder.GNSSDecoderState, t)
    data = decoder.data
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - data.t_0e)
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
    t_0e = data.t_0e
    time_from_ephemeris_reference_epoch = fold_week_crossover(t - t_0e)
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
    # Singularity-free form of ν̇. The textbook expression carries a `sin E / sin ν`
    # ratio, which is `0/0 → NaN` at perigee (`E = ν = 0`) and apogee (`E = ν = π`) —
    # points every satellite sweeps through twice per orbit. Substituting the identity
    # `sin ν = √(1−e²)·sin E / (1 − e·cos E)` cancels the ratio into the finite
    # `(1 − e·cos E)/√(1−e²)`, leaving `ν̇ = Ė·(1 + e·cos ν)/√(1−e²)`: identical away
    # from those points, finite everywhere.
    true_anomaly_dot =
        eccentric_anomaly_dot * (1.0 + data.e * cos(true_anomaly)) / sqrt(1.0 - data.e^2)
    argument_of_latitude = true_anomaly + data.ω
    argument_of_latitude_correction =
        data.C_us * sin(2 * argument_of_latitude) +
        data.C_uc * cos(2 * argument_of_latitude)
    radius_correction =
        data.C_rs * sin(2 * argument_of_latitude) +
        data.C_rc * cos(2 * argument_of_latitude)
    inclination_correction =
        data.C_is * sin(2 * argument_of_latitude) +
        data.C_ic * cos(2 * argument_of_latitude)
    corrected_argument_of_latitude = argument_of_latitude + argument_of_latitude_correction
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

    # The longitude of the ascending node, and the frame the orbital plane is then
    # projected into, are where the BeiDou GEO satellites part company with everything
    # else. A non-GEO satellite is projected straight into the earth-fixed frame, so the
    # earth-rotation term `−ω_e·t_k` is folded into Ω_k here. A BeiDou GEO satellite is
    # instead projected into a custom inertial frame — Ω_k carries no `−ω_e·t_k` — and
    # rotated into ECEF afterwards by `_rotate_beidou_geo`; see that function for why.
    geo = is_geo_orbit(decoder)
    corrected_longitude_of_ascending_node =
        data.Ω_0 + (el.Ω_dot - (geo ? 0.0 : constants.Ω_dot_e)) *
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
    geo ? _rotate_beidou_geo(position, velocity, constants.Ω_dot_e,
        time_from_ephemeris_reference_epoch) : (position = position, velocity = velocity)
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
    # Folded modulo the week, like the vector loop's `pseudorange_from_tows` and
    # every other time difference in this package: the inputs are seconds-of-week
    # counts, and the week wrap does not pass through all of them at once. Within
    # one constellation the wrap sweeps through the transmit-time spread (tens of
    # milliseconds, once a week); in a mixed solve the scale alignment
    # (`calc_time_scale_offsets`) widens the straddle to the full defined offset —
    # a BeiDou time reads `SOW + 14` on the GPS count, which exceeds 604800 for
    # the first 14 s of every GPS week while the GPS times have already wrapped.
    # Unfolded, every difference across the wrap is off by a week — 1.8e14 m.
    reference_times = map(time -> fold_week_crossover(t_ref - time), times)
    pseudoranges = reference_times .* SPEED_OF_LIGHT
    return pseudoranges, t_ref
end