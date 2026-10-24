# The atmospheric delays of a whole epoch, ionospheric and tropospheric together, as the
# solve removes them from the pseudoranges.

"""
    predict_atmospheric_delays(user_position, measurements, model::IonosphericModel,
                               time_of_week, day_of_year, enable_tropospheric_correction)
        -> Vector{Float64}

Per-satellite ionospheric + tropospheric delay (metres), to be subtracted from the
pseudoranges, seen from `user_position` (ECEF, metres: an `ECEF`, an `SVector{3}`, or any
vector whose first three elements are the position, such as a least-squares state vector
`[x, y, z, tc₁, …]`). `model` carries the constellation-wide ionospheric parameters from
[`select_ionospheric_parameters`](@ref) (`IonosphericModel(nothing)` skips the
ionosphere); the troposphere uses the blind Saastamoinen zenith delays mapped by the
Niell mapping functions, whose seasonal term takes `day_of_year` (see
[`tropospheric_delay`](@ref)), unless `enable_tropospheric_correction` is `false`.

A single corrected solve is enough: the delays depend on position only through the
satellite elevation/azimuth (and, for the troposphere, the user height), and
∂delay/∂position is negligible over the metre-level position uncertainty (a 15 m
shift moves the elevation by ~1e-5°), so delays predicted at a nearby position are
accurate to well under a millimetre — no iterate-to-convergence needed. The user
geodetic coordinates and the ENU transform depend only on the user position, so they
are computed once and reused across satellites.

The allocating form of [`predict_atmospheric_delays!`](@ref).
"""
predict_atmospheric_delays(user_position, measurements, model, time_of_week,
    day_of_year, enable_tropospheric_correction) = predict_atmospheric_delays!(
    Vector{Float64}(undef, length(measurements)), user_position, measurements, model,
    time_of_week, day_of_year, enable_tropospheric_correction)

"""
    predict_atmospheric_delays!(delays, user_position, measurements,
                                model::IonosphericModel, time_of_week, day_of_year,
                                enable_tropospheric_correction)
        -> delays

[`predict_atmospheric_delays`](@ref) into `delays` (resized to the satellite count).
Allocates nothing.
"""
function predict_atmospheric_delays!(
    delays,
    user_position,
    measurements,
    model,
    time_of_week,
    day_of_year,
    enable_tropospheric_correction,
)
    # The per-satellite loop specialises once per ionospheric model, behind this barrier:
    # the field is a five-member `Union`, one more than Julia splits by itself, so the
    # split is spelled out.
    predict(parameters) = _predict_atmospheric_delays!(delays, user_position, measurements,
        parameters, time_of_week, day_of_year, enable_tropospheric_correction)
    parameters = model.parameters
    parameters isa KlobucharParameters && return predict(parameters)
    parameters isa BeiDouKlobucharParameters && return predict(parameters)
    parameters isa NTCMGParameters && return predict(parameters)
    parameters isa BDGIMParameters && return predict(parameters)
    return predict(nothing)
end

function _predict_atmospheric_delays!(
    delays,
    user_position,
    measurements,
    parameters,
    time_of_week,
    day_of_year,
    enable_tropospheric_correction,
)
    delays = resize!(delays, length(measurements))
    user_pos = ECEF(user_position[1], user_position[2], user_position[3])
    user_lla = LLAfromECEF(wgs84)(user_pos)
    enu_from_ecef = ENUfromECEF(user_pos, wgs84)
    for (j, measurement) in enumerate(measurements)
        elevation, azimuth = calc_elevation_and_azimuth(enu_from_ecef, measurement.position)
        iono = ionospheric_delay(
            parameters,
            measurement.center_frequency,
            elevation,
            azimuth,
            user_lla,
            time_of_week,
        )
        tropo =
            enable_tropospheric_correction ?
            tropospheric_delay(elevation, user_lla, day_of_year) : 0.0
        delays[j] = iono + tropo
    end
    delays
end
