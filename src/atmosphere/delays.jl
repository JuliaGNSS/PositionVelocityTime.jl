# The atmospheric delays of a whole epoch, ionospheric and tropospheric together, as the
# solve removes them from the pseudoranges.

"""
    predict_atmospheric_delays(user_position, measurements, model::IonosphericModel,
                               time_of_week, day_of_year, enable_tropospheric_correction)
        -> Vector{Float64}

Per-satellite ionospheric plus tropospheric delay (m), to subtract from the pseudoranges,
seen from `user_position` (ECEF, m; any vector whose first three elements are the
position, a state vector `[x, y, z, tc₁, …]` included). `model` holds the ionospheric
parameters ([`select_ionospheric_parameters`](@ref); `IonosphericModel(nothing)` skips the
ionosphere), and `day_of_year` feeds [`tropospheric_delay`](@ref).

Delays predicted at a nearby position are good enough: a 15 m shift moves the elevation by
~1e-5°, well under a millimetre of delay, so one corrected solve needs no iteration.
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

[`predict_atmospheric_delays`](@ref) into `delays`, resized to the satellite count.
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
    # A barrier that specialises the loop per model: Julia splits a `Union` of at most four
    # members by itself, and this one has five.
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
