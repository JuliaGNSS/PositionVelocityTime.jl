# The atmospheric delays of a whole epoch, ionospheric and tropospheric together, as the
# solve removes them from the pseudoranges.

"""
    predict_atmospheric_delays(user_pos, measurements, model::IonosphericModel,
                               epoch::TAITime, enable_tropospheric_correction)
        -> Vector{Float64}

Per-satellite ionospheric plus tropospheric delay (m), to subtract from the pseudoranges,
seen from `user_pos` (ECEF, m; any vector whose first three elements are the
position, a state vector `[x, y, z, tc₁, …]` included). `model` holds the ionospheric
parameters ([`select_ionospheric_parameters`](@ref); `IonosphericModel(nothing)` skips the
ionosphere), and `epoch` is the TAI epoch the delays are predicted for, of which
[`tropospheric_delay`](@ref) takes the day of year, off its TAI calendar label.

Delays predicted at a nearby position are good enough: a 15 m shift moves the elevation by
~1e-5°, well under a millimetre of delay, so one corrected solve needs no iteration.
"""
predict_atmospheric_delays(user_pos, measurements, model, epoch,
    enable_tropospheric_correction) = predict_atmospheric_delays!(
    Vector{Float64}(undef, length(measurements)), user_pos, measurements, model,
    epoch, enable_tropospheric_correction)

"""
    predict_atmospheric_delays!(delays, user_pos, measurements,
                                model::IonosphericModel, epoch::TAITime,
                                enable_tropospheric_correction)
        -> delays

[`predict_atmospheric_delays`](@ref) into `delays`, resized to the satellite count.
Allocates nothing.
"""
function predict_atmospheric_delays!(
    delays,
    user_pos,
    measurements,
    model,
    epoch::TAITime,
    enable_tropospheric_correction,
)
    # A barrier that specialises the loop per model: Julia splits a `Union` of at most four
    # members by itself, and this one has five.
    predict(parameters) = _predict_atmospheric_delays!(delays, user_pos, measurements,
        parameters, epoch, enable_tropospheric_correction)
    parameters = model.parameters
    parameters isa KlobucharParameters && return predict(parameters)
    parameters isa BeiDouKlobucharParameters && return predict(parameters)
    parameters isa NeQuickGParameters && return predict(parameters)
    parameters isa BDGIMParameters && return predict(parameters)
    return predict(nothing)
end

function _predict_atmospheric_delays!(
    delays,
    user_pos,
    measurements,
    parameters,
    epoch,
    enable_tropospheric_correction,
)
    delays = resize!(delays, length(measurements))
    user_ecef = ECEF(user_pos[1], user_pos[2], user_pos[3])
    user_lla = LLAfromECEF(wgs84)(user_ecef)
    enu_from_ecef = ENUfromECEF(user_ecef, wgs84)
    doy = dayofyear(DateTime(epoch))
    for (j, measurement) in enumerate(measurements)
        elevation, azimuth = calc_elevation_and_azimuth(enu_from_ecef, measurement.position)
        iono = ionospheric_delay(
            parameters,
            measurement.center_frequency,
            elevation,
            azimuth,
            user_lla,
            epoch,
        )
        tropo =
            enable_tropospheric_correction ?
            tropospheric_delay(elevation, user_lla, doy) : 0.0
        delays[j] = iono + tropo
    end
    delays
end
