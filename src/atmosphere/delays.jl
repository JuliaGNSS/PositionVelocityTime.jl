# The atmospheric delays of a whole epoch, ionospheric and tropospheric together, as the
# solve removes them from the pseudoranges.

"""
    predict_atmospheric_delays(user_pos, measurements, parameters,
                               reference_time, doy, enable_tropospheric_correction) -> Vector{Float64}

Per-satellite ionospheric + tropospheric delay (metres), to be subtracted from the
pseudoranges, seen from `user_pos` (ECEF, metres; any vector whose first three
elements are the position, a state vector `[x, y, z, tc₁, …]` included). `parameters` is
the constellation-wide set of ionospheric parameters
from [`select_ionospheric_parameters`](@ref) (`nothing` skips the ionosphere); the
troposphere uses the blind Saastamoinen zenith delays mapped by the Niell mapping
functions, whose seasonal term takes the day of year `doy` (see
[`tropospheric_delay`](@ref)), unless `enable_tropospheric_correction` is `false`.

A single corrected solve is enough: the delays depend on position only through the
satellite elevation/azimuth (and, for the troposphere, the user height), and
∂delay/∂position is negligible over the metre-level position uncertainty (a 15 m
shift moves the elevation by ~1e-5°), so delays predicted at a nearby position are
accurate to well under a millimetre — no iterate-to-convergence needed. The user
geodetic coordinates and the ENU transform depend only on `user_pos`, so they are
computed once and reused across satellites.
"""
Base.@nospecializeinfer function predict_atmospheric_delays(
    user_pos,
    measurements,
    @nospecialize(parameters),
    reference_time,
    doy,
    enable_tropospheric_correction,
)::Vector{Float64}
    # Taken `@nospecialize`d, so that with a single method a caller holding `parameters`
    # as `Any` resolves statically to this one unspecialised body, where an `Any`
    # argument to a specialising function would be a dynamic dispatch (which a
    # `juliac --trim` build rejects). The narrowing to each model's concrete type is
    # spelled out here instead; the solver itself does not come through here, it passes
    # an `IonosphericModel` straight to `predict_atmospheric_delays!`.
    #
    # The return annotation is load-bearing too: without it the `Any` argument makes the
    # return type `Any`.
    predict(model) = predict_atmospheric_delays!(
        Vector{Float64}(undef, length(measurements)), user_pos, measurements,
        IonosphericModel(model), reference_time, doy, enable_tropospheric_correction)
    parameters isa KlobucharParams && return predict(parameters)
    parameters isa BeiDouKlobucharParams && return predict(parameters)
    parameters isa NTCMGParams && return predict(parameters)
    parameters isa BDGIMParams && return predict(parameters)
    isnothing(parameters) && return predict(nothing)
    throw(ArgumentError("not a set of ionospheric parameters: $(typeof(parameters))"))
end

"""
    predict_atmospheric_delays!(delays, user_pos, measurements, model::IonosphericModel,
                                reference_time, doy, enable_tropospheric_correction)
        -> delays

[`predict_atmospheric_delays`](@ref) into `delays` (resized to the satellite count), with
the ionospheric parameters carried by an [`IonosphericModel`](@ref). Allocates nothing.
"""
function predict_atmospheric_delays!(
    delays,
    user_pos,
    measurements,
    model,
    reference_time,
    doy,
    enable_tropospheric_correction,
)
    # The per-satellite loop specialises once per ionospheric model, behind this barrier:
    # the field is a five-member `Union`, one more than Julia splits by itself, so the
    # split is spelled out.
    predict(parameters) = _predict_atmospheric_delays!(delays, user_pos, measurements,
        parameters, reference_time, doy, enable_tropospheric_correction)
    parameters = model.parameters
    parameters isa KlobucharParams && return predict(parameters)
    parameters isa BeiDouKlobucharParams && return predict(parameters)
    parameters isa NTCMGParams && return predict(parameters)
    parameters isa BDGIMParams && return predict(parameters)
    return predict(nothing)
end

function _predict_atmospheric_delays!(
    delays,
    user_pos,
    measurements,
    parameters,
    reference_time,
    doy,
    enable_tropospheric_correction,
)
    delays = resize!(delays, length(measurements))
    user_ecef = ECEF(user_pos[1], user_pos[2], user_pos[3])
    user_lla = LLAfromECEF(wgs84)(user_ecef)
    enu_from_ecef = ENUfromECEF(user_ecef, wgs84)
    for (j, measurement) in enumerate(measurements)
        elevation, azimuth = calc_elevation_and_azimuth(enu_from_ecef, measurement.position)
        iono = ionospheric_delay(
            parameters,
            measurement.center_frequency,
            elevation,
            azimuth,
            user_lla,
            reference_time,
        )
        tropo =
            enable_tropospheric_correction ? tropospheric_delay(elevation, user_lla, doy) :
            0.0
        delays[j] = iono + tropo
    end
    delays
end
