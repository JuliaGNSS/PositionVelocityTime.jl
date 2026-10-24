# Ionospheric correction from the broadcast coefficients alone: one model for the whole
# solve, chosen by `select_ionospheric_parameters` (see there for the order and why), and
# applied to every satellite and band. The Klobuchar variants give a delay at a reference
# carrier (L1, B1I) rescaled by 1/f²; NTCM-G and BDGIM give TEC, converted at the
# satellite's own carrier.

"""
    KlobucharParameters(α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3)

The eight Klobuchar ionospheric coefficients decoded from a GPS navigation message
(LNAV, CNAV or CNAV-2 — every civil signal broadcasts the same set), in IS-GPS-200 SI
units (seconds and seconds·semicircle⁻ⁿ). The field names mirror GNSSDecoder's GPS
data containers (`α_0…α_3`, `β_0…β_3`).
"""
struct KlobucharParameters
    α_0::Float64
    α_1::Float64
    α_2::Float64
    α_3::Float64
    β_0::Float64
    β_1::Float64
    β_2::Float64
    β_3::Float64
end

"""
    BeiDouKlobucharParameters(α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3)

The eight coefficients of BeiDou's Klobuchar variant, decoded from the legacy D1/D2
message (B1I/B3I), in BDS-SIS-ICD-B1I-3.0 §5.2.4.7 SI units (seconds and seconds·π⁻ⁿ,
per semicircle like the GPS set). A type of its own because the ICD prescribes a
different algorithm (`beidou_klobuchar_group_delay`), with the delay defined at
B1I (1561.098 MHz) rather than L1.
"""
struct BeiDouKlobucharParameters
    α_0::Float64
    α_1::Float64
    α_2::Float64
    α_3::Float64
    β_0::Float64
    β_1::Float64
    β_2::Float64
    β_3::Float64
end

"""
    NTCMGParameters(a_i0, a_i1, a_i2, week::Integer)

The broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`
(sfu, sfu/deg, sfu/deg²) decoded from an E1B navigation message, together with the
GST week number needed to derive the day of year and universal time for NTCM-G.
"""
struct NTCMGParameters
    a_i0::Float64
    a_i1::Float64
    a_i2::Float64
    week::Int
end

"""
    BDGIMParameters(α_1, …, α_9, week::Integer)

The nine broadcast BDGIM coefficients decoded from a BDS-3 B-CNAV message
(B1C/B2a/B2b), in the ICD's units of TECU, together with the BDT week number.

The names are the ICD's `α_1 … α_9` (BDS-SIS-ICD-B1C-1.0 Table 7-10), GNSSDecoder's
`α_bdgim_1 … α_bdgim_9`, which arrive with each scale factor and sign applied
(including `α_5`'s −2⁻³). `week` dates the epoch: the model's time argument is a
Modified Julian Date.
"""
struct BDGIMParameters
    α_1::Float64
    α_2::Float64
    α_3::Float64
    α_4::Float64
    α_5::Float64
    α_6::Float64
    α_7::Float64
    α_8::Float64
    α_9::Float64
    week::Int
end

"""
    get_klobuchar_parameters(decoder)
        -> Union{KlobucharParameters,BeiDouKlobucharParameters,Nothing}

Klobuchar α/β decoded from a GPS navigation message (LNAV, CNAV or CNAV-2), or, as a
[`BeiDouKlobucharParameters`](@ref), from the BeiDou legacy D1/D2 message; `nothing` if
not decoded yet or the decoder carries no Klobuchar set.
"""
get_klobuchar_parameters(decoder) = nothing
function get_klobuchar_parameters(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGPSData},
)
    d = decoder.data
    # Decoded together (subframe 4, page 18), but a partly populated decoder returns
    # `nothing` rather than throwing.
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return KlobucharParameters(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

# BeiDou's variant (see `BeiDouKlobucharParameters`). The B-CNAV messages broadcast
# BDGIM instead and report no Klobuchar set.
function get_klobuchar_parameters(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
)
    d = decoder.data
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return BeiDouKlobucharParameters(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

"""
    get_ntcm_g_parameters(decoder) -> Union{NTCMGParameters,Nothing}

NTCM-G Effective Ionisation Level coefficients decoded from a Galileo navigation
message (I/NAV on E1-B or E5b, F/NAV on E5a), or `nothing` if they or the week number
are not decoded yet or the decoder carries none (E6-B's C/NAV broadcasts none).
"""
get_ntcm_g_parameters(decoder) = nothing
function get_ntcm_g_parameters(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGalileoEphemerisData},
)
    d = decoder.data
    (isnothing(d.a_i0) || isnothing(d.a_i1) || isnothing(d.a_i2) || isnothing(d.WN)) &&
        return nothing
    return NTCMGParameters(d.a_i0, d.a_i1, d.a_i2, d.WN)
end

"""
    get_bdgim_parameters(decoder) -> Union{BDGIMParameters,Nothing}

The nine BDGIM coefficients decoded from a BDS-3 B-CNAV message (B1C, B2a or B2b, all
broadcasting the same set), or `nothing` if they or the week number are not decoded yet
or the decoder carries none (the legacy D1/D2 message broadcasts Klobuchar).
"""
get_bdgim_parameters(decoder) = nothing

# One method for B-CNAV1 (BDS-SIS-ICD-B1C-1.0 §7.8, subframe 3 page type 1), B-CNAV2
# (BDS-SIS-ICD-B2a-1.0 §7.8, MT30) and B-CNAV3 (BDS-SIS-ICD-B2b-1.0 §7.7, MT30). `WN`
# comes in another subframe/message type than the coefficients, so having one without
# the other is an ordinary transient state: hence the guard rather than an unwrap.
function get_bdgim_parameters(
    decoder::GNSSDecoder.GNSSDecoderState{<:AbstractBeiDouCNAVData},
)
    d = decoder.data
    any(
        isnothing,
        (
            d.α_bdgim_1,
            d.α_bdgim_2,
            d.α_bdgim_3,
            d.α_bdgim_4,
            d.α_bdgim_5,
            d.α_bdgim_6,
            d.α_bdgim_7,
            d.α_bdgim_8,
            d.α_bdgim_9,
            d.WN,
        ),
    ) && return nothing
    return BDGIMParameters(
        d.α_bdgim_1,
        d.α_bdgim_2,
        d.α_bdgim_3,
        d.α_bdgim_4,
        d.α_bdgim_5,
        d.α_bdgim_6,
        d.α_bdgim_7,
        d.α_bdgim_8,
        d.α_bdgim_9,
        d.WN,
    )
end

# The four candidate coefficient sets one epoch can yield. A struct rather than a tuple
# of what was found: a tuple's type changes with its contents, and the resulting `Union`
# of tuple types is boxed when the per-group results are collected; this struct is one
# concrete type with inline `isbits` fields.
struct IonosphericCandidates
    gps_klobuchar::Union{Nothing,KlobucharParameters}
    beidou_klobuchar::Union{Nothing,BeiDouKlobucharParameters}
    bdgim::Union{Nothing,BDGIMParameters}
    ntcm_g::Union{Nothing,NTCMGParameters}
end

const NO_IONOSPHERIC_CANDIDATES = IonosphericCandidates(nothing, nothing, nothing, nothing)

# Fold one decoder into the candidates. Each set is constellation-wide, so the first
# decoder that carries one wins.
function update_ionospheric_candidates(candidates, decoder)
    (; gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g) = candidates
    if decoder.data isa GNSSDecoder.AbstractBeiDouData
        # D1/D2 broadcasts Klobuchar, B-CNAV BDGIM; each accessor answers `nothing`
        # for the other family.
        isnothing(beidou_klobuchar) &&
            (beidou_klobuchar = get_klobuchar_parameters(decoder))
        isnothing(bdgim) && (bdgim = get_bdgim_parameters(decoder))
    else
        isnothing(gps_klobuchar) && (gps_klobuchar = get_klobuchar_parameters(decoder))
    end
    isnothing(ntcm_g) && (ntcm_g = get_ntcm_g_parameters(decoder))
    IonosphericCandidates(gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g)
end

# Combine two groups' candidates, the earlier group's set winning.
merge_ionospheric_candidates(a, b) = IonosphericCandidates(
    first_available(a.gps_klobuchar, b.gps_klobuchar),
    first_available(a.beidou_klobuchar, b.beidou_klobuchar),
    first_available(a.bdgim, b.bdgim),
    first_available(a.ntcm_g, b.ntcm_g),
)
first_available(x, y) = isnothing(x) ? y : x

# Fold the per-group candidates, earliest first; recursive so that it unrolls.
merge_all_ionospheric_candidates(::Tuple{}) = NO_IONOSPHERIC_CANDIDATES
merge_all_ionospheric_candidates(candidates::Tuple) = merge_ionospheric_candidates(
    first(candidates),
    merge_all_ionospheric_candidates(Base.tail(candidates)),
)

# The order of preference documented on `select_ionospheric_parameters`, top rung first.
function select_from_ionospheric_candidates(candidates)
    (; gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g) = candidates
    !isnothing(ntcm_g) && return ntcm_g
    !isnothing(bdgim) && return bdgim
    !isnothing(gps_klobuchar) && return gps_klobuchar
    return beidou_klobuchar          # BeiDouKlobucharParameters, or nothing if neither
end

"""
    IonosphericModel(parameters)

The ionospheric `parameters` of one solve (see [`select_ionospheric_parameters`](@ref)) in
one concrete type, so that the solver compiles once whichever model was selected and
passing it allocates nothing: the bare `Union` would be boxed.
"""
struct IonosphericModel
    parameters::Union{
        Nothing,
        KlobucharParameters,
        BeiDouKlobucharParameters,
        NTCMGParameters,
        BDGIMParameters,
    }
end

"""
    select_ionospheric_parameters(groups)
        -> Union{KlobucharParameters,BeiDouKlobucharParameters,NTCMGParameters,
                 BDGIMParameters,Nothing}

Scan the satellite decoders of one epoch's [`SignalGroups`](@ref) and pick the
single set of ionospheric parameters to correct the whole solve with, in a fixed order of
preference:

 1. NTCM-G, if Galileo Effective Ionisation Level coefficients were decoded.
 2. BDGIM, if a BDS-3 B-CNAV (B1C/B2a/B2b) coefficient set was decoded.
 3. Klobuchar from GPS, if the GPS α/β were decoded.
 4. BeiDou's own Klobuchar variant from its legacy D1/D2 message (B1I/B3I).
 5. `nothing` — no coefficients at all, so no correction.

The two global TEC models (NTCM-G, and BDGIM's spherical-harmonic expansion) come
first because they are the better ones: Klobuchar is a half-cosine approximation that
IS-GPS-200 only undertakes to remove ~50 % of the delay with. Among equals the order is
a tie-break (NTCM-G kept ahead of BDGIM, GPS ahead of BeiDou), but it is fixed, so the
choice does not depend on the order the satellites arrive in. Each set is global to its
constellation, so the first decoder carrying it is used, in row order (see
[`calc_pvt`](@ref)).

This standalone form scans every satellite of `groups`, healthy or not;
[`collect_measurements`](@ref) folds the same scan into the collection pass, where only
satellites that clear the health gate are seen.
"""
function select_ionospheric_parameters(groups)
    # Over the groups' `Tuple`, not the `NamedTuple`, as in `collect_measurements`.
    candidates = map(values(_normalize_signal_groups(groups))) do group
        group_candidates = NO_IONOSPHERIC_CANDIDATES
        for state in group.satellites
            group_candidates =
                update_ionospheric_candidates(group_candidates, state.decoder)
        end
        group_candidates
    end
    select_from_ionospheric_candidates(merge_all_ionospheric_candidates(candidates))
end

"""
    ionospheric_delay(parameters, center_frequency, elevation, azimuth, lla, time_of_week)
        -> Float64

Slant ionospheric group delay in metres for one satellite, using the
constellation-wide `parameters` returned by [`select_ionospheric_parameters`](@ref):

- `::Nothing` → `0.0` (no coefficients were decoded).
- [`KlobucharParameters`](@ref) → Klobuchar model (IS-GPS-200, Fig. 20-4).
- [`BeiDouKlobucharParameters`](@ref) → BeiDou's Klobuchar variant
  (BDS-SIS-ICD-B1I-3.0 §5.2.4.7).
- [`NTCMGParameters`](@ref) → NTCM-G model.
- [`BDGIMParameters`](@ref) → BDGIM (BDS-SIS-ICD-B1C-1.0 §7.8.2).

`elevation`/`azimuth` (radians) give the line of sight from the user geodetic position
`lla` (a `Geodesy.LLA`), precomputed so that an epoch shares one geodetic conversion with
[`tropospheric_delay`](@ref): `LLAfromECEF(wgs84)(user)` and
[`calc_elevation_and_azimuth`](@ref)`(ENUfromECEF(user, wgs84), sat)`. `time_of_week` is
in seconds on any constellation's scale: BDT's 14 s moves every model's time argument by
orders of magnitude less than its own error. `center_frequency` is the carrier in Hz, as
a plain `Float64` (as on [`SatelliteMeasurement`](@ref)).
"""
ionospheric_delay(::Nothing, center_frequency, elevation, azimuth, lla, time_of_week) = 0.0

function ionospheric_delay(
    parameters::KlobucharParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    time_of_week,
)
    # IS-GPS-200 works in semicircles: lat/lon in deg/180, elevation/azimuth in rad/π.
    l1_seconds = klobuchar_group_delay(
        lla.lat / 180,
        lla.lon / 180,
        elevation / π,
        azimuth / π,
        time_of_week,
        (parameters.α_0, parameters.α_1, parameters.α_2, parameters.α_3),
        (parameters.β_0, parameters.β_1, parameters.β_2, parameters.β_3),
    )
    # The delay is defined at L1 (IS-GPS-200); rescale by 1/f² to this carrier.
    return SPEED_OF_LIGHT * l1_seconds *
           (ustrip(Hz, get_center_frequency(GPSL1CA)) / center_frequency)^2
end

function ionospheric_delay(
    parameters::BeiDouKlobucharParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    time_of_week,
)
    # BDS-SIS-ICD-B1I-3.0 §5.2.4.7 works in radians, unlike IS-GPS-200's semicircles.
    b1i_seconds = beidou_klobuchar_group_delay(
        deg2rad(lla.lat),
        deg2rad(lla.lon),
        elevation,
        azimuth,
        time_of_week,
        (parameters.α_0, parameters.α_1, parameters.α_2, parameters.α_3),
        (parameters.β_0, parameters.β_1, parameters.β_2, parameters.β_3),
    )
    # The delay is defined at B1I (the ICD's I_B1I); rescale by 1/f² to this carrier.
    return SPEED_OF_LIGHT * b1i_seconds *
           (ustrip(Hz, get_center_frequency(GNSSSignals.BeiDouB1I)) / center_frequency)^2
end

function ionospheric_delay(
    parameters::NTCMGParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    time_of_week,
)
    doy, ut = _galileo_doy_and_ut(parameters.week, time_of_week)
    (; a_i0, a_i1, a_i2) = parameters
    stec = ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2) # TECU
    f = center_frequency
    # Eq. 1: group delay [m] = 40.3 / f² · STEC, with STEC in electrons/m² (1 TECU = 1e16).
    return 40.3 / f^2 * stec * 1.0e16
end

function ionospheric_delay(
    parameters::BDGIMParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    time_of_week,
)
    mjd = _beidou_modified_julian_date(parameters.week, time_of_week)
    (; α_1, α_2, α_3, α_4, α_5, α_6, α_7, α_8, α_9) = parameters
    α = (α_1, α_2, α_3, α_4, α_5, α_6, α_7, α_8, α_9)
    stec = bdgim_stec(elevation, azimuth, lla, mjd, α)   # TECU
    f = center_frequency
    # Eq. 7-6: T_ion [m] = 40.28e16/f² · STEC (TECU), `f` this signal's own carrier.
    # BDGIM's 40.28, not NTCM-G's 40.3: each model keeps its ICD's constant.
    return 40.28 / f^2 * stec * 1.0e16
end
