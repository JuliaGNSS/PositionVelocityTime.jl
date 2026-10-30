# Ionospheric correction from the broadcast coefficients alone: one model for the whole
# solve, chosen by `select_ionospheric_parameters` (see there for the order and why), and
# applied to every satellite and band. The Klobuchar variants give a delay at a reference
# carrier (L1, B1I) rescaled by 1/f²; NTCM-G and BDGIM give TEC, converted at the
# satellite's own carrier. The coefficient sets and their readers are GNSSDecoder's
# (`KlobucharParameters`, …, `get_klobuchar_parameters`, …); the models are here.

# The four candidate coefficient sets one epoch can yield. A struct rather than a tuple
# of what was found: a tuple's type changes with its contents, and the resulting `Union`
# of tuple types is boxed when the per-group results are collected; this struct is one
# concrete type with inline `isbits` fields.
struct IonosphericCandidates
    gps_klobuchar::Union{Nothing,KlobucharParameters}
    beidou_klobuchar::Union{Nothing,BeiDouKlobucharParameters}
    bdgim::Union{Nothing,BDGIMParameters}
    ntcm_g::Union{Nothing,NeQuickGParameters}
end

const NO_IONOSPHERIC_CANDIDATES = IonosphericCandidates(nothing, nothing, nothing, nothing)

# Fold one decoder into the candidates. Each set is constellation-wide, so the first
# decoder that carries one wins.
function update_ionospheric_candidates(candidates, decoder)
    (; gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g) = candidates
    # GPS and BeiDou D1/D2 both broadcast a Klobuchar set, told apart by its type. Every
    # accessor answers `nothing` on a decoder that carries no such set.
    klobuchar = get_klobuchar_parameters(decoder)
    if klobuchar isa BeiDouKlobucharParameters
        isnothing(beidou_klobuchar) && (beidou_klobuchar = klobuchar)
    elseif klobuchar isa KlobucharParameters
        isnothing(gps_klobuchar) && (gps_klobuchar = klobuchar)
    end
    isnothing(bdgim) && (bdgim = get_bdgim_parameters(decoder))
    isnothing(ntcm_g) && (ntcm_g = get_nequick_g_parameters(decoder))
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
        NeQuickGParameters,
        BDGIMParameters,
    }
end

"""
    select_ionospheric_parameters(groups)
        -> Union{KlobucharParameters,BeiDouKlobucharParameters,NeQuickGParameters,
                 BDGIMParameters,Nothing}

Scan the satellite decoders of one epoch's signal groups (see [`calc_pvt`](@ref)) and pick the
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
constellation, so the first decoder carrying it is used, in the order of the measurements
(see [`calc_pvt`](@ref)).

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
    ionospheric_delay(parameters, center_frequency, elevation, azimuth, lla, epoch::TAITime)
        -> Float64

Slant ionospheric group delay in metres for one satellite, using the
constellation-wide `parameters` returned by [`select_ionospheric_parameters`](@ref):

- `::Nothing` → `0.0` (no coefficients were decoded).
- `KlobucharParameters` → Klobuchar model (IS-GPS-200, Fig. 20-4).
- `BeiDouKlobucharParameters` → BeiDou's Klobuchar variant (BDS-SIS-ICD-B1I-3.0
  §5.2.4.7).
- `NeQuickGParameters` → NTCM-G model, which takes the Effective Ionisation Level
  coefficients Galileo broadcasts for NeQuick-G.
- `BDGIMParameters` → BDGIM (BDS-SIS-ICD-B1C-1.0 §7.8.2).

`elevation`/`azimuth` (radians) give the line of sight from the user geodetic position
`lla` (a `Geodesy.LLA`), precomputed so that an epoch shares one geodetic conversion with
[`tropospheric_delay`](@ref): `LLAfromECEF(wgs84)(user)` and
[`calc_elevation_and_azimuth`](@ref)`(ENUfromECEF(user, wgs84), sat)`. `epoch` is the
TAI epoch the delay is predicted for: the Klobuchar models take its time of week on GPS
Time or BeiDou Time, NTCM-G its day of year and universal time and BDGIM its Modified
Julian Date, the latter two read off the TAI calendar label (37 s from UTC, orders of
magnitude below the models' error). `center_frequency` is the carrier in Hz, as a plain
`Float64` (as on [`PreparedObservation`](@ref)).
"""
ionospheric_delay(::Nothing, center_frequency, elevation, azimuth, lla, epoch) = 0.0

function ionospheric_delay(
    p::KlobucharParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    epoch::TAITime,
)
    # IS-GPS-200 works in semicircles: lat/lon in deg/180, elevation/azimuth in rad/π.
    l1_seconds = klobuchar_group_delay(
        lla.lat / 180,
        lla.lon / 180,
        elevation / π,
        azimuth / π,
        get_time_of_week(GPST(), epoch),
        (p.α_0, p.α_1, p.α_2, p.α_3),
        (p.β_0, p.β_1, p.β_2, p.β_3),
    )
    # The delay is defined at L1 (IS-GPS-200); rescale by 1/f² to this carrier.
    return SPEED_OF_LIGHT * l1_seconds *
           (ustrip(Hz, get_center_frequency(GPSL1CA)) / center_frequency)^2
end

function ionospheric_delay(
    p::BeiDouKlobucharParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    epoch::TAITime,
)
    # BDS-SIS-ICD-B1I-3.0 §5.2.4.7 works in radians, unlike IS-GPS-200's semicircles.
    b1i_seconds = beidou_klobuchar_group_delay(
        deg2rad(lla.lat),
        deg2rad(lla.lon),
        elevation,
        azimuth,
        get_time_of_week(BDT(), epoch),
        (p.α_0, p.α_1, p.α_2, p.α_3),
        (p.β_0, p.β_1, p.β_2, p.β_3),
    )
    # The delay is defined at B1I (the ICD's I_B1I); rescale by 1/f² to this carrier.
    return SPEED_OF_LIGHT * b1i_seconds *
           (ustrip(Hz, get_center_frequency(GNSSSignals.BeiDouB1I)) / center_frequency)^2
end

function ionospheric_delay(
    p::NeQuickGParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    epoch::TAITime,
)
    # The day of year and universal time (h) off the TAI calendar label; its 37 s lead on
    # UTC is negligible here.
    t = DateTime(epoch)
    doy, ut = dayofyear(t), Dates.value(t - floor(t, Day)) / 3_600_000
    (; a_i0, a_i1, a_i2) = p
    stec = ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2) # TECU
    f = center_frequency
    # Eq. 1: group delay [m] = 40.3 / f² · STEC, with STEC in electrons/m² (1 TECU = 1e16).
    return 40.3 / f^2 * stec * 1.0e16
end

function ionospheric_delay(
    p::BDGIMParameters,
    center_frequency,
    elevation,
    azimuth,
    lla,
    epoch::TAITime,
)
    # The Modified Julian Date off the TAI calendar label; its 37 s lead on UTC is
    # negligible for a date.
    mjd = (epoch - TAITime(DateTime(1858, 11, 17))) / 86_400
    (; α_1, α_2, α_3, α_4, α_5, α_6, α_7, α_8, α_9) = p
    α = (α_1, α_2, α_3, α_4, α_5, α_6, α_7, α_8, α_9)
    stec = bdgim_stec(elevation, azimuth, lla, mjd, α)   # TECU
    f = center_frequency
    # Eq. 7-6: T_ion [m] = 40.28e16/f² · STEC (TECU), `f` this signal's own carrier.
    # BDGIM's 40.28, not NTCM-G's 40.3: each model keeps its ICD's constant.
    return 40.28 / f^2 * stec * 1.0e16
end
