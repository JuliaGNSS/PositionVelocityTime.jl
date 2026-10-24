# ===========================================================================
#  Ionospheric correction — constellation-wide model selection
#
#  The broadcast ionospheric coefficients are global to a GNSS (every GPS
#  satellite broadcasts the same Klobuchar α/β; every Galileo satellite the same
#  NTCM-G Effective Ionisation Level a_i0..a_i2; every BDS-3 satellite the same
#  BDGIM α_1..α_9). So rather than choosing a model per satellite, a single model
#  is chosen for the whole solve from whatever the healthy decoders have
#  delivered, and applied to *all* satellites, in this fixed order of preference:
#
#      NTCM-G (Galileo) > BDGIM (BDS-3) > Klobuchar (GPS) > Klobuchar (BeiDou)
#
#  with no coefficients decoded at all meaning no correction (0 m for every sat).
#  The order is by model class, not by constellation. NTCM-G and BDGIM are both
#  global TEC models — a fitted empirical TEC field and a spherical-harmonic
#  expansion of one — whereas Klobuchar is an eight-coefficient half-cosine
#  approximation that IS-GPS-200 only undertakes to remove ~50 % of the delay
#  with; that is the substantive step, BDGIM ahead of either Klobuchar set. The
#  remaining rungs are tie-breaks rather than claims: NTCM-G stays first so a
#  Galileo-bearing epoch behaves exactly as it did before BDGIM existed, and GPS
#  stays ahead of BeiDou among the Klobuchar sets as it already did. What matters
#  about all four is that they are a fixed total order, so the chosen model does
#  not flip with the order the receiver happens to hand satellites over in.
#
#  Klobuchar reaches this package from two constellations: GPS on every civil
#  signal, and BeiDou on the legacy D1/D2 message — as BeiDou's own variant of
#  the algorithm, referenced to B1I rather than L1 (see
#  `BeiDouKlobucharParameters`). BDS-3's B-CNAV messages (B1C/B2a/B2b) instead
#  broadcast BDGIM, so a BDS-3-only epoch is corrected by BDGIM and no longer
#  needs a B1I/B3I or GPS satellite to be tracked alongside it.
#
#  Only data actually decoded from the navigation message is used; there are no
#  user-supplied fallback coefficients. Whichever model is chosen is applied across
#  constellations and frequency bands (L1/E1, L2, L5/E5a, E6, …), not just the band
#  it was broadcast on, but the two model families reach a given satellite's carrier
#  differently. The two Klobuchar variants return a delay *at one reference carrier*
#  (L1 for the GPS model, B1I for the BeiDou one) which is then rescaled by 1/f² to the satellite's
#  actual carrier. NTCM-G and BDGIM return TEC, which is frequency-free, and the
#  satellite's own carrier frequency enters the metre conversion directly — no
#  reference frequency exists to rescale from.
# ===========================================================================

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

The eight coefficients of BeiDou's own Klobuchar variant, decoded from the legacy
D1/D2 message (B1I/B3I), in BDS-SIS-ICD-B1I-3.0 §5.2.4.7 SI units (seconds and
seconds·π⁻ⁿ — per semicircle, exactly like the GPS set). A type of its own rather
than a [`KlobucharParameters`](@ref) because the ICD prescribes a different algorithm,
not merely different numbers — see [`beidou_klobuchar_group_delay`](@ref) — and
defines the resulting delay at the B1I carrier (1561.098 MHz) rather than L1.
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

The names are the ICD's own `α_1 … α_9` (BDS-SIS-ICD-B1C-1.0 Table 7-10);
GNSSDecoder spells them `α_bdgim_1 … α_bdgim_9` on its BeiDou containers to keep
them apart from the legacy D1/D2 Klobuchar `α_0 … α_3`, and already applies each
one's scale factor and sign — including `α_5`'s negative scale factor of −2⁻³ —
so they arrive here as plain signed TECU.

`week` is needed for the same reason `NTCMGParameters` carries one: the model's
time argument is a date (a Modified Julian Date), which a time of week alone
cannot supply. Unlike the Klobuchar set there is no reference frequency, because
BDGIM yields TEC rather than a delay at a particular carrier.
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

Klobuchar α/β decoded from a GPS navigation message (LNAV `GPSL1CAData`, CNAV
`GPSCNAVData` on L5/L2C, or CNAV-2 `GPSL1C_DData`) or from the BeiDou legacy
D1/D2 message (`BeiDouDNAVData` on B1I/B3I), or `nothing` if
they have not been broadcast yet or the decoder carries no Klobuchar set. The same
single-frequency Klobuchar model is broadcast on all GPS civil signals; BeiDou
broadcasts its own set, returned as a [`BeiDouKlobucharParameters`](@ref) because its
ICD prescribes its own variant of the algorithm, not just its own numbers.
"""
get_klobuchar_parameters(decoder) = nothing
function get_klobuchar_parameters(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGPSData},
)
    d = decoder.data
    # All eight coefficients must be present: they are decoded together (subframe 4,
    # page 18), but guard each so a partially-populated decoder returns `nothing`
    # rather than throwing when a `nothing` hits a `Float64` field.
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return KlobucharParameters(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

# The legacy D1/D2 message broadcasts an eight-coefficient Klobuchar-style set,
# but BDS-SIS-ICD-B1I-3.0 §5.2.4.7 prescribes BeiDou's own variant of the
# algorithm, not a restatement of IS-GPS-200 — exact pierce-point geometry,
# geographic rather than geomagnetic latitude, a period clamped from above too,
# and the delay referenced to B1I rather than L1 — so the set gets a type of its
# own (see [`BeiDouKlobucharParameters`](@ref) and
# [`beidou_klobuchar_group_delay`](@ref)). The BDS-3 B-CNAV messages instead
# broadcast BDGIM, reached through `get_bdgim_parameters` below — those containers
# therefore report no Klobuchar set.
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
message (I/NAV on E1-B or E5b, or F/NAV on E5a — all broadcast `a_i0…a_i2`), or
`nothing` if they (or the week number) have not been decoded yet or the decoder
carries no such set. Galileo E6-B is in the latter group: C/NAV broadcasts no
ionospheric coefficients, so it falls to the generic method rather than being
dispatched here on a field it does not have.
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

The nine BDGIM coefficients decoded from a BDS-3 B-CNAV message — B-CNAV1 on B1C,
B-CNAV2 on B2a, B-CNAV3 on B2b, all three broadcasting the same set — or
`nothing` if they, or the week number the model needs to
date its epoch, have not been decoded yet, or the decoder carries no such set. The
legacy D1/D2 message (`BeiDouDNAVData`, on B1I/B3I) is in the latter group: it
broadcasts Klobuchar instead, so it falls to this generic method rather than being
dispatched on fields it does not have.
"""
get_bdgim_parameters(decoder) = nothing

# Nine coefficients broadcast identically by all three B-CNAV messages — B-CNAV1
# on B1C (BDS-SIS-ICD-B1C-1.0 §7.8, subframe 3 page type 1), B-CNAV2 on B2a
# (BDS-SIS-ICD-B2a-1.0 §7.8, message type 30) and B-CNAV3 on B2b
# (BDS-SIS-ICD-B2b-1.0 §7.7, message type 30) — so one method on the union covers
# all three rather than three identical ones. GNSSDecoder parses them with the
# shared `beidou_bdgim_block`, which applies each scale factor and sign, so
# nothing is rescaled here.
#
# `WN` is required alongside the coefficients: BDGIM's time argument is a Modified
# Julian Date (the sun's mean longitude and the two-hourly prediction epoch), which
# a time of week alone cannot date. It is decoded in a different subframe/message
# type from the coefficients on every one of the three signals, so a decoder that
# has one without the other is the ordinary transient state, not a corrupt one —
# hence the guard rather than an unwrap.
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

# The four candidate coefficient sets one epoch can yield. A struct of `Union`-with-
# `nothing` fields rather than a plain tuple of whatever was found: the accumulation
# changes which sets are present as it goes, and a tuple would change its *type* with
# them — a `Union` of tuple types that is boxed (an allocation per group, per epoch) when
# the per-group results are collected. The struct is one concrete type whatever it holds,
# and its `isbits` `Union` fields are stored inline. The accumulation stays inside the
# per-group loop of `collect_group!`, where the decoder type is concrete and every
# accessor below is statically dispatched.
struct IonosphericCandidates
    gps_klobuchar::Union{Nothing,KlobucharParameters}
    beidou_klobuchar::Union{Nothing,BeiDouKlobucharParameters}
    bdgim::Union{Nothing,BDGIMParameters}
    ntcm_g::Union{Nothing,NTCMGParameters}
end

const NO_IONOSPHERIC_CANDIDATES = IonosphericCandidates(nothing, nothing, nothing, nothing)

# Fold one decoder into the candidates. Each set is global to its constellation, so
# the first decoder that carries one wins and later ones are not even asked.
function update_ionospheric_candidates(candidates, decoder)
    (; gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g) = candidates
    if decoder.data isa GNSSDecoder.AbstractBeiDouData
        # The two BeiDou message families broadcast different models: the legacy
        # D1/D2 one Klobuchar, the BDS-3 B-CNAV ones BDGIM. Both accessors are
        # asked; each answers `nothing` for the family that is not its own.
        isnothing(beidou_klobuchar) &&
            (beidou_klobuchar = get_klobuchar_parameters(decoder))
        isnothing(bdgim) && (bdgim = get_bdgim_parameters(decoder))
    else
        isnothing(gps_klobuchar) && (gps_klobuchar = get_klobuchar_parameters(decoder))
    end
    isnothing(ntcm_g) && (ntcm_g = get_ntcm_g_parameters(decoder))
    IonosphericCandidates(gps_klobuchar, beidou_klobuchar, bdgim, ntcm_g)
end

# Combine two groups' candidates, keeping the earlier group's set where both have one
# — the same first-one-wins rule the within-group fold uses, extended across groups in
# group order.
merge_ionospheric_candidates(a, b) = IonosphericCandidates(
    first_available(a.gps_klobuchar, b.gps_klobuchar),
    first_available(a.beidou_klobuchar, b.beidou_klobuchar),
    first_available(a.bdgim, b.bdgim),
    first_available(a.ntcm_g, b.ntcm_g),
)
first_available(x, y) = isnothing(x) ? y : x

# Fold a tuple of per-group candidates, earliest group first. Recursive over the tuple
# rather than `reduce`d, so that it unrolls over the groups.
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

The ionospheric `parameters` of one solve (see [`select_ionospheric_parameters`](@ref)),
held in a single concrete type. The parameters themselves are a `Union` of four coefficient
sets and `nothing`; the solver takes this wrapper instead so that it compiles once whichever
model an epoch selected, and so that passing the model costs no allocation — a `Union` value
handed to an unspecialised argument is boxed, whereas a field of a concrete struct is stored
inline.
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

The two global TEC models come first because they are the better ones: BDGIM is a
spherical-harmonic TEC expansion, whereas Klobuchar is the eight-coefficient
half-cosine approximation IS-GPS-200 only undertakes to remove ~50 % of the delay
with. The rungs among equals are tie-breaks rather than claims — NTCM-G stays
ahead of BDGIM so a Galileo-bearing epoch behaves exactly as it did before, and
GPS stays ahead of BeiDou among the Klobuchar sets as it already did — but every
rung is fixed rather than data-dependent, so the chosen model does not flip with
the order the receiver happens to hand satellites over in. The coefficients are
global to a constellation, so the first decoder that carries each set is used —
"first" in group order × within-group order, the same order everything else about the
epoch is ordered in.

Every satellite of `groups` is scanned, healthy or not — the caller decides which
satellites to offer, exactly as it decides which to solve with. [`calc_pvt`](@ref) does
not call this at all: [`collect_measurements`](@ref) folds the same scan into its own
pass over the satellites, which is the one place they are already being visited, and
there the scan sees only the satellites that clear the health gate. This is the
standalone form, for a consumer running its own estimator over the same measurement
model.
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

The line of sight is given by the satellite `elevation`/`azimuth` (radians) and
the user geodetic position `lla` (a `Geodesy.LLA`); `time_of_week` is the measurement's
system time of week in seconds. Which constellation's scale that is does not matter:
GPST and GST coincide to within the sub-microsecond GGTO, and BDT trails them by 14 s,
which moves Klobuchar's diurnal phase by 0.016 % of a day, NTCM-G's universal time
by the same, and BDGIM's solar longitude by 0.06 milliradians — orders of magnitude
below any of the three models' own error. The geometry is taken precomputed — and
shared across satellites and with [`tropospheric_delay`](@ref) — so a whole-epoch
correction does the user geodetic conversion only once. Derive the geometry from
ECEF with `LLAfromECEF(wgs84)(user)` and
[`calc_elevation_and_azimuth`](@ref)`(ENUfromECEF(user, wgs84), sat)`.

`center_frequency` is the ranging signal's carrier in Hz as a plain `Float64`
(`ustrip(Hz, get_center_frequency(signal))`, precomputed on every
[`SatelliteMeasurement`](@ref)) — every model needs the carrier and nothing else about
the signal, which is why the row carries the frequency rather than the signal.
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
    # The Klobuchar broadcast coefficients define the group delay at the GPS L1
    # frequency (IS-GPS-200). The ionospheric delay scales as 1/f², so rescale it
    # to this satellite's actual carrier frequency.
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
    # The coefficients define the group delay along the B1I propagation path (the
    # ICD's I_B1I) — 1561.098 MHz, not L1. The delay scales as 1/f², so rescale
    # from B1I to this satellite's actual carrier.
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
    # Eq. 7-6: T_ion [m] = M_F · 40.28e16/f² · VTEC, i.e. 40.28e16/f² · STEC with STEC
    # in TECU. `f` is this satellite's own carrier, exactly as the ICD says ("the
    # carrier frequency of the current signal"); there is no reference frequency to
    # rescale from, unlike the Klobuchar branch. Note the constant is BDGIM's own
    # 40.28, not NTCM-G's 40.3 — a 0.05 % difference, kept ICD-exact per model.
    return 40.28 / f^2 * stec * 1.0e16
end
