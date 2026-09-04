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
#  `BeiDouKlobucharParams`). BDS-3's B-CNAV messages (B1C/B2a/B2b) instead
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
    KlobucharParams(α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3)

The eight Klobuchar ionospheric coefficients decoded from a GPS L1 navigation
message, in IS-GPS-200 SI units (seconds and seconds·semicircle⁻ⁿ). The field
names mirror `GNSSDecoder.GPSL1CAData` (`α_0…α_3`, `β_0…β_3`).
"""
struct KlobucharParams
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
    BeiDouKlobucharParams(α_0, α_1, α_2, α_3, β_0, β_1, β_2, β_3)

The eight coefficients of BeiDou's own Klobuchar variant, decoded from the legacy
D1/D2 message (B1I/B3I), in BDS-SIS-ICD-B1I-3.0 §5.2.4.7 SI units (seconds and
seconds·π⁻ⁿ — per semicircle, exactly like the GPS set). A type of its own rather
than a [`KlobucharParams`](@ref) because the ICD prescribes a different algorithm,
not merely different numbers — see [`beidou_klobuchar_group_delay`](@ref) — and
defines the resulting delay at the B1I carrier (1561.098 MHz) rather than L1.
"""
struct BeiDouKlobucharParams
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
    NTCMGParams(a_i0, a_i1, a_i2, week_number::Integer)

The broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`
(sfu, sfu/deg, sfu/deg²) decoded from an E1B navigation message, together with the
GST week number needed to derive the day of year and universal time for NTCM-G.
"""
struct NTCMGParams
    a_i0::Float64
    a_i1::Float64
    a_i2::Float64
    week_number::Int
end

"""
    BDGIMParams(α_1, …, α_9, week_number::Integer)

The nine broadcast BDGIM coefficients decoded from a BDS-3 B-CNAV message
(B1C/B2a/B2b), in the ICD's units of TECu, together with the BDT week number.

The names are the ICD's own `α_1 … α_9` (BDS-SIS-ICD-B1C-1.0 Table 7-10);
GNSSDecoder spells them `α_bdgim_1 … α_bdgim_9` on its BeiDou containers to keep
them apart from the legacy D1/D2 Klobuchar `α_0 … α_3`, and already applies each
one's scale factor and sign — including `α_5`'s negative scale factor of −2⁻³ —
so they arrive here as plain signed TECu.

`week_number` is needed for the same reason `NTCMGParams` carries one: the model's
time argument is a date (a Modified Julian Date), which a time of week alone
cannot supply. Unlike the Klobuchar set there is no reference frequency, because
BDGIM yields TEC rather than a delay at a particular carrier.
"""
struct BDGIMParams
    α_1::Float64
    α_2::Float64
    α_3::Float64
    α_4::Float64
    α_5::Float64
    α_6::Float64
    α_7::Float64
    α_8::Float64
    α_9::Float64
    week_number::Int
end

"""
    klobuchar_params(decoder) -> Union{KlobucharParams,BeiDouKlobucharParams,Nothing}

Klobuchar α/β decoded from a GPS navigation message (LNAV `GPSL1CAData`, CNAV
`GPSCNAVData` on L5/L2C, or CNAV-2 `GPSL1C_DData`) or from the BeiDou legacy
D1/D2 message (`BeiDouDNAVData` on B1I/B3I), or `nothing` if
they have not been broadcast yet or the decoder carries no Klobuchar set. The same
single-frequency Klobuchar model is broadcast on all GPS civil signals; BeiDou
broadcasts its own set, returned as a [`BeiDouKlobucharParams`](@ref) because its
ICD prescribes its own variant of the algorithm, not just its own numbers.
"""
klobuchar_params(decoder) = nothing
function klobuchar_params(decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGPSData})
    d = decoder.data
    # All eight coefficients must be present: they are decoded together (subframe 4,
    # page 18), but guard each so a partially-populated decoder returns `nothing`
    # rather than throwing when a `nothing` hits a `Float64` field.
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return KlobucharParams(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

# The legacy D1/D2 message broadcasts an eight-coefficient Klobuchar-style set,
# but BDS-SIS-ICD-B1I-3.0 §5.2.4.7 prescribes BeiDou's own variant of the
# algorithm, not a restatement of IS-GPS-200 — exact pierce-point geometry,
# geographic rather than geomagnetic latitude, a period clamped from above too,
# and the delay referenced to B1I rather than L1 — so the set gets a type of its
# own (see [`BeiDouKlobucharParams`](@ref) and
# [`beidou_klobuchar_group_delay`](@ref)). The BDS-3 B-CNAV messages instead
# broadcast BDGIM, reached through `bdgim_params` below — those containers
# therefore report no Klobuchar set.
function klobuchar_params(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
)
    d = decoder.data
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return BeiDouKlobucharParams(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end

"""
    ntcm_g_params(decoder) -> Union{NTCMGParams,Nothing}

NTCM-G Effective Ionisation Level coefficients decoded from a Galileo navigation
message (I/NAV on E1-B or E5b, or F/NAV on E5a — all broadcast `a_i0…a_i2`), or
`nothing` if they (or the week number) have not been decoded yet or the decoder
carries no such set. Galileo E6-B is in the latter group: C/NAV broadcasts no
ionospheric coefficients, so it falls to the generic method rather than being
dispatched here on a field it does not have.
"""
ntcm_g_params(decoder) = nothing
function ntcm_g_params(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.AbstractGalileoEphemerisData},
)
    d = decoder.data
    (isnothing(d.a_i0) || isnothing(d.a_i1) || isnothing(d.a_i2) || isnothing(d.WN)) &&
        return nothing
    return NTCMGParams(d.a_i0, d.a_i1, d.a_i2, d.WN)
end

"""
    bdgim_params(decoder) -> Union{BDGIMParams,Nothing}

The nine BDGIM coefficients decoded from a BDS-3 B-CNAV message — B-CNAV1 on B1C,
B-CNAV2 on B2a, B-CNAV3 on B2b, all three broadcasting the same set — or
`nothing` if they, or the week number the model needs to
date its epoch, have not been decoded yet, or the decoder carries no such set. The
legacy D1/D2 message (`BeiDouDNAVData`, on B1I/B3I) is in the latter group: it
broadcasts Klobuchar instead, so it falls to this generic method rather than being
dispatched on fields it does not have.
"""
bdgim_params(decoder) = nothing

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
function bdgim_params(decoder::GNSSDecoder.GNSSDecoderState{<:AbstractBeiDouCNAVData})
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
    return BDGIMParams(
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

"""
    select_ionospheric_correction(states)
        -> Union{KlobucharParams,BeiDouKlobucharParams,NTCMGParams,BDGIMParams,Nothing}

Scan all (healthy) satellite decoders and pick the single ionospheric correction
to apply to the whole solve, in a fixed order of preference:

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
global to a constellation, so the first decoder that carries each set is used.
"""
function select_ionospheric_correction(states)
    gps_klobuchar = nothing
    beidou_klobuchar = nothing
    bdgim = nothing
    ntcm_g = nothing
    for state in states
        if state.decoder.data isa GNSSDecoder.AbstractBeiDouData
            # The two BeiDou message families broadcast different models: the legacy
            # D1/D2 one Klobuchar, the BDS-3 B-CNAV ones BDGIM. Both accessors are
            # asked; each answers `nothing` for the family that is not its own.
            beidou_klobuchar === nothing &&
                (beidou_klobuchar = klobuchar_params(state.decoder))
            bdgim === nothing && (bdgim = bdgim_params(state.decoder))
        else
            gps_klobuchar === nothing && (gps_klobuchar = klobuchar_params(state.decoder))
        end
        ntcm_g === nothing && (ntcm_g = ntcm_g_params(state.decoder))
    end
    # The order of preference documented above, top rung first.
    ntcm_g !== nothing && return ntcm_g
    bdgim !== nothing && return bdgim
    gps_klobuchar !== nothing && return gps_klobuchar
    return beidou_klobuchar          # BeiDouKlobucharParams, or nothing if neither
end

"""
    ionospheric_delay(correction, system, elevation, azimuth, lla, time_of_week) -> Float64

Slant ionospheric group delay in metres for one satellite (`system` is the
satellite's GNSS, used for its carrier frequency), using the constellation-wide
`correction` returned by [`select_ionospheric_correction`](@ref):

- `::Nothing` → `0.0` (no coefficients were decoded).
- [`KlobucharParams`](@ref) → Klobuchar model (IS-GPS-200, Fig. 20-4).
- [`BeiDouKlobucharParams`](@ref) → BeiDou's Klobuchar variant (BDS-SIS-ICD-B1I-3.0 §5.2.4.7).
- [`NTCMGParams`](@ref) → NTCM-G model.
- [`BDGIMParams`](@ref) → BDGIM (BDS-SIS-ICD-B1C-1.0 §7.8.2).

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
[`_elevation_azimuth`](@ref)`(ENUfromECEF(user, wgs84), sat)`.
"""
ionospheric_delay(::Nothing, system, elevation, azimuth, lla, time_of_week) = 0.0

function ionospheric_delay(p::KlobucharParams, system, elevation, azimuth, lla, time_of_week)
    # IS-GPS-200 works in semicircles: lat/lon in deg/180, elevation/azimuth in rad/π.
    l1_seconds = klobuchar_group_delay(
        lla.lat / 180,
        lla.lon / 180,
        elevation / π,
        azimuth / π,
        time_of_week,
        (p.α_0, p.α_1, p.α_2, p.α_3),
        (p.β_0, p.β_1, p.β_2, p.β_3),
    )
    # The Klobuchar broadcast coefficients define the group delay at the GPS L1
    # frequency (IS-GPS-200). The ionospheric delay scales as 1/f², so rescale it
    # to this satellite's actual carrier frequency.
    f = get_center_frequency(system)
    return SPEED_OF_LIGHT * l1_seconds * (get_center_frequency(GPSL1CA) / f)^2
end

function ionospheric_delay(
    p::BeiDouKlobucharParams,
    system,
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
        (p.α_0, p.α_1, p.α_2, p.α_3),
        (p.β_0, p.β_1, p.β_2, p.β_3),
    )
    # The coefficients define the group delay along the B1I propagation path (the
    # ICD's I_B1I) — 1561.098 MHz, not L1. The delay scales as 1/f², so rescale
    # from B1I to this satellite's actual carrier.
    f = get_center_frequency(system)
    return SPEED_OF_LIGHT * b1i_seconds *
           (get_center_frequency(GNSSSignals.BeiDouB1I) / f)^2
end

function ionospheric_delay(p::NTCMGParams, system, elevation, azimuth, lla, time_of_week)
    doy, ut = _galileo_doy_and_ut(p.week_number, time_of_week)
    stec = ntcm_g_stec(elevation, azimuth, lla, doy, ut, p.a_i0, p.a_i1, p.a_i2) # TECU
    f = ustrip(Hz, get_center_frequency(system))
    # Eq. 1: group delay [m] = 40.3 / f² · STEC, with STEC in electrons/m² (1 TECU = 1e16).
    return 40.3 / f^2 * stec * 1.0e16
end

function ionospheric_delay(p::BDGIMParams, system, elevation, azimuth, lla, time_of_week)
    mjd = _bdgim_modified_julian_date(p.week_number, time_of_week)
    α = (p.α_1, p.α_2, p.α_3, p.α_4, p.α_5, p.α_6, p.α_7, p.α_8, p.α_9)
    stec = bdgim_stec(elevation, azimuth, lla, mjd, α)   # TECU
    f = ustrip(Hz, get_center_frequency(system))
    # Eq. 7-6: T_ion [m] = M_F · 40.28e16/f² · VTEC, i.e. 40.28e16/f² · STEC with STEC
    # in TECU. `f` is this satellite's own carrier, exactly as the ICD says ("the
    # carrier frequency of the current signal"); there is no reference frequency to
    # rescale from, unlike the Klobuchar branch. Note the constant is BDGIM's own
    # 40.28, not NTCM-G's 40.3 — a 0.05 % difference, kept ICD-exact per model.
    return 40.28 / f^2 * stec * 1.0e16
end

"""
    klobuchar_group_delay(φ_u, λ_u, E, A, gps_time, α, β)

Klobuchar single-frequency ionospheric group delay for GPS L1 (IS-GPS-200N,
Fig. 20-4), returned in **seconds**. All angles are in **semicircles**:
`φ_u`/`λ_u` are the user geodetic latitude/longitude, `E`/`A` the satellite
elevation/azimuth. `gps_time` is GPS system time in seconds. `α`/`β` are the
4-element Klobuchar coefficient tuples (SI units).
"""
function klobuchar_group_delay(φ_u, λ_u, E, A, gps_time, α, β)
    # Earth-centred angle between user and ionospheric pierce point (semicircles)
    ψ = 0.0137 / (E + 0.11) - 0.022
    # Geodetic latitude of the ionospheric pierce point (IPP), clamped per ICD
    φ_i = clamp(φ_u + ψ * cos(A * π), -0.416, 0.416)
    # Geodetic longitude of the IPP
    λ_i = λ_u + ψ * sin(A * π) / cos(φ_i * π)
    # Geomagnetic latitude of the IPP
    φ_m = φ_i + 0.064 * cos((λ_i - 1.617) * π)
    # Local time at the IPP (seconds), wrapped to [0, 86400)
    t = mod(4.32e4 * λ_i + gps_time, 86400.0)
    # Obliquity / slant factor
    F = 1.0 + 16.0 * (0.53 - E)^3
    # Amplitude (s) and period (s) of the cosine model, with ICD floors
    AMP = max(α[1] + φ_m * (α[2] + φ_m * (α[3] + φ_m * α[4])), 0.0)
    PER = max(β[1] + φ_m * (β[2] + φ_m * (β[3] + φ_m * β[4])), 72000.0)
    x = 2π * (t - 50400.0) / PER
    return abs(x) < 1.57 ? F * (5.0e-9 + AMP * (1 - x^2 / 2 + x^4 / 24)) : F * 5.0e-9
end

"""
    beidou_klobuchar_group_delay(φ_u, λ_u, E, A, t_E, α, β)

BeiDou's Klobuchar-variant ionospheric group delay along the B1I propagation
path (BDS-SIS-ICD-B1I-3.0 §5.2.4.7), returned in **seconds**. All angles are in
**radians**, the ICD's working unit: `φ_u`/`λ_u` are the user geodetic
latitude/longitude, `E`/`A` the satellite elevation/azimuth. `t_E` is the BDT
second of week. `α`/`β` are the 4-element coefficient tuples (SI units; their
`s·π⁻ⁿ` scale factors make the polynomial argument `φ_M/π` — semicircles,
exactly like the GPS set).

The same eight-coefficient half-cosine idea as [`klobuchar_group_delay`](@ref),
but not the same algorithm — the ICD replaces IS-GPS-200's fitted
approximations with exact spherical geometry and drops the geomagnetic frame:

- The Earth-central angle to the ionospheric pierce point is the exact
  `π/2 − E − arcsin(R/(R+h)·cos E)` at R = 6378 km, h = 375 km, and the
  obliquity the matching exact `1/√(1 − (R/(R+h)·cos E)²)`, instead of the GPS
  model's `0.0137/(E+0.11) − 0.022` and `1 + 16(0.53−E)³` fits.
- The latitude feeding the A₂/A₄ polynomials is the pierce point's *geographic*
  latitude from the spherical triangle — no `+0.064·cos(λ−1.617)` geomagnetic
  conversion, no ±0.416-semicircle clamp — and its longitude uses `arcsin`
  rather than a flat-Earth division.
- The period A₄ is clamped from above at 172800 s as well as from below at
  72000 s, and the cosine is evaluated exactly rather than by the two-term
  Taylor expansion inside |x| < 1.57.

The 5·10⁻⁹ s night floor and the |t − 50400| < A₄/4 day-time window are those
of the GPS model.
"""
function beidou_klobuchar_group_delay(φ_u, λ_u, E, A, t_E, α, β)
    # Sine of the pierce point's geocentric zenith angle: R/(R+h)·cosE, with the
    # ICD's constants R = 6378 km (mean Earth radius) and h = 375 km (ionosphere
    # height). Both the pierce-point geometry and the obliquity flow from it.
    sinz = 6378.0 / (6378.0 + 375.0) * cos(E)
    # Earth-central angle between user and pierce point (radians), exact
    ψ = π / 2 - E - asin(sinz)
    # Geographic — not geomagnetic — latitude and longitude of the pierce point
    φ_M = asin(sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(A))
    λ_M = λ_u + asin(sin(ψ) * sin(A) / cos(φ_M))
    # Local time at the pierce point (seconds), wrapped to [0, 86400)
    t = mod(t_E + λ_M * 43200.0 / π, 86400.0)
    # Amplitude (s) and period (s) of the cosine model, with the ICD's clamps —
    # the period is bounded on both sides, unlike the GPS model's floor
    φ = φ_M / π
    A_2 = max(α[1] + φ * (α[2] + φ * (α[3] + φ * α[4])), 0.0)
    A_4 = clamp(β[1] + φ * (β[2] + φ * (β[3] + φ * β[4])), 72000.0, 172800.0)
    # Vertical delay I′z (s): a true cosine within a quarter period of 14:00 local
    # time, the night floor outside it
    x = 2π * (t - 50400.0) / A_4
    I_z = abs(t - 50400.0) < A_4 / 4 ? 5.0e-9 + A_2 * cos(x) : 5.0e-9
    # Slant delay along the B1I path: the exact obliquity, from the same geometry
    return I_z / sqrt(1.0 - sinz^2)
end

"""
    _elevation_azimuth(enu_from_ecef::ENUfromECEF, sat_position) -> (elevation, azimuth)

Elevation and azimuth (radians) of `sat_position` (ECEF) in the local East-North-Up
frame defined by `enu_from_ecef = ENUfromECEF(user_position, wgs84)`. Azimuth is
measured clockwise from North. The transform is taken precomputed so it can be built
once per user position and reused across satellites.
"""
function _elevation_azimuth(enu_from_ecef::ENUfromECEF, sat_position)
    sat_enu = get_sat_enu(enu_from_ecef, ECEF(sat_position))
    elevation = sat_enu.ϕ
    # `SphericalFromCartesian` measures θ counter-clockwise from East (the ENU +x
    # axis); the ionospheric models use azimuth measured clockwise from North, i.e.
    # π/2 − θ. Only cos/sin of the azimuth are used downstream, so the wrap is moot.
    azimuth = π / 2 - sat_enu.θ
    return elevation, azimuth
end

# ===========================================================================
#  NTCM-G — Galileo single-frequency ionospheric model
#  (European GNSS (Galileo) NTCM-G Ionospheric Model Description, Issue 1.0,
#   May 2022). Driven by the broadcast Effective Ionisation Level coefficients
#   a_i0, a_i1, a_i2. Returns slant TEC (TECU); equation numbers below refer to
#   that document.
# ===========================================================================

# NTCM-G model coefficients k1..k12 (Table 3)
const _NTCM_K = (
    0.92519,
    0.16951,
    0.00443,
    0.06626,
    0.00899,
    0.21289,
    -0.15414,
    -0.38439,
    1.14023,
    1.20556,
    1.41808,
    0.13985,
)
const _NTCM_RE = 6371.0                 # Earth mean radius [km] (Table 2)
const _NTCM_HI = 450.0                  # ionospheric pierce point height [km]
const _NTCM_GNP_LAT = deg2rad(79.74)    # geomagnetic North pole latitude
const _NTCM_GNP_LON = deg2rad(-71.78)   # geomagnetic North pole longitude

# Effective Ionisation Level Azpar [sfu] from the broadcast coefficients (Eq. 2).
function _azpar(a_i0, a_i1, a_i2)
    radicand = a_i0^2 + 1633.33 * a_i1^2 + 4802000.0 * a_i2^2 + 3266.67 * a_i0 * a_i2
    return sqrt(max(radicand, 0.0))
end

# Ionospheric pierce point geographic latitude/longitude [rad] (Eq. 24-26),
# given user geodetic lat/lon and satellite elevation/azimuth [rad].
function _pierce_point(φ_u, λ_u, elevation, azimuth)
    ψ = π / 2 - elevation - asin(_NTCM_RE / (_NTCM_RE + _NTCM_HI) * cos(elevation))
    # Clamp to asin's domain: the φ_pp argument is a unit dot-product that can
    # overshoot ±1 by a rounding ulp, and the λ_pp argument genuinely diverges as
    # the pierce point approaches a pole (cos(φ_pp) → 0).
    φ_pp = asin(clamp(sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth), -1.0, 1.0))
    λ_pp = λ_u + asin(clamp(sin(ψ) * sin(azimuth) / cos(φ_pp), -1.0, 1.0))
    return φ_pp, λ_pp
end

# Sun's declination [rad] for the day of year (Eq. 28).
_sun_declination(doy) = deg2rad(23.44) * sin(deg2rad(0.9856 * (doy - 80.7)))

# Modified Single Layer Model mapping function (Eq. 32-33).
function _mslm_mapping_function(elevation)
    sinz = _NTCM_RE / (_NTCM_RE + _NTCM_HI) * sin(0.9782 * (π / 2 - elevation))
    return 1.0 / sqrt(1.0 - sinz^2)
end

"""
    ntcm_g_vtec(φ_pp, λ_pp, doy, ut, azpar) -> Float64

Vertical TEC in TECU at the ionospheric pierce point (geographic latitude/longitude
`φ_pp`/`λ_pp` in radians) for day of year `doy`, universal time `ut` in hours, and
Effective Ionisation Level `azpar` in solar flux units. Implements the NTCM-G model
`VTEC = F1·F2·F3·F4·F5` (Eq. 3-15).
"""
function ntcm_g_vtec(φ_pp, λ_pp, doy, ut, azpar)
    k = _NTCM_K
    lt = ut + rad2deg(λ_pp) / 15                      # local time [h] (Eq. 27)
    δ = _sun_declination(doy)

    # Solar zenith angle dependence (Eq. 29-30)
    cosχ3 = cos(φ_pp - δ) + 0.4                        # cosχ*** (PF1 = 0.4)
    cosχ2 = cos(φ_pp - δ) - (2 / π) * φ_pp * sin(δ)    # cosχ**

    # F1 — local-time dependency (Eq. 4-7)
    V_D = 2π * (lt - 14) / 24
    V_SD = 2π * lt / 12
    V_TD = 2π * lt / 8
    F1 =
        cosχ3 +
        cosχ2 * (
            k[1] * cos(V_D) +
            k[2] * cos(V_SD) +
            k[3] * sin(V_SD) +
            k[4] * cos(V_TD) +
            k[5] * sin(V_TD)
        )

    # F2 — seasonal dependency (Eq. 8-10)
    V_A = 2π * (doy - 18) / 365.25
    V_SA = 4π * (doy - 6) / 365.25
    F2 = 1 + k[6] * cos(V_A) + k[7] * cos(V_SA)

    # Geomagnetic latitude of the pierce point (Eq. 31)
    φ_m = asin(
        sin(φ_pp) * sin(_NTCM_GNP_LAT) +
        cos(φ_pp) * cos(_NTCM_GNP_LAT) * cos(λ_pp - _NTCM_GNP_LON),
    )

    # F3 — geomagnetic field dependency (Eq. 11), φ_m in radians
    F3 = 1 + k[8] * cos(φ_m)

    # F4 — equatorial (Appleton) anomaly dependency (Eq. 12-14), φ_m in degrees
    φ_m_deg = rad2deg(φ_m)
    EC1 = -(φ_m_deg - 16.0)^2 / (2 * 12.0^2)
    EC2 = -(φ_m_deg + 10.0)^2 / (2 * 13.0^2)
    F4 = 1 + k[9] * exp(EC1) + k[10] * exp(EC2)

    # F5 — solar activity dependency (Eq. 15)
    F5 = k[11] + k[12] * azpar

    return max(F1 * F2 * F3 * F4 * F5, 0.0)
end

"""
    ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2) -> Float64

Slant TEC in TECU along the user→satellite line of sight using NTCM-G, for the
satellite `elevation`/`azimuth` (radians) seen from the user geodetic position
`lla` (a `Geodesy.LLA`), day of year `doy`, universal time `ut` (hours), and the
broadcast Galileo Effective Ionisation Level coefficients `a_i0`/`a_i1`/`a_i2`.
"""
function ntcm_g_stec(elevation, azimuth, lla, doy, ut, a_i0, a_i1, a_i2)
    φ_pp, λ_pp = _pierce_point(deg2rad(lla.lat), deg2rad(lla.lon), elevation, azimuth)
    vtec = ntcm_g_vtec(φ_pp, λ_pp, doy, ut, _azpar(a_i0, a_i1, a_i2))
    return _mslm_mapping_function(elevation) * vtec
end

# Day of year and universal time (hours) from a Galileo System Time (GST) week
# number and time of week [s]. The GST epoch is taken from GNSSSignals
# (`get_system_start_time(GST())` = 1999-08-21T23:59:47 UTC); GST is continuous
# (offset from UTC by leap seconds, ~18 s), negligible for the day-of-year / UT
# inputs of NTCM-G.
function _galileo_doy_and_ut(week_number, time_of_week)
    epoch = get_system_start_time(GST())
    t = epoch +
        Millisecond(round(Int, (week_number * SECONDS_PER_WEEK + time_of_week) * 1000))
    return dayofyear(t), hour(t) + minute(t) / 60 + second(t) / 3600 + millisecond(t) / 3.6e6
end

# ===========================================================================
#  BDGIM — BeiDou Global Ionospheric delay correction Model
#  (BDS-SIS-ICD-B1C-1.0 §7.8, Eq. 7-6 … 7-17 and Tables 7-10/7-11/7-12; the
#   same model appears verbatim as BDS-SIS-ICD-B2a-1.0 §7.8 Eq. 7-6 … 7-17 and
#   BDS-SIS-ICD-B2b-1.0 §7.7 Eq. 7-5 … 7-16, where the equation and table
#   numbers are each one lower.) Equation numbers below are the B1C ones.
#
#  A modified spherical-harmonics VTEC model: nine broadcast coefficients α_1…α_9
#  weight nine normalised Legendre harmonics of the ionospheric pierce point's
#  *solar-fixed geomagnetic* coordinates, on top of a purely predictive part A_0
#  built from 17 further harmonics whose amplitudes are 12-period Fourier series
#  with coefficients fixed in the ICD (not broadcast). VTEC is mapped to the line
#  of sight by a single-layer obliquity factor. Returns TECU; the metre conversion
#  lives in `ionospheric_delay`.
# ===========================================================================

# Suggested parameter values, BDS-SIS-ICD-B1C-1.0 §7.8.2 (the paragraph after
# Eq. 7-17). Re and H_ion in km — only their ratio enters, so the unit is free.
const _BDGIM_RE = 6378.0                    # mean radius of the Earth [km]
const _BDGIM_HION = 400.0                   # ionospheric single-layer shell altitude [km]
const _BDGIM_POLE_LAT = 80.27 / 180 * π     # north magnetic pole geographic latitude [rad]
const _BDGIM_POLE_LON = -72.58 / 180 * π    # north magnetic pole geographic longitude [rad]

# Degree/order (n_i, m_i) of the nine broadcast harmonics, Table 7-11. A negative
# m selects the sine rather than the cosine term (Eq. 7-11).
const _BDGIM_NM =
    ((0, 0), (1, 0), (1, 1), (1, -1), (2, 0), (2, 1), (2, -1), (2, 2), (2, -2))

# Degree/order (n_j, m_j) of the seventeen harmonics of the predictive part,
# read from the header row of Table 7-12. (Eq. 7-14 says "Table 7-11", which
# lists only the nine broadcast pairs — an ICD cross-reference slip; the j-indexed
# pairs exist only in Table 7-12's header.)
const _BDGIM_NM_PREDICTED = (
    (3, 0), (3, 1), (3, -1), (3, 2), (3, -2), (3, 3), (3, -3),   # degree 3
    (4, 0), (4, 1), (4, -1), (4, 2), (4, -2),                    # degree 4
    (5, 0), (5, 1), (5, -1), (5, 2), (5, -2),                    # degree 5
)

# Prediction periods T_k [days], k = 1…12, from the rightmost column of Table 7-12.
const _BDGIM_PERIODS =
    (1.0, 0.5, 0.33, 14.6, 27.0, 121.6, 182.51, 365.25, 4028.71, 2014.35, 1342.90, 1007.18)

# ---------------------------------------------------------------------------
# Table 7-12, the non-broadcast coefficients [TECu]: the constant row a_{0,j}
# and the twelve (a_{k,j}, b_{k,j}) pairs, each 17 entries wide (j = 1…17).
#
# TRANSCRIPTION, DOUBLE-CHECKED. Text extraction mangles the ICD's tables, so
# this table was taken three ways and the three agreed entry for entry: the text
# layer of BDS-SIS-ICD-B1C-1.0 Table 7-12, the *rendered page image* of that same
# table (the ICD PDF, printed page 50, read as four overlapping crops so every
# digit was legible), and the text layer of the independently typeset
# BDS-SIS-ICD-B2b-1.0 (2020) Table 7-10, which restates the same table.
# `test/ionosphere.jl` pins row, column and absolute-sum checksums so that a typo
# introduced by a later edit fails a test rather than quietly shifting TEC.
# ---------------------------------------------------------------------------
const _BDGIM_A0 = (
    -0.61, -1.31, -2.0, -0.03, 0.15, -0.48, -0.4, 2.28, -0.16, -0.21, -0.1, -0.13, 0.21,
    0.68, 1.06, 0.0, -0.12,
)

const _BDGIM_A = (
    (-0.51, -0.43, 0.34, -0.01, 0.17, 0.02, -0.06, 0.3, 0.44, -0.28, -0.31, -0.17, 0.04,
        0.39, -0.12, 0.12, 0.0),
    (-0.06, -0.05, 0.06, 0.17, 0.15, 0.0, 0.11, -0.05, -0.16, 0.02, 0.11, 0.04, 0.12, 0.07,
        0.02, -0.14, -0.14),
    (0.01, -0.03, 0.01, -0.01, 0.05, -0.03, 0.05, -0.03, -0.01, 0.0, -0.08, -0.04, 0.0,
        -0.02, -0.03, 0.0, -0.03),
    (-0.01, 0.0, 0.01, 0.0, 0.01, 0.0, -0.01, -0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (0.0, 0.0, 0.03, 0.01, 0.02, 0.01, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (-0.19, -0.02, 0.12, -0.1, 0.06, 0.0, -0.02, -0.08, -0.02, -0.07, 0.01, 0.03, 0.15,
        0.06, -0.05, -0.03, -0.1),
    (-0.18, 0.06, -0.55, -0.02, 0.09, -0.08, 0.0, 0.86, -0.18, -0.05, -0.07, 0.04, 0.14,
        -0.03, 0.37, -0.11, -0.12),
    (1.09, -0.14, -0.21, 0.52, 0.27, 0.0, 0.11, 0.17, 0.23, 0.35, -0.05, 0.02, -0.6, 0.02,
        0.01, 0.27, 0.32),
    (-0.34, -0.09, -1.22, 0.05, 0.15, -0.29, -0.17, 1.58, -0.06, -0.15, 0.0, 0.13, 0.28,
        -0.08, 0.62, -0.01, -0.04),
    (-0.13, 0.07, -0.37, 0.05, 0.06, -0.11, -0.07, 0.46, 0.0, -0.04, 0.01, 0.07, 0.09,
        -0.05, 0.15, -0.01, 0.01),
    (-0.06, 0.13, -0.07, 0.03, 0.02, -0.05, -0.05, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
    (-0.03, 0.08, -0.01, 0.04, 0.01, -0.02, -0.02, -0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
)

const _BDGIM_B = (
    (0.23, -0.2, -0.31, 0.16, -0.03, 0.02, 0.04, 0.18, 0.34, 0.45, 0.19, -0.25, -0.12, 0.18,
        0.4, -0.09, 0.21),
    (0.02, -0.08, -0.06, -0.11, 0.15, -0.14, 0.01, 0.01, 0.04, -0.14, -0.05, 0.08, 0.08,
        -0.01, 0.01, 0.11, -0.12),
    (0.0, -0.02, -0.03, -0.05, -0.01, -0.07, -0.03, -0.01, 0.02, -0.01, 0.03, -0.1, 0.01,
        0.05, -0.01, 0.04, 0.0),
    (0.0, -0.02, 0.01, 0.0, -0.01, 0.01, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0),
    (0.01, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (-0.09, 0.07, 0.03, 0.06, 0.09, 0.01, 0.02, 0.0, -0.04, -0.02, -0.01, 0.01, -0.1, 0.0,
        -0.01, 0.02, 0.05),
    (0.15, -0.31, 0.13, 0.05, -0.09, -0.03, 0.06, -0.36, 0.08, 0.05, 0.06, -0.02, -0.05,
        0.06, -0.2, 0.04, 0.07),
    (0.5, -0.08, -0.38, 0.36, 0.14, 0.04, 0.0, 0.25, 0.17, 0.27, -0.03, -0.03, -0.32, -0.1,
        0.2, 0.1, 0.3),
    (0.0, -0.11, -0.22, 0.01, 0.02, -0.03, -0.01, 0.49, -0.03, -0.02, 0.01, 0.02, 0.04,
        -0.04, 0.16, -0.02, -0.01),
    (0.05, 0.03, 0.07, 0.02, -0.01, 0.03, 0.02, -0.04, -0.01, -0.01, 0.02, 0.03, 0.02,
        -0.04, -0.04, -0.01, 0.0),
    (0.03, -0.02, 0.04, -0.01, -0.03, 0.02, 0.01, 0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
    (0.04, -0.02, -0.04, 0.0, -0.01, 0.0, 0.01, 0.07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0),
)

"""
    _bdgim_legendre(n, m, x) -> Float64

Classic, un-normalised associated Legendre function `P_{n,m}(x)` (Eq. 7-13),
evaluated by the ICD's own recursion rather than by a closed form, so the code
reads against the document: the sectoral seed `P_{n,n} = (2n−1)!!·(1−x²)^(n/2)`,
the one-step case `P_{m+1,m} = x·(2m+1)·P_{m,m}`, and the degree recursion
`P_{n,m} = [(2n−1)·x·P_{n−1,m} − (n+m−1)·P_{n−2,m}] / (n−m)`. `n ≥ m ≥ 0`; the
model needs no degree above 5.
"""
function _bdgim_legendre(n::Int, m::Int, x)
    # (2m−1)!! = (2m−1)(2m−3)…1, an empty product (= 1) at m = 0, so P_{0,0} = 1.
    double_factorial = 1.0
    for i in 1:m
        double_factorial *= 2i - 1
    end
    p_mm = double_factorial * (1 - x^2)^(m / 2)
    n == m && return p_mm
    p_prev = p_mm                       # P_{m,m}
    p_cur = x * (2m + 1) * p_mm         # P_{m+1,m}
    for k in (m+2):n
        p_prev, p_cur = p_cur, ((2k - 1) * x * p_cur - (k + m - 1) * p_prev) / (k - m)
    end
    return p_cur
end

# Normalisation N_{n,m} of Eq. 7-12: √[(n−m)!·(2n+1)·(2 − δ_{0,m}) / (n+m)!], with
# δ_{0,m} = 1 for m = 0 and 0 otherwise. Max degree 5 keeps (n+m)! ≤ 10! exact in Int.
_bdgim_normalization(n::Int, m::Int) =
    sqrt(factorial(n - m) * (2n + 1) * (m == 0 ? 1 : 2) / factorial(n + m))

# One normalised spherical-harmonic term of Eq. 7-11 / Eq. 7-14, at the solar-fixed
# geomagnetic latitude/longitude φ′/λ′ [rad]: P̃_{|n|,|m|}(sin φ′) times cos(m·λ′)
# for m ≥ 0 and sin(−m·λ′) for m < 0.
function _bdgim_harmonic(n::Int, m::Int, sinφ′, λ′)
    p = _bdgim_normalization(abs(n), abs(m)) * _bdgim_legendre(abs(n), abs(m), sinφ′)
    return m ≥ 0 ? p * cos(m * λ′) : p * sin(-m * λ′)
end

"""
    _bdgim_pierce_point(φ_u, λ_u, elevation, azimuth) -> (φ_g, λ_g)

Geographic latitude/longitude [rad] of the Earth projection of the ionospheric
pierce point (Eq. 7-7 and 7-8), for a user at geodetic `φ_u`/`λ_u` [rad] seeing a
satellite at `elevation`/`azimuth` [rad]. The Earth-central angle is
`ψ = π/2 − E − arcsin(Re/(Re+H_ion)·cos E)` at the ICD's shell height.
"""
function _bdgim_pierce_point(φ_u, λ_u, elevation, azimuth)
    ψ = π / 2 - elevation - asin(_BDGIM_RE / (_BDGIM_RE + _BDGIM_HION) * cos(elevation))
    # Clamp only against a rounding ulp past ±1; the argument is a unit dot product.
    φ_g = asin(
        clamp(
            sin(φ_u) * cos(ψ) + cos(φ_u) * sin(ψ) * cos(azimuth),
            -1.0,
            1.0,
        ),
    )
    # Eq. 7-8 writes a single-argument arctan of a quotient. Taken literally that
    # folds the longitude offset into (−π/2, π/2), which is wrong whenever the
    # denominator is negative — the pierce point is then on the far side of the
    # pole from the user. The printed numerator and denominator are exactly the
    # textbook great-circle destination pair (cos φ_g · sin Δλ, cos φ_g · cos Δλ)
    # up to the positive factor cos φ_g, so the two-argument `atan(y, x)` of the
    # pair as printed recovers the full-circle longitude the model needs.
    λ_g = λ_u + atan(
        sin(ψ) * sin(azimuth) * cos(φ_u),
        cos(ψ) - sin(φ_u) * sin(φ_g),
    )
    return φ_g, λ_g
end

"""
    _bdgim_solar_geomagnetic(φ_g, λ_g, mjd) -> (φ′, λ′)

Latitude/longitude [rad] of the pierce point in the solar-fixed geomagnetic frame
(Eq. 7-9 and 7-10), from its geographic `φ_g`/`λ_g` [rad] and the Modified Julian
Date `mjd` of the epoch. The pole is the ICD's fixed one (80.27°N, 72.58°W), and
the sun's mean geographic longitude is `S_lon = π·(1 − 2·(t − ⌊t⌋))` — the
sub-solar meridian sweeping west once per MJD day, from 180°E at MJD midnight
through Greenwich at 12:00 UT (⌊t⌋ being MJD midnight, and `int()` in the ICD).

Only the longitude is referred to the sun; `φ′ = φ_m` (Eq. 7-10).
"""
function _bdgim_solar_geomagnetic(φ_g, λ_g, mjd)
    Δλ = λ_g - _BDGIM_POLE_LON
    φ_m = asin(
        clamp(
            sin(_BDGIM_POLE_LAT) * sin(φ_g) +
            cos(_BDGIM_POLE_LAT) * cos(φ_g) * cos(Δλ),
            -1.0,
            1.0,
        ),
    )
    # As in Eq. 7-8, the quotient's two-argument arctan is the intended one: the
    # geomagnetic longitude runs over the full circle, and the printed pair is the
    # textbook (sin λ_m, cos λ_m) pair up to a positive scale, so `atan(y, x)` of
    # it as printed is exact.
    λ_m = atan(
        cos(φ_g) * sin(Δλ) * cos(_BDGIM_POLE_LAT),
        sin(_BDGIM_POLE_LAT) * sin(φ_m) - sin(φ_g),
    )
    s_lon = π * (1 - 2 * (mjd - floor(mjd)))
    # The subtracted term is the same Eq. 7-9 longitude evaluated at the sun's own
    # position (φ_g = 0, λ_g = S_lon), written out.
    λ′ = λ_m - atan(
        sin(s_lon - _BDGIM_POLE_LON),
        sin(_BDGIM_POLE_LAT) * cos(s_lon - _BDGIM_POLE_LON),
    )
    return φ_m, λ′
end

"""
    _bdgim_reference_epoch(mjd) -> Float64

The prediction reference epoch `t_p` (Eq. 7-15) for an epoch at Modified Julian
Date `mjd`: the odd hour (01:00, 03:00, …, 23:00) of that day nearest the epoch,
returned in days (MJD). The odd hours are two hours apart and span the whole day,
so the nearest one is always within an hour and always on the same day — no wrap
is needed. `t_p` therefore steps in two-hour blocks, which quantises the
prediction *amplitudes* β_j; VTEC itself still varies continuously in between,
because the solar-fixed longitude λ′ tracks the true epoch rather than `t_p`.
"""
function _bdgim_reference_epoch(mjd)
    day = floor(mjd)
    hour_of_day = (mjd - day) * 24
    return day + (2 * floor(hour_of_day / 2) + 1) / 24
end

"""
    _bdgim_prediction_amplitudes(t_p) -> NTuple{17,Float64}

The seventeen prediction amplitudes `β_j` (Eq. 7-15) at reference epoch `t_p`
(days, MJD): `β_j = a_{0,j} + Σ_{k=1}^{12} [a_{k,j}·cos(ω_k·t_p) + b_{k,j}·sin(ω_k·t_p)]`
with `ω_k = 2π/T_k` and the non-broadcast coefficients of Table 7-12. `t_p` is an
absolute MJD, so the series is a continuous function of it — no day-of-year or
epoch wrap enters, and the ~4028.71-day term is what carries the solar cycle.
"""
function _bdgim_prediction_amplitudes(t_p)
    # The twelve (sin, cos) pairs do not depend on j, so they are formed once and
    # shared by all seventeen amplitudes rather than recomputed 17 times.
    trig = ntuple(k -> sincos(2π / _BDGIM_PERIODS[k] * t_p), Val(12))
    return ntuple(Val(17)) do j
        β = _BDGIM_A0[j]
        for k in eachindex(trig)
            s, c = trig[k]
            β += _BDGIM_A[k][j] * c + _BDGIM_B[k][j] * s
        end
        β
    end
end

"""
    bdgim_vtec(φ_g, λ_g, mjd, α) -> Float64

Vertical TEC in TECU at the ionospheric pierce point whose geographic
latitude/longitude are `φ_g`/`λ_g` (radians), at Modified Julian Date `mjd`, from
the nine broadcast BDGIM coefficients `α` (a 9-tuple, TECu). Implements
`VTEC = A_0 + Σ_{i=1}^{9} α_i·A_i` (Eq. 7-16) with the broadcast harmonics `A_i`
of Eq. 7-11 and the predictive part `A_0` of Eq. 7-14/7-15.

The result is floored at zero. That floor is *not* in the ICD, which prescribes no
handling of a non-positive VTEC; but a negative vertical electron content is
unphysical and would enter [`ionospheric_delay`](@ref) as a group *advance*, so it
is clamped exactly as the NTCM-G branch above clamps its own product. It is
reachable rather than theoretical: A_0 spans degrees 3-5 only and so has no level
of its own — it is a ±10-15 TECu structure term, negative over much of the globe —
and a broadcast set too small to lift it (an all-zero α is a legal message) leaves
the sum below zero.
"""
function bdgim_vtec(φ_g, λ_g, mjd, α)
    φ′, λ′ = _bdgim_solar_geomagnetic(φ_g, λ_g, mjd)
    sinφ′ = sin(φ′)
    # A_0, the predictive part: seventeen harmonics weighted by the Fourier
    # amplitudes β_j at the nearest odd hour (Eq. 7-14).
    β = _bdgim_prediction_amplitudes(_bdgim_reference_epoch(mjd))
    a_0 = sum(
        β[j] * _bdgim_harmonic(_BDGIM_NM_PREDICTED[j]..., sinφ′, λ′) for
        j in eachindex(β)
    )
    # The broadcast part: nine harmonics weighted by the decoded α_i (Eq. 7-11).
    broadcast_part = sum(
        α[i] * _bdgim_harmonic(_BDGIM_NM[i]..., sinφ′, λ′) for i in eachindex(α)
    )
    return max(a_0 + broadcast_part, 0.0)
end

"""
    bdgim_mapping_function(elevation) -> Float64

BDGIM's single-layer obliquity factor `M_F` (Eq. 7-17), converting vertical to
slant TEC for a satellite at `elevation` (radians):
`1/√(1 − (Re/(Re+H_ion)·cos E)²)`. It is 1 at the zenith and grows to ~2.95 at the
horizon.
"""
bdgim_mapping_function(elevation) =
    1 / sqrt(1 - (_BDGIM_RE / (_BDGIM_RE + _BDGIM_HION) * cos(elevation))^2)

"""
    bdgim_stec(elevation, azimuth, lla, mjd, α) -> Float64

Slant TEC in TECU along the user→satellite line of sight using BDGIM, for the
satellite `elevation`/`azimuth` (radians) seen from the user geodetic position
`lla` (a `Geodesy.LLA`), at Modified Julian Date `mjd`, with the nine broadcast
coefficients `α` (a 9-tuple, TECu). This is `M_F · VTEC` — steps (1) through (5)
of BDS-SIS-ICD-B1C-1.0 §7.8.2; step (6), the metre conversion, is
[`ionospheric_delay`](@ref).
"""
function bdgim_stec(elevation, azimuth, lla, mjd, α)
    φ_g, λ_g = _bdgim_pierce_point(deg2rad(lla.lat), deg2rad(lla.lon), elevation, azimuth)
    return bdgim_mapping_function(elevation) * bdgim_vtec(φ_g, λ_g, mjd, α)
end

# Modified Julian Date of the epoch given by a BDT week number and time of week [s].
# BDGIM's only time arguments are the sun's mean longitude (a fraction of an MJD
# day, Eq. 7-10) and the prediction reference epoch t_p (Eq. 7-15), so what is
# needed is a date, not a time scale: `get_system_start_time(BDT())` is the UTC
# calendar label of BDT week 0 (2006-01-01T00:00:00), and adding BDT seconds to it
# ignores the leap seconds inserted since — ≤ 5 s of error on an argument whose
# fastest term has a two-hour granularity and a one-day period.
const _MJD_EPOCH = DateTime(1858, 11, 17)   # MJD 0, by definition
function _bdgim_modified_julian_date(week_number, time_of_week)
    epoch_mjd = Dates.value(get_system_start_time(BDT()) - _MJD_EPOCH) / 86_400_000
    return epoch_mjd + (week_number * SECONDS_PER_WEEK + time_of_week) / 86_400
end
