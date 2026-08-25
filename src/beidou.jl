# ===========================================================================
#  BeiDou (BDS) support
#
#  Everything that differs between BeiDou and the GPS/Galileo messages this
#  package already handled, gathered in one file: the ephemeris flavour, the GEO
#  reference frame, the time of week, the per-signal group delays, the BDT–GPST
#  offset, and the Klobuchar coefficients.
#
#  Notably absent: any accessor bridging field names. The clock polynomial and
#  the two reference epochs are spelled `a_f0…a_f2`, `t_0c` and `t_0e` on every
#  BeiDou container, the same as on the GPS and Galileo ones, so the generic
#  `clock_bias` / `clock_reference_time` / `ephemeris_reference_time` methods
#  cover BeiDou with nothing added here. That is a property of GNSSDecoder rather
#  than of this file — those containers used three different spellings between
#  them until the decoder settled on one (see its CONTEXT.md, "Field naming").
#
#  Five signals are covered, in two message families (GNSSDecoder ≥ 3.14):
#
#    - **B1I / B3I** — the legacy D1 (MEO/IGSO) and D2 (GEO) message, one shared
#      `BeiDouDNAVData`. Directly-broadcast Keplerian elements, exactly like GPS
#      LNAV and Galileo I/NAV, so the generic `orbital_elements` covers it.
#    - **B1C / B2a / B2b** — the BDS-3 B-CNAV1/2/3 messages. Quasi-Keplerian
#      elements like GPS CNAV, but off BeiDou's own reference values and with an
#      outright `Ω̇` rather than GPS's `ΔΩ̇`.
#
#  BDS-3 also broadcasts BDGIM ionospheric coefficients (`α_1…α_9`) on the
#  B-CNAV messages. That model is not implemented here; a BeiDou-only solve
#  falls back to the legacy D1/D2 Klobuchar set (`klobuchar_params` below) when
#  a B1I or B3I satellite is present, and to no ionospheric correction
#  otherwise. See `select_ionospheric_correction`.
# ===========================================================================

"""
    BeiDouCNAVData

The three BDS-3 B-CNAV messages — B-CNAV1 on B1C, B-CNAV2 on B2a, B-CNAV3 on
B2b. They share a quasi-Keplerian ephemeris (`A_REF + ΔA`, `Δn_0`, `Ȧ`, …), so
[`orbital_elements`](@ref) dispatches on the union rather than on each in turn.

The legacy D1/D2 message (`BeiDouDNAVData`, on B1I and B3I) is deliberately not
part of this: it broadcasts `sqrt_A` and `Δn` directly and so takes the generic
Keplerian path.
"""
const BeiDouCNAVData = Union{
    GNSSDecoder.BeiDouB1CData,
    GNSSDecoder.BeiDouB2aData,
    GNSSDecoder.BeiDouB2bData,
}

# ---- Ephemeris --------------------------------------------------------------

# Broadcast orbit-type codes (BDS-SIS-ICD-B1C-1.0 Table 7-3): 1 = GEO, 2 = IGSO,
# 3 = MEO, 0 reserved.
const BEIDOU_SAT_TYPE_GEO = 1
const BEIDOU_SAT_TYPE_IGSO = 2
const BEIDOU_SAT_TYPE_MEO = 3

"""
Reference semi-major axis for a BeiDou B-CNAV `sat_type` (BDS-SIS-ICD-B1C-1.0
Table 7-3 and its B2a/B2b counterparts): 27 906 100 m for MEO, 42 162 200 m for
IGSO and GEO. `sat_type` is 1 = GEO, 2 = IGSO, 3 = MEO, and 0 is reserved.

The reserved value throws rather than picking one. It cannot arrive through
[`calc_pvt`](@ref) — all three B-CNAV messages gate
`is_decoding_completed_for_positioning` on the orbit type being one of the three
defined ones, so such a satellite is filtered out before it is propagated — and
the two references differ by 14 256 100 m, so a satellite that has not said which
one its `ΔA` is relative to has not given its orbit at all. There is no
defensible guess to make, which is the same reason the group-delay terms the
decoder guarantees are read unwrapped: a state that cannot occur should surface
as an error, not as a plausible-looking number.
"""
function beidou_reference_semi_major_axis(sat_type)
    sat_type == BEIDOU_SAT_TYPE_MEO && return 27_906_100.0
    sat_type in (BEIDOU_SAT_TYPE_GEO, BEIDOU_SAT_TYPE_IGSO) && return 42_162_200.0
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
function orbital_elements(data::BeiDouCNAVData, μ, t_k)
    A = beidou_reference_semi_major_axis(data.sat_type) + data.ΔA
    n = sqrt(μ / A^3) + data.Δn_0 + 0.5 * data.Δn_0_dot * t_k
    (A = A, sqrt_A = sqrt(A), A_dot = data.A_dot, n = n, Ω_dot = data.Ω_dot)
end

# ---- The GEO reference frame ------------------------------------------------

"""
Inclination of the BeiDou GEO auxiliary frame relative to BDCS, −5° in radians
(BDS-SIS-ICD-B1I-3.0 §5.2.4.12).
"""
const BEIDOU_GEO_TILT = -5.0 * π / 180

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
    st, ct = sincos(BEIDOU_GEO_TILT)
    R_z = @SMatrix [cφ sφ 0.0; -sφ cφ 0.0; 0.0 0.0 1.0]
    R_z_dot = Ω_dot_e * @SMatrix [-sφ cφ 0.0; -cφ -sφ 0.0; 0.0 0.0 0.0]
    R_x = @SMatrix [1.0 0.0 0.0; 0.0 ct st; 0.0 -st ct]
    tilted_position = R_x * position
    (
        position = R_z * tilted_position,
        velocity = R_z_dot * tilted_position + R_z * (R_x * velocity),
    )
end

# ---- Group delay ------------------------------------------------------------
#
# BeiDou references its broadcast clock polynomial to the B3I signal, on both the
# legacy and the BDS-3 messages, and then publishes a group delay per signal to
# get from there to the signal a range was actually generated on
# (BDS-SIS-ICD-B1I-3.0 §5.2.4.10, BDS-SIS-ICD-B1C-1.0 §7.6 and its B2a/B2b
# counterparts). B3I itself therefore needs no correction at all — the one
# ranging signal in this package whose clock is already referred to it.
#
# The data components carry an additional inter-signal correction on top of their
# band's pilot group delay (`ISC_B1Cd`, `ISC_B2ad`), exactly like the GPS L1C and
# L5 ISCs above.
#
# `group_delay_term` appears on B2a's terms and on no others, which is not an
# oversight but the decoder's gating read back. `is_decoding_completed_for_positioning`
# requires a signal's own single-band correction wherever it rides in the block the
# gate already waits for, and excludes it where requiring it would mean waiting on a
# message the ICD does not schedule:
#
#   - D1/D2 gates `T_GD1`, the correction for the one band it is decoded on that has
#     a ranging signal here (not `T_GD2` — see the B1I method below); B1C gates
#     `T_GD_B1Cp` and `ISC_B1Cd` (bits
#     546-569 of the same CRC-protected subframe 2 as the ephemeris); B2b gates
#     `T_GD_B2bI` (same MT30 as the clock). All read unwrapped — a `nothing` there
#     means a broken decoder, and taking it as zero would hide that as a metre of
#     range.
#   - B2a alone is excluded: `T_GD_B2ap`, `ISC_B2ad` and `T_GD_B1Cp` ride in MT30
#     only, while the clock and IODC are in all of MT30-34, and BDS-SIS-ICD-B2a-1.0
#     §6.2 schedules only the MT10/11 pair. Gating on MT30 would gate on an interval
#     nothing bounds, so those really can be absent and are taken as zero.

# Legacy D1/D2 on B1I: T_GD1 is the B1I group delay differential. The message also
# carries T_GD2, for BDS-2's B2I, and nothing reads it: B2I has no signal definition
# in GNSSSignals and will not get one (JuliaGNSS/GNSSSignals.jl#156 — 13 of the 15
# remaining BDS-2 satellites were decommissioned in April 2026, leaving two IGSO
# transmitters and no successor, since B2a and B2b replaced it on BDS-3). So this is
# a field with no consumer rather than one awaiting a ranging signal.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
    ::GNSSSignals.BeiDouB1I,
    t,
) = t - decoder.data.T_GD1
# ... and on B3I: none, the clock is already B3I-referenced.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
    ::GNSSSignals.BeiDouB3I,
    t,
) = t

# B-CNAV1 (B1C) carries the B1C pilot group delay plus the data-component ISC. It
# also broadcasts the B2a pilot's group delay, but no method reads it: a B2a range is
# generated on the B2a component and so carries a B2a decoder, never a B1C one.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB1CData},
    ::GNSSSignals.BeiDouB1C_D,
    t,
) = t - decoder.data.T_GD_B1Cp - decoder.data.ISC_B1Cd
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB1CData},
    ::GNSSSignals.BeiDouB1C_P,
    t,
) = t - decoder.data.T_GD_B1Cp
# B-CNAV2 (B2a) mirrors it: the B2a pilot group delay plus the B2a data ISC. It too
# carries the other band's pilot delay (`T_GD_B1Cp`) that nothing here reads.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2aData},
    ::GNSSSignals.BeiDouB2aI,
    t,
) =
    t - group_delay_term(decoder.data.T_GD_B2ap) -
    group_delay_term(decoder.data.ISC_B2ad)
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2aData},
    ::GNSSSignals.BeiDouB2aQ,
    t,
) = t - group_delay_term(decoder.data.T_GD_B2ap)
# B-CNAV3 (B2b) carries only its own: B2b_I has no pilot and no sibling ISC.
correct_by_group_delay(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouB2bData},
    ::GNSSSignals.BeiDouB2bI,
    t,
) = t - decoder.data.T_GD_B2bI

# There is deliberately no catch-all method, and the methods above cover exactly the
# ranging signals of the band each message is decoded on — data component and pilot.
# Anything else is a cross-band pairing, which cannot arise: a range is generated on
# the band its own data component was decoded from. Returning `t` for such a call, as
# an earlier version of this file did, is the one answer that cannot be told apart
# from a correct zero correction; a `MethodError` says what actually happened. GPS has
# never had such a fallback either.
#
# D1/D2 spans two bands because one message is broadcast on both B1I and B3I, so both
# of those are same-band pairings with their own decoder — as are GPS CNAV's L5 and
# L2C signals, for the same reason.

# ---- BDT–GPST offset (BGTO) --------------------------------------------------
#
# Nothing here. The BeiDou counterpart of the Galileo GGTO used to need ~90 lines
# to normalise three broadcast shapes — D1/D2's epochless two-term polynomial,
# B2a/B2b's `GNSS_ID`-gated three-term one, and B1C's dictionary keyed by GNSS ID
# — plus the `GNSS_ID == 1` trap (`0` means "unavailable", not GPS) and, on the
# Galileo side, the 6-bit reference week. GNSSDecoder's `get_time_offset` answers
# all of it behind one call, so `gpst_offset_available` and `calc_gpst_offset` in
# `src/sat_time.jl` are one method each and neither mentions a constellation.

# ---- Ionosphere -------------------------------------------------------------

# The legacy D1/D2 message broadcasts an eight-coefficient Klobuchar-style set,
# but BDS-SIS-ICD-B1I-3.0 §5.2.4.7 prescribes BeiDou's own variant of the
# algorithm, not a restatement of IS-GPS-200 — exact pierce-point geometry,
# geographic rather than geomagnetic latitude, a period clamped from above too,
# and the delay referenced to B1I rather than L1 — so the set gets a type of its
# own (see `BeiDouKlobucharParams` and `beidou_klobuchar_group_delay`). The
# BDS-3 B-CNAV messages instead broadcast BDGIM (`α_1…α_9`), a different model
# that is not implemented here — those containers therefore report no Klobuchar
# set.
function klobuchar_params(
    decoder::GNSSDecoder.GNSSDecoderState{<:GNSSDecoder.BeiDouDNAVData},
)
    d = decoder.data
    any(isnothing, (d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)) &&
        return nothing
    return BeiDouKlobucharParams(d.α_0, d.α_1, d.α_2, d.α_3, d.β_0, d.β_1, d.β_2, d.β_3)
end
