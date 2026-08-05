# ===========================================================================
#  A-priori measurement uncertainty — C/N₀ and elevation → variance
#
#  Ordinary least squares gives every satellite the same say in the fix, which is
#  only right when every pseudorange carries the same noise. It does not: code
#  measurement noise is a strong function of C/N₀ (the tracking loops' thermal
#  jitter), and the residual model errors — multipath and whatever the atmospheric
#  models leave behind — grow towards the horizon. A marginal satellite therefore
#  arrives with roughly an order of magnitude more noise than a strong one, and with
#  five or six satellites it can visibly drag a fix the others would have pinned down.
#
#  The standard remedy is weighted least squares with *a-priori* weights `w = 1/σ²`
#  from a variance model, which is what this file provides: one variance for the
#  pseudorange (the position solve) and one for the range rate (the velocity solve),
#  since code and carrier noise follow different loops. The weights are a-priori by
#  design — residual-based reweighting (IRLS) would be the wrong tool here, because a
#  four-to-six-satellite epoch against 4+ unknowns has barely any degrees of freedom,
#  so its residuals are shaped by the geometry and by the fit itself rather than by
#  measurement quality, and it would happily down-weight a perfectly good satellite.
#
#  Only the *ratios* of the variances affect the estimate; their absolute scale shows
#  up solely in the reported formal accuracy (see [`FormalAccuracy`](@ref)).
# ===========================================================================

# Non-thermal part of the user-equivalent range error: residual satellite
# ephemeris/clock error plus what the ionospheric and tropospheric models leave
# behind, as a single elevation-independent floor. This is the term that keeps the
# weights finite and the weight spread bounded when one satellite reports a very high
# C/N₀ — a strong satellite is not infinitely precise, it is limited by this.
const _RESIDUAL_UERE_SIGMA = 1.5  # m

# Elevation-mapped term `a/sin(El)`: multipath and residual atmosphere, both of which
# grow along the lengthening slant path towards the horizon. `a` is the zenith value,
# so this contributes 0.5 m straight up and ~5.7 m at 5° elevation. The RTKLIB
# baseline model has the same `a²+b²/sin²El` shape (`varerr` in `pntpos.c`), with the
# code-noise term there scaled up from the carrier-phase one instead of derived from
# C/N₀ as below.
const _ELEVATION_SIGMA = 0.5  # m

# Code-loop parameters behind the C/N₀ term: a 1 Hz single-sided noise bandwidth and a
# one-chip early–late spacing, the conventional reference configuration for the
# coherent-DLL thermal-jitter expression (Kaplan & Hegarty ch. 5, Groves ch. 9)
#
#     σ_code = (c / f_chip) · √(B_n · d / 2) · 10^(−CN0[dBHz]/20)   metres
#
# whose square is the `b²·10^(−CN0/10)` term of the standard model — with `b` derived
# per signal from its chipping rate rather than tuned, so a wider-band signal is
# credited its narrower correlation peak automatically. For GPS L1 C/A
# (`c/f_chip` ≈ 293 m) this gives ≈ 1.2 m at 45 dBHz and ≈ 6.5 m at 30 dBHz, the
# order-of-magnitude spread that motivates weighting at all. A receiver whose loops
# differ can bypass the model with `SatelliteState.pseudorange_variance`.
const _CODE_LOOP_BANDWIDTH = 1.0   # Hz
const _EARLY_LATE_SPACING = 1.0    # chips

# Robustness floor and ceiling on σ. A C/N₀ estimate inherits its estimator's bias and
# spread, and the elevation term diverges towards the horizon, so σ is clamped: no
# single reading can hand one satellite near-infinite leverage (the floor), and a
# satellite whose C/N₀ reads absurdly low is heavily down-weighted but never fully
# discarded (the ceiling). The bounds also bound the conditioning cost of weighting:
# `cond(√W·H) ≤ (σ_max/σ_min)·cond(H)`, i.e. at most two orders of magnitude here.
const _PSEUDORANGE_SIGMA_BOUNDS = (0.5, 50.0)  # m

# Uniform σ assumed when *no* satellite carries any uncertainty information. The solve
# is then plain OLS (uniform weights change nothing), so this value never moves the
# estimate — it only sets the scale of the reported [`FormalAccuracy`](@ref), which in
# that case is the familiar "DOP × UERE". 5 m is a nominal single-frequency UERE.
const _NOMINAL_PSEUDORANGE_SIGMA = 5.0  # m

# Range-rate (Doppler) model, for the velocity/clock-drift solve. Doppler noise follows
# the carrier loop, not the DLL, so it gets its own variance: the coherent FLL
# thermal-jitter expression
#
#     σ_rangerate = λ / (2π·T) · √(4·B_n / (C/N₀))   metres/second
#
# with a 1 Hz loop bandwidth and a 20 ms coherent integration (one GPS L1 C/A data
# symbol — the longest integration available before bit sync is exploited, and the
# common choice once it is). For L1 that is ≈ 0.10 m/s at 30 dBHz and ≈ 0.02 m/s at
# 45 dBHz. The floor stands for what is left after the loop: satellite-velocity and
# clock-drift model error.
const _RANGE_RATE_SIGMA_FLOOR = 0.01     # m/s
const _CARRIER_LOOP_BANDWIDTH = 1.0      # Hz
const _CARRIER_INTEGRATION_TIME = 0.02   # s
const _RANGE_RATE_SIGMA_BOUNDS = (0.005, 5.0)  # m/s

"""
    _is_usable_cn0(cn0) -> Bool

Whether `cn0` is a C/N₀ reading the variance models can use: present, finite and
positive. `nothing` is the "not reported" case, and a non-positive level is not a
plausible measurement of a satellite that is being tracked — `0.0dBHz` is in
particular what `Tracking`'s moment-based estimator returns before it has integrated
a single prompt — so both are treated as *unknown* rather than as an unusably weak
signal, leaving the C/N₀ term out of the variance instead of driving it to the
ceiling.
"""
_is_usable_cn0(::Nothing) = false
_is_usable_cn0(cn0) = isfinite(ustrip(cn0)) && cn0 > 0.0dBHz

"""
    _linear_cn0(cn0) -> Float64

The carrier-to-noise-density ratio as a linear ratio in Hz (e.g. `45.0dBHz` →
`3.16e4`), the form the thermal-jitter expressions take it in. Unitful's
logarithmic-unit conversion does the `10^(CN0/10)` itself.
"""
_linear_cn0(cn0) = ustrip(Hz, uconvert(Hz, cn0))

"""
    has_measurement_uncertainty(state::SatelliteState) -> Bool

Whether `state` carries the information needed to weight it: an explicit
`pseudorange_variance`, or a usable `cn0` (present, finite and positive — see
`_is_usable_cn0`).
[`calc_pvt`](@ref) weights an epoch's solve only when at least one satellite does;
with none it stays ordinary least squares, so a caller that reports neither keeps
exactly the unweighted behaviour.
"""
has_measurement_uncertainty(state::SatelliteState) =
    !isnothing(state.pseudorange_variance) || _is_usable_cn0(state.cn0)

"""
    pseudorange_variance(system, cn0, elevation) -> Float64
    pseudorange_variance(state::SatelliteState, elevation) -> Float64

A-priori variance of a pseudorange measurement in m², the reciprocal of its
least-squares weight. `system` is the ranging signal (its chipping rate sets the
code-noise scale), `cn0` its carrier-to-noise-density ratio as a `dBHz` level (or
`nothing` when not reported) and `elevation` the satellite elevation in radians (or
`nothing` when no position is available yet, as at a cold start):

    σ² = σ²_UERE + (a / sin El)² + (c/f_chip)² · (B_n·d/2) · 10^(−CN0/10)

The three terms are the residual non-thermal UERE, the elevation-mapped
multipath/atmosphere term, and the code-loop thermal jitter; each of the latter two
is dropped when its input is `nothing`. σ is clamped to
$(_PSEUDORANGE_SIGMA_BOUNDS) m, so no single reading can dominate the fix and none is
silently discarded. Elevations below `_LOW_ELEVATION_THRESHOLD` (≈ 2.87°, including
negative ones) are treated as that bound — the same guard the tropospheric mapping
uses, and the reason a grazing satellite ends up de-weighted rather than trusted.

The `SatelliteState` form is what [`calc_pvt`](@ref) calls: it returns the state's
own `pseudorange_variance` when the caller supplied one (the model is then bypassed
entirely, units stripped to m²) and otherwise evaluates the model above. A supplied
variance is clamped to the same bounds, so a degenerate one (zero, negative, `NaN`)
cannot turn into an infinite weight and take over the fix.

$SIGNATURES
"""
function pseudorange_variance(system, cn0, elevation)
    σ² = _RESIDUAL_UERE_SIGMA^2
    if !isnothing(elevation)
        sin_el = sin(max(elevation, _LOW_ELEVATION_THRESHOLD))
        σ² += (_ELEVATION_SIGMA / sin_el)^2
    end
    if _is_usable_cn0(cn0)
        chip_length = SPEEDOFLIGHT / ustrip(Hz, get_code_frequency(system))
        σ² +=
            chip_length^2 * _EARLY_LATE_SPACING * _CODE_LOOP_BANDWIDTH / 2 /
            _linear_cn0(cn0)
    end
    _bound_pseudorange_variance(σ²)
end

pseudorange_variance(state::SatelliteState, elevation) =
    isnothing(state.pseudorange_variance) ?
    pseudorange_variance(state.system, state.cn0, elevation) :
    _bound_pseudorange_variance(ustrip(m^2, state.pseudorange_variance))

# `clamp` on the variance rather than on σ, and `NaN`-safe: `clamp` would pass a `NaN`
# straight through (every comparison is false), and a `NaN` weight poisons the whole
# solve, so map it to the ceiling — the treatment a meaningless reading deserves.
_bound_pseudorange_variance(σ²) =
    isnan(σ²) ? _PSEUDORANGE_SIGMA_BOUNDS[2]^2 :
    clamp(σ², _PSEUDORANGE_SIGMA_BOUNDS[1]^2, _PSEUDORANGE_SIGMA_BOUNDS[2]^2)

"""
    range_rate_variance(system, cn0) -> Float64
    range_rate_variance(state::SatelliteState) -> Float64

A-priori variance of a Doppler-derived range-rate measurement in (m/s)², the
reciprocal of its weight in the velocity and clock-drift solve:

    σ² = σ²_floor + (λ / 2πT)² · 4·B_n · 10^(−CN0/10)

`λ` is the carrier wavelength of `system` and the C/N₀ term is the coherent
frequency-loop thermal jitter (see the constants above); it is dropped when `cn0` is
not usable, which leaves the variance uniform and the velocity solve unweighted. σ is
clamped to $(_RANGE_RATE_SIGMA_BOUNDS) m/s, for the reasons given for the pseudorange
bounds. Unlike the pseudorange, there is no per-satellite override field: the Doppler
is not the measurement a caller is likely to have its own error model for.

$SIGNATURES
"""
function range_rate_variance(system, cn0)
    σ² = _RANGE_RATE_SIGMA_FLOOR^2
    if _is_usable_cn0(cn0)
        λ = SPEEDOFLIGHT / ustrip(Hz, get_center_frequency(system))
        σ² +=
            (λ / (2π * _CARRIER_INTEGRATION_TIME))^2 * 4 * _CARRIER_LOOP_BANDWIDTH /
            _linear_cn0(cn0)
    end
    clamp(σ², _RANGE_RATE_SIGMA_BOUNDS[1]^2, _RANGE_RATE_SIGMA_BOUNDS[2]^2)
end

range_rate_variance(state::SatelliteState) = range_rate_variance(state.system, state.cn0)

"""
    predict_pseudorange_variances(ξ, states, sat_positions) -> Vector{Float64}

Per-satellite pseudorange variances (m²) evaluated about the least-squares state
vector `ξ = [x, y, z, tc₁, …]`, whose first three elements give the user ECEF
position the elevations are taken from. Like `predict_atmospheric_delays`
this only needs a metre-accurate position — ∂σ/∂position over that uncertainty is
utterly negligible (a 15 m shift moves an elevation by ~1e-5°) — so one prediction
about the previous or bootstrap solution serves the solve that follows, and the ENU
transform is built once per epoch rather than per satellite.

$SIGNATURES
"""
function predict_pseudorange_variances(ξ, states, sat_positions)
    enu_from_ecef = ENUfromECEF(ECEF(ξ[1], ξ[2], ξ[3]), wgs84)
    map(states, sat_positions) do state, sat_pos
        elevation, _ = _elevation_azimuth(enu_from_ecef, sat_pos)
        pseudorange_variance(state, elevation)
    end
end
