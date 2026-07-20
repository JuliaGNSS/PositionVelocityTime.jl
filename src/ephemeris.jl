"""
    Ephemeris

Concrete (`Float64`-only) snapshot of one satellite's decoded broadcast orbit — the
Keplerian elements, their perturbation corrections and the propagation constants —
extracted once from a `GNSSDecoderState` and then used by the orbit kernels
([`calc_satellite_position_and_velocity`](@ref), `calc_eccentric_anomaly`,
`calc_relativistic_correction`). The broadcast clock model lives in the companion
[`ClockModel`](@ref), which composes an `Ephemeris` (the relativistic clock term is
orbit-dependent); a decoder whose clock is not (yet) decoded can still be propagated.

The decoder's data fields are all `Union{Nothing, T}` (they fill in as the navigation
message decodes), and the Kepler propagation reads ~20 of them in long arithmetic
chains — too many small unions for the compiler to split, so every intermediate gets
boxed. Extracting the fields into this concrete struct behind a function barrier makes
the downstream math allocation-free and several times faster; `calc_pvt` performs the
extraction once per satellite per solve.

The two broadcast ephemeris parameterisations are normalised at extraction (see
[`orbital_terms`](@ref)): the directly-broadcast Keplerian elements (GPS LNAV, Galileo
I/NAV, F/NAV) and the quasi-Keplerian CNAV/CNAV-2 deltas both reduce to a semi-major
axis `A` (+ rate `A_dot`), a mean motion `n_0` at the reference epoch (+ half-rate
`n_dot_half`, zero for the directly-broadcast case) and an effective `Ω_dot`, so a
single set of kernels serves all four navigation messages. The propagation constants
(Earth rotation rate `Ω_dot_e`, relativistic `F`) are copied from the decoder so the
kernels need only this struct.

Constructing from a decoder that has not decoded all the orbit fields throws an
`ArgumentError`; the fields required are a subset of the ones guaranteed by
`is_decoding_completed_for_positioning`.
"""
struct Ephemeris
    t_0e::Float64
    A::Float64
    sqrt_A::Float64
    A_dot::Float64
    n_0::Float64
    n_dot_half::Float64
    e::Float64
    M_0::Float64
    ω::Float64
    i_0::Float64
    i_dot::Float64
    Ω_0::Float64
    Ω_dot::Float64
    C_us::Float64
    C_uc::Float64
    C_rs::Float64
    C_rc::Float64
    C_is::Float64
    C_ic::Float64
    Ω_dot_e::Float64
    F::Float64
end

function Ephemeris(decoder::GNSSDecoder.GNSSDecoderState)
    data = decoder.data
    constants = decoder.constants
    orb = orbital_terms(data, constants.μ)
    # `something` doubles as the Union{Nothing,T} → T function barrier: each argument
    # is concretely inferred (so nothing is boxed) and a missing field throws instead
    # of poisoning the arithmetic downstream.
    Ephemeris(
        Float64(something(data.t_0e)),
        orb.A,
        orb.sqrt_A,
        orb.A_dot,
        orb.n_0,
        orb.n_dot_half,
        something(data.e),
        something(data.M_0),
        something(data.ω),
        something(data.i_0),
        something(data.i_dot),
        something(data.Ω_0),
        orb.Ω_dot,
        something(data.C_us),
        something(data.C_uc),
        something(data.C_rs),
        something(data.C_rc),
        something(data.C_is),
        something(data.C_ic),
        constants.Ω_dot_e,
        constants.F,
    )
end

"""
    ClockModel

Concrete snapshot of one satellite's broadcast clock correction: the clock polynomial
(`a_f0`, `a_f1`, `a_f2` at `t_0c`), the signal-dependent group delay folded into the
constant `group_delay` (seconds, subtracted from the clock correction), and the
[`Ephemeris`](@ref) the relativistic clock term is computed from. Used by the clock
kernels (`correct_clock`, `calc_satellite_clock_drift`, `calc_corrected_time`);
`calc_pvt` extracts one per satellite per solve and reuses its `ephemeris` for the
orbit propagation, so the `Union{Nothing, …}` decoder fields are read exactly once
(see [`Ephemeris`](@ref) for why that matters).
"""
struct ClockModel
    ephemeris::Ephemeris
    a_f0::Float64
    a_f1::Float64
    a_f2::Float64
    t_0c::Float64
    group_delay::Float64
end

# `correct_by_group_delay(decoder, system, Δt)` computes `Δt - group_delay` (with the
# per-message/per-signal dispatch and missing-term handling), so evaluating it at zero
# recovers the constant. Reusing it keeps this extraction and the directly-tested
# group-delay logic a single source of truth.
function ClockModel(decoder::GNSSDecoder.GNSSDecoderState, system::AbstractGNSSSignal)
    data = decoder.data
    ClockModel(
        Ephemeris(decoder),
        something(data.a_f0),
        something(data.a_f1),
        something(data.a_f2),
        Float64(something(data.t_0c)),
        -correct_by_group_delay(decoder, system, 0.0),
    )
end

"""
    orbital_terms(data, μ) -> (; A, sqrt_A, A_dot, n_0, n_dot_half, Ω_dot)

Normalise the Keplerian terms that differ between the directly-broadcast ephemerides
(GPS LNAV `GPSL1CAData`, Galileo I/NAV `GalileoE1BData` / F/NAV `GalileoE5aData`) and
the quasi-Keplerian GPS CNAV/CNAV-2 ephemerides (`GPSCNAVData` for L5 and L2C,
`GPSL1C_DData` for L1C):

- `A`: semi-major axis at `t_0e` (m); `sqrt_A²` for the directly-broadcast case,
  `A_REF + ΔA` for CNAV/CNAV-2
- `sqrt_A`: its square root (m^½); the broadcast `sqrt_A` directly (no round-trip)
  for the directly-broadcast case, `√A` for CNAV/CNAV-2 (which carry no `sqrt_A`)
- `A_dot`: its rate (m/s; `0` for the directly-broadcast case)
- `n_0`: mean motion at `t_0e` (rad/s), including the broadcast correction
  (`Δn` / `Δn_0`)
- `n_dot_half`: half the mean-motion rate `½ Δṅ_0` (rad/s²; `0` for the
  directly-broadcast case), so the mean motion at time-from-ephemeris `t_k` is
  `n_0 + n_dot_half · t_k` for every message type
- `Ω_dot`: effective rate of right ascension (rad/s); broadcast directly, or
  `Ω̇_REF + ΔΩ̇` for CNAV/CNAV-2
"""
function orbital_terms(data::GNSSDecoder.AbstractGNSSData, μ)
    sqrt_A = something(data.sqrt_A)
    (
        A = sqrt_A^2,
        sqrt_A = sqrt_A,
        A_dot = 0.0,
        n_0 = sqrt(μ) / sqrt_A^3 + something(data.Δn),
        n_dot_half = 0.0,
        Ω_dot = something(data.Ω_dot),
    )
end
function orbital_terms(data::GPSModernNavData, μ)
    # Quasi-Keplerian reference values from the CNAV user algorithm (IS-GPS-200N;
    # identical in IS-GPS-705J and IS-GPS-800J); the broadcast fields are deltas off
    # these.
    A_REF = 26_559_710.0        # m
    Ω_dot_REF = -2.6e-9 * π     # rad/s (-2.6e-9 semicircles/s)
    A = A_REF + something(data.ΔA)
    (
        A = A,
        sqrt_A = sqrt(A),
        A_dot = something(data.A_dot),
        n_0 = sqrt(μ / A^3) + something(data.Δn_0),
        n_dot_half = 0.5 * something(data.Δn_0_dot),
        Ω_dot = Ω_dot_REF + something(data.ΔΩ_dot),
    )
end
