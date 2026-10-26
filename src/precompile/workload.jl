# Precompile workload (PrecompileTools). The first `calc_pvt` of a session costs
# 4.7 s of compilation on a workstation and three to four times that on an
# embedded ARM host, and a live receiver pays it at the moment it has its first
# four ephemerides — the position fix it has been waiting a minute for arrives
# seconds late while every tracking loop sits unattended (GNSSReceiver.jl#107).
#
# What is specialised on the navigation-data type is the *collection* pass
# (`collect_measurements`): the decoder accessors, the propagator, the clock
# polynomial and the broadcast time offsets. The solver behind the flat
# `SatelliteMeasurement` row is one body shared by every data type and every
# constellation mix, so it is compiled by the first solve here and reused by all
# the others — which is what lets this workload be a list of data types rather
# than a list of their combinations.
#
# So: one solve per type this package dispatches on — GPS LNAV (`GPSL1CAData`),
# GPS CNAV on both its signals (`GPSCNAVData` on L5I and L2CM), GPS CNAV-2
# (`GPSL1C_DData`) and Galileo I/NAV and F/NAV (`GalileoINAVData`,
# `GalileoE5aData`) — plus, once, the solver paths a data type does not vary:
# warm start, corrections disabled, a `Vector`-backed group and a two-group
# mixed constellation.
#
# The satellites are the test suite's own fixtures (`test/fixtures.jl`): decoder
# states captured over Aachen on 2021-05-31 with the code phases, carrier
# Dopplers and phases observed alongside them. Only GPS L1 C/A and Galileo E1B
# were recorded, so the other data types are expressed from those — see
# `_precompile_states`. Every configuration below solves to the fixture's own
# position, which is what makes this a real solve rather than a shape.
using PrecompileTools: @setup_workload, @compile_workload

# Quasi-Keplerian reference values of the CNAV user algorithm; the broadcast
# fields are deltas off these (IS-GPS-200N, mirrored in `orbital_elements`).
const _PRECOMPILE_A_REF = 26_559_710.0
const _PRECOMPILE_Ω_DOT_REF = -2.6e-9 * π

# GPS LNAV ephemeris → the quasi-Keplerian CNAV / CNAV-2 form. The LNAV fields are
# `Union{Nothing,…}` until decoded; `something` narrows them to their value, which
# keeps these keyword calls statically dispatched in the trim check (`test/trim`).
_precompile_modern_nav(T, d; extra...) = T(;
    t_0e = something(d.t_0e),
    ΔA = something(d.sqrt_A)^2 - _PRECOMPILE_A_REF,
    A_dot = 0.0,
    Δn_0 = something(d.Δn),
    Δn_0_dot = 0.0,
    M_0 = something(d.M_0),
    e = something(d.e),
    ω = something(d.ω),
    Ω_0 = something(d.Ω_0),
    i_0 = something(d.i_0),
    ΔΩ_dot = something(d.Ω_dot) - _PRECOMPILE_Ω_DOT_REF,
    i_dot = something(d.i_dot),
    C_is = something(d.C_is),
    C_ic = something(d.C_ic),
    C_rs = something(d.C_rs),
    C_rc = something(d.C_rc),
    C_us = something(d.C_us),
    C_uc = something(d.C_uc),
    t_0c = something(d.t_0c),
    a_f0 = something(d.a_f0),
    a_f1 = something(d.a_f1),
    a_f2 = something(d.a_f2),
    T_GD = something(d.T_GD),
    extra...,
)

# CNAV carries the full week number where LNAV broadcasts it modulo 1024; the
# fixtures are GPS week 2160 (2021-05-31).
_precompile_cnav(d) = _precompile_modern_nav(
    GNSSDecoder.GPSCNAVData,
    d;
    TOW = something(d.TOW),
    WN = 2160,
    l1_health = false,
    l2_health = false,
    l5_health = false,
)

# CNAV-2 counts time as a two-hour interval plus an 18 s time of interval
# (IS-GPS-800J §3.5.3), which `get_tow` reassembles.
_precompile_cnav2(d) = _precompile_modern_nav(
    GNSSDecoder.GPSL1C_DData,
    d;
    ITOW = something(d.TOW) ÷ 7200,
    toi = (something(d.TOW) % 7200) ÷ 18,
    WN = 2160,
    l1c_health = false,
)

# Galileo F/NAV carries the same ephemeris as I/NAV, with the E5a health and
# group delay in place of the E1B/E5b ones.
_precompile_fnav(d) = GNSSDecoder.GalileoE5aData(;
    WN = d.WN,
    TOW = d.TOW,
    t_0e = d.t_0e,
    M_0 = d.M_0,
    e = d.e,
    sqrt_A = d.sqrt_A,
    Ω_0 = d.Ω_0,
    i_0 = d.i_0,
    ω = d.ω,
    i_dot = d.i_dot,
    Ω_dot = d.Ω_dot,
    Δn = d.Δn,
    C_uc = d.C_uc,
    C_us = d.C_us,
    C_rc = d.C_rc,
    C_rs = d.C_rs,
    C_ic = d.C_ic,
    C_is = d.C_is,
    t_0c = d.t_0c,
    a_f0 = d.a_f0,
    a_f1 = d.a_f1,
    a_f2 = d.a_f2,
    IOD_nav1 = d.IOD_nav1,
    E5a_SHS = GNSSDecoder.signal_ok,
    E5a_DVS = GNSSDecoder.navigation_data_valid,
    BGD_E1_E5a = d.BGD_E1_E5a,
)

"""
    _precompile_states(system, table, make_data, source_system)

The fixture satellites as [`SatelliteState`](@ref)s on `system`, with their
navigation data converted by `make_data`.

The observables were recorded on `source_system` (GPS L1 C/A or Galileo E1B).
Expressing them on another signal of the same satellite has to preserve the
transmit epoch, which [`calc_uncorrected_time`](@ref) builds as `TOW + n_bits/data_rate + (code_phase reduced to one data symbol)/chip_rate`. So the
sub-TOW offset is converted as a whole: taken off the source signal, then split
again into whole data symbols and a code phase in the target's own units.
Scaling the two terms independently would be wrong — a signal whose data symbol
is shorter has to carry part of the code phase in its bit count instead.
"""
function _precompile_states(system, table, make_data, source_system)
    map(table) do s
        data = make_data(s.data)
        target = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.GNSSDecoderState(system, s.prn);
            raw_data = data,
            data = data,
        )
        source = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.GNSSDecoderState(source_system, s.prn);
            raw_data = s.data,
            data = s.data,
        )
        source_rate = ustrip(Hz, get_data_frequency(source))
        source_chips = ustrip(Hz, get_code_frequency(source_system))
        offset =
            s.num_bits / source_rate +
            mod(s.code_phase, source_chips / source_rate) / source_chips
        rate = ustrip(Hz, get_data_frequency(target))
        chips = ustrip(Hz, get_code_frequency(system))
        num_bits = floor(Int, offset * rate)
        SatelliteState(;
            decoder = GNSSDecoder.GNSSDecoderState(
                target;
                num_bits_after_valid_syncro_sequence = num_bits,
            ),
            system,
            code_phase = (offset - num_bits / rate) * chips,
            carrier_doppler = s.carrier_doppler,
            carrier_phase = s.carrier_phase * ustrip(Hz, get_center_frequency(system)) /
                            ustrip(Hz, get_center_frequency(source_system)),
        )
    end
end

@setup_workload begin
    gps(system, make_data = identity) =
        _precompile_states(system, _PRECOMPILE_GPS_L1CA_STATES, make_data, GPSL1CA())
    galileo(system, make_data = identity) =
        _precompile_states(system, _PRECOMPILE_GALILEO_E1B_STATES, make_data, GalileoE1B())
    l1ca = gps(GPSL1CA())
    e1b = galileo(GalileoE1B())
    # A `Dictionary`-backed group keyed by PRN — the shape a receiver carrying its
    # satellites per signal hands over, and what Tracking.jl builds.
    group(signal, states) =
        SignalGroup(signal, Dictionary([s.decoder.prn for s in states], states))
    @compile_workload begin
        for signal_group in (
            group(GPSL1CA(), l1ca),
            group(GPSL5I(), gps(GPSL5I(), _precompile_cnav)),
            group(GPSL2CM(), gps(GPSL2CM(), _precompile_cnav)),
            group(GPSL1C_D(), gps(GPSL1C_D(), _precompile_cnav2)),
            group(GalileoE1B(), e1b),
            group(GalileoE5aI(), galileo(GalileoE5aI(), _precompile_fnav)),
        )
            calc_pvt(signal_group; approximate_year = 2021)
        end
        # The paths that do not vary with the navigation-data type, compiled once on
        # GPS L1 C/A and shared from there: the warm-start branch of the least-squares
        # solve, the correction-free branch, the `Vector`-backed group, and the
        # multi-group shape of a mixed-constellation epoch.
        l1ca_group = group(GPSL1CA(), l1ca)
        pvt = calc_pvt(l1ca_group; approximate_year = 2021)
        calc_pvt(l1ca_group, pvt; approximate_year = 2021)
        calc_pvt(
            l1ca_group;
            approximate_year = 2021,
            enable_ionospheric_correction = false,
            enable_tropospheric_correction = false,
        )
        calc_pvt(SignalGroup(GPSL1CA(), l1ca); approximate_year = 2021)
        mixed = (gps = l1ca_group, galileo = group(GalileoE1B(), e1b))
        mixed_pvt = calc_pvt(mixed; approximate_year = 2021)
        calc_pvt(mixed, mixed_pvt; approximate_year = 2021)
        # The in-place entry point shares the solver compiled above; only its own thin
        # wrapper and the solution copy of an unsolvable epoch are new.
        workspace = PVTWorkspace()
        in_place = calc_pvt!(PVTSolution(), workspace, mixed, mixed_pvt;
            approximate_year = 2021)
        in_place = calc_pvt!(in_place, workspace, mixed, in_place; approximate_year = 2021)
        calc_pvt!(PVTSolution(), workspace, SignalGroup(GPSL1CA(), l1ca[1:3]), in_place;
            approximate_year = 2021)
    end
end
