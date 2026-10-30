# Precompile workload (PrecompileTools). The first `calc_pvt` of a session costs 4.7 s of
# compilation on a workstation and three to four times that on an embedded ARM host, paid
# by a live receiver just as its first fix is due (GNSSReceiver.jl#107).
#
# Only the collection pass specialises on the navigation-data type; the solver behind the
# measurements is compiled once (see collection.jl). So the workload is one solve per data
# type (GPS LNAV, CNAV on L5I and L2CM, CNAV-2, Galileo I/NAV and F/NAV), plus once each
# the paths a data type does not vary: warm start, corrections disabled, a
# `Vector`-backed group and a mixed constellation.
#
# The satellites are the test fixtures (`test/fixtures.jl`, Aachen, 2021-05-31), recorded
# on GPS L1 C/A and Galileo E1B only; the other data types are expressed from them (see
# `_precompile_states`), and every configuration solves to the fixture's position.
using PrecompileTools: @setup_workload, @compile_workload

# The CNAV reference values (IS-GPS-200N) the broadcast deltas are off.
const _PRECOMPILE_A_REF = 26_559_710.0
const _PRECOMPILE_Ω_DOT_REF = -2.6e-9 * π

# GPS LNAV ephemeris → the quasi-Keplerian CNAV / CNAV-2 form. `something` narrows the
# fields, keeping these keyword calls static for the trim check (`test/trim`).
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
# (IS-GPS-800J §3.5.3), which GNSSDecoder's `get_broadcast_time_of_week` reassembles.
_precompile_cnav2(d) = _precompile_modern_nav(
    GNSSDecoder.GPSCNAV2Data,
    d;
    ITOW = something(d.TOW) ÷ 7200,
    toi = (something(d.TOW) % 7200) ÷ 18,
    WN = 2160,
    l1c_health = false,
)

# Galileo F/NAV carries the same ephemeris as I/NAV, with the E5a health and
# group delay in place of the E1B/E5b ones.
_precompile_fnav(d) = GNSSDecoder.GalileoFNAVData(;
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
    _precompile_states(signal, table, make_data, source_signal)

The fixture satellites as [`SatelliteState`](@ref)s on `signal`, their navigation data
converted by `make_data`. The observables, recorded on `source_signal`, are re-expressed
on `signal` preserving the transmit epoch (`TOW + n_bits/bit_rate + code_phase/chip_rate`,
see `calc_uncorrected_transmit_time`): the sub-TOW offset is converted as a whole, then
split into the target's bits and code phase, since a shorter bit carries part
of the code phase in its bit count.
"""
function _precompile_states(signal, table, make_data, source_signal)
    map(table) do s
        data = make_data(s.data)
        target = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.GNSSDecoderState(signal, s.prn);
            raw_data = data,
            data = data,
        )
        source = GNSSDecoder.GNSSDecoderState(
            GNSSDecoder.GNSSDecoderState(source_signal, s.prn);
            raw_data = s.data,
            data = s.data,
        )
        source_rate = ustrip(Hz, get_bit_rate(source))
        source_chips = ustrip(Hz, get_chip_rate(source_signal))
        offset =
            s.num_bits / source_rate +
            mod(s.code_phase, source_chips / source_rate) / source_chips
        rate = ustrip(Hz, get_bit_rate(target))
        chips = ustrip(Hz, get_chip_rate(signal))
        num_bits = floor(Int, offset * rate)
        SatelliteState(;
            decoder = GNSSDecoder.GNSSDecoderState(
                target;
                num_bits_since_time_reference = num_bits,
            ),
            signal,
            code_phase = (offset - num_bits / rate) * chips,
            carrier_doppler = s.carrier_doppler,
            carrier_phase = s.carrier_phase * ustrip(Hz, get_center_frequency(signal)) /
                            ustrip(Hz, get_center_frequency(source_signal)),
        )
    end
end

@setup_workload begin
    gps(signal, make_data = identity) =
        _precompile_states(signal, _PRECOMPILE_GPS_L1CA_STATES, make_data, GPSL1CA())
    galileo(signal, make_data = identity) =
        _precompile_states(signal, _PRECOMPILE_GALILEO_E1B_STATES, make_data, GalileoE1B())
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
        # The paths that do not vary with the data type, compiled once on GPS L1 C/A.
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
        # The in-place entry point; only its wrapper and the unsolvable-epoch copy are new.
        workspace = PVTWorkspace()
        in_place = calc_pvt!(PVTSolution(), workspace, mixed, mixed_pvt;
            approximate_year = 2021)
        in_place = calc_pvt!(in_place, workspace, mixed, in_place; approximate_year = 2021)
        calc_pvt!(PVTSolution(), workspace, SignalGroup(GPSL1CA(), l1ca[1:3]), in_place;
            approximate_year = 2021)
    end
end
