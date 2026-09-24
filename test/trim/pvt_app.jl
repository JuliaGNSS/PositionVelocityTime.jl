# Entry point of the `juliac --trim=safe` check (see `check.jl`): the solves of the
# package's precompile workload, on its own fixture satellites, printed so the
# trimmed executable's output can be compared against a regular Julia session.
using PositionVelocityTime, GNSSSignals, GNSSDecoder
using PositionVelocityTime:
    SignalGroup, _precompile_states, _precompile_cnav, _precompile_cnav2, _precompile_fnav,
    _PRECOMPILE_GPS_L1CA_STATES, _PRECOMPILE_GALILEO_E1B_STATES
using Unitful: Hz

gps(system, make_data = identity) = SignalGroup(
    system, _precompile_states(system, _PRECOMPILE_GPS_L1CA_STATES, make_data, GPSL1CA()))
galileo(system, make_data = identity) = SignalGroup(
    system, _precompile_states(system, _PRECOMPILE_GALILEO_E1B_STATES, make_data, GalileoE1B()))

# The fixtures hold no BeiDou satellites, so this group's is undecoded and drops out as
# unhealthy — the solution is that of the other groups — but the solve is compiled, and
# so verified, for the BeiDou message family and its ionospheric model.
beidou(system) = SignalGroup(system, [SatelliteState(;
    decoder = GNSSDecoderState(system, 1), system, code_phase = 0.0, carrier_doppler = 0.0Hz)])

# One `print` per value: a long `print(io, xs...)` is not specialised on its
# arguments' types, which leaves the call dynamic.
function report(io, name, pvt)
    print(io, name, ":")
    for x in (pvt.position..., pvt.velocity..., pvt.time.second, pvt.time.fraction,
        pvt.dop.GDOP, length(pvt.sats))
        print(io, " ", x)
    end
    for (system, bias) in pvt.inter_system_biases
        print(io, " ", get_time_system_id(system), "=", bias.val)
    end
    println(io)
end

function solve(io, name, groups)
    cold = calc_pvt(groups; approximate_year = 2021)
    report(io, name * " cold", cold)
    report(io, name * " warm", calc_pvt(groups, cold; approximate_year = 2021))
    report(io, name * " uncorrected", calc_pvt(groups; approximate_year = 2021,
        enable_ionospheric_correction = false, enable_tropospheric_correction = false))
    # The in-place solve, warm-started from (and overwriting) its own solution.
    in_place = calc_pvt!(PVTSolution(), PVTWorkspace(), groups, cold; approximate_year = 2021)
    report(io, name * " in place", calc_pvt!(in_place, PVTWorkspace(), groups, in_place;
        approximate_year = 2021))
end

function (@main)(ARGS)
    io = Core.stdout
    l1ca = gps(GPSL1CA())
    e1b = galileo(GalileoE1B())
    solve(io, "GPS L1 C/A", l1ca)
    solve(io, "GPS L5", gps(GPSL5I(), _precompile_cnav))
    solve(io, "GPS L2C", gps(GPSL2CM(), _precompile_cnav))
    solve(io, "GPS L1C", gps(GPSL1C_D(), _precompile_cnav2))
    solve(io, "Galileo E1B", e1b)
    solve(io, "Galileo E5a", galileo(GalileoE5aI(), _precompile_fnav))
    solve(io, "GPS + Galileo", (gps = l1ca, galileo = e1b))
    solve(io, "GPS + Galileo + BeiDou", (gps = l1ca, galileo = e1b,
        b1i = beidou(BeiDouB1I()), b1c = beidou(BeiDouB1C_D()), b2a = beidou(BeiDouB2aI())))
    return 0
end
