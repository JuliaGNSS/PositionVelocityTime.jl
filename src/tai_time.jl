"""
    TAITimes

A minimal, trim-safe TAI epoch, [`TAITime`](@ref), standing in for AstroTime.jl's
`TAIEpoch` as the type of `PositionVelocityTime.PVTSolution`'s `time`.

AstroTime's `__init__` cannot be verified by a `juliac --trim=safe` build (it registers
`Dates` format tokens and, through EarthOrientation and RemoteFiles, looks for data files
and download tools). `TAITime` stores an epoch as `TAIEpoch` does, so the two convert
exactly; with AstroTime loaded, a package extension provides `TAIEpoch(t)` and
`TAITime(epoch)`. Once AstroTime can be trimmed, this module can be dropped for it.
"""
module TAITimes

using Dates: Dates, DateTime, Millisecond

export TAITime

# J2000, as a TAI calendar label: the origin `second` counts from, as in AstroTime.
const J2000 = DateTime(2000, 1, 1, 12)

"""
    TAITime(second::Integer, fraction::Real)
    TAITime(datetime::Dates.DateTime)

An epoch on International Atomic Time (TAI): `second` whole seconds after J2000
(`2000-01-01T12:00:00` TAI) plus `fraction ∈ [0, 1)` of a second, AstroTime's `TAIEpoch`
layout. A `fraction` outside `[0, 1)` is normalised into it.

`TAITime(datetime)` reads a `DateTime` as a TAI calendar label (millisecond
resolution); `Dates.DateTime(t)` gives the TAI calendar label back, rounded to the
millisecond. The difference of two `TAITime`s is in seconds, `t + Δ` shifts `t` by
`Δ` seconds, and `≈` compares like AstroTime's epochs do.

# Fields
- `second::Int64`: whole TAI seconds since J2000.
- `fraction::Float64`: the fraction of the second, in `[0, 1)`.
"""
struct TAITime
    second::Int64
    fraction::Float64
    function TAITime(second::Integer, fraction::Real)
        whole = floor(Int64, fraction)
        rest = Float64(fraction) - whole
        # A fraction just below a whole number rounds up to exactly `1.0` in the
        # subtraction (`-1e-17 - (-1)` is `1.0` in `Float64`), which would escape the
        # `[0, 1)` invariant; carry it into the whole seconds instead.
        rest == 1.0 && (whole += 1; rest = 0.0)
        new(Int64(second) + whole, rest)
    end
end

function TAITime(datetime::DateTime)
    milliseconds = Dates.value(datetime - J2000)
    TAITime(fld(milliseconds, 1000), mod(milliseconds, 1000) / 1000)
end

Dates.DateTime(t::TAITime) =
    J2000 + Millisecond(1000 * t.second + round(Int64, 1000 * t.fraction))

Base.:-(a::TAITime, b::TAITime) = (a.second - b.second) + (a.fraction - b.fraction)

function Base.:+(t::TAITime, Δ::Real)
    whole = floor(Int64, Δ)
    TAITime(t.second + whole, t.fraction + (Δ - whole))
end
Base.:-(t::TAITime, Δ::Real) = t + (-Δ)

Base.isless(a::TAITime, b::TAITime) = isless((a.second, a.fraction), (b.second, b.fraction))

# As AstroTime's `isapprox` for epochs: on the seconds since J2000 as one `Float64`.
Base.isapprox(a::TAITime, b::TAITime; kwargs...) =
    isapprox(a.second + a.fraction, b.second + b.fraction; kwargs...)

Base.show(io::IO, t::TAITime) = print(io, DateTime(t), " TAI")

end
