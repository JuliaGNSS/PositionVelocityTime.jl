# Builds `pvt_app.jl` into a `juliac --trim=safe` executable and checks that it prints
# what the same app prints in a regular Julia session. Needs Julia 1.12+ and
# the JuliaC app (`pkg> app add JuliaC`):
#
#     julia --project=test/trim -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/trim test/trim/check.jl
#
# The trim verifier rejects any call it cannot resolve statically, so a build failure
# here names the dynamic dispatch that crept into the solve.
juliac = something(
    Sys.which("juliac"),
    joinpath(first(DEPOT_PATH), "bin", Sys.iswindows() ? "juliac.bat" : "juliac"),
)
isfile(juliac) || error("juliac not found; install it with `pkg> app add JuliaC`")

app = joinpath(@__DIR__, "pvt_app.jl")
build = mktempdir()
run(`$juliac --output-exe pvt_app --project $(@__DIR__) --bundle $build --trim=safe
    --experimental $app`)

executable = joinpath(build, "bin", Sys.iswindows() ? "pvt_app.exe" : "pvt_app")
trimmed = read(`$executable`, String)
regular = read(`$(Base.julia_cmd()) --project=$(@__DIR__) $app`, String)

# Field by field: text exactly, numbers to within 1e-12 relative or a per-quantity absolute
# tolerance. The two are not bit for bit the same everywhere: on a GitHub runner the
# velocities differed in their last two or three digits, while the same Julia version
# matched exactly on a local machine. Such last-bit differences can also decide whether the
# position solve takes one more sub-millimetre step: two warm starts of the same epoch
# already end up to ~3e-5 m apart in one session, so the regular session and the trimmed
# executable may land as far apart too. Positions and the inter-system biases (m) are
# therefore compared to 1e-4 m, and velocities, which follow them through the range-rate
# model, to 1e-7 m/s; times, DOPs and counts stay at 1e-9. A trimming error shows up as
# a wrong value, not as a difference of this size.
const POSITION_ATOL = 1e-4  # m
const VELOCITY_ATOL = 1e-7  # m/s
const DEFAULT_ATOL = 1e-9

function same_number(a, b, atol)
    a == b && return true
    x, y = tryparse(Float64, a), tryparse(Float64, b)
    x !== nothing && y !== nothing && isapprox(x, y; rtol = 1e-12, atol)
end

# A line of `pvt_app.jl` is `name: x y z vx vy vz second fraction GDOP satellites`,
# followed by one `system=bias` per inter-system bias.
value_atol(i) = i <= 3 ? POSITION_ATOL : i <= 6 ? VELOCITY_ATOL : DEFAULT_ATOL
function same_value(i, a, b)
    occursin('=', a) || return same_number(a, b, value_atol(i))
    system_a, bias_a = split(a, '='; limit = 2)
    system_b, bias_b = split(b, '='; limit = 2)
    system_a == system_b && same_number(bias_a, bias_b, POSITION_ATOL)
end
function same_line(a, b)
    a == b && return true
    parts_a, parts_b = split(a, ": "; limit = 2), split(b, ": "; limit = 2)
    length(parts_a) == length(parts_b) == 2 && parts_a[1] == parts_b[1] || return false
    values_a, values_b = split(parts_a[2], ' '), split(parts_b[2], ' ')
    length(values_a) == length(values_b) &&
        all(i -> same_value(i, values_a[i], values_b[i]), eachindex(values_a))
end
function same_output(a, b)
    lines_a, lines_b = split(a, '\n'), split(b, '\n')
    length(lines_a) == length(lines_b) && all(splat(same_line), zip(lines_a, lines_b))
end

print(trimmed)
if !same_output(trimmed, regular)
    println("\nThe regular session printed instead:\n", regular)
    error("the trimmed executable's output differs from the regular session's")
end
println("Trimmed executable ($(filesize(executable) ÷ 1024) KiB) matches the regular ",
    "session on all $(count(==('\n'), trimmed)) solves.")
