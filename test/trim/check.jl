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

# Field by field: text exactly, numbers to within 1e-12 relative or 1e-9 absolute (metres,
# m/s and seconds). The two are not bit for bit the same everywhere: on a GitHub runner
# the velocities differed in their last two or three digits, while the same Julia
# version matched exactly on a local machine.
function same_field(a, b)
    a == b && return true
    x, y = tryparse(Float64, a), tryparse(Float64, b)
    x !== nothing && y !== nothing && isapprox(x, y; rtol = 1e-12, atol = 1e-9)
end
function same_line(a, b)
    fields_a, fields_b = split(a, (' ', '=')), split(b, (' ', '='))
    length(fields_a) == length(fields_b) && all(splat(same_field), zip(fields_a, fields_b))
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
