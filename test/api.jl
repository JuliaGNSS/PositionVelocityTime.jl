# The API is what the reference documents: its bare names are the exports, its qualified
# ones unexported names bound by hand. The internals documented for reference, after
# "The unexported internals", and the trim-safe stand-ins are not part of it.
@testset "the API reference lists the exports bare and the rest qualified" begin
    PVT = PositionVelocityTime
    bare, qualified = Symbol[], Symbol[]
    internals, in_docs = false, false
    for line in eachline(joinpath(@__DIR__, "..", "docs", "src", "api.md"))
        startswith(line, "## ") && (internals = line == "## Trimmed Executables")
        startswith(line, "The unexported internals") && (internals = true)
        line == "```@docs" && (in_docs = true; continue)
        line == "```" && (in_docs = false; continue)
        (in_docs && !internals && !startswith(line, "Base.")) || continue
        if startswith(line, "PositionVelocityTime.")
            push!(qualified, Symbol(chopprefix(line, "PositionVelocityTime.")))
        else
            push!(bare, Symbol(line))
        end
    end
    # Without `all`, `names` lists the exports.
    exported = setdiff(names(PVT), [:PositionVelocityTime])
    @test sort(bare) == sort(exported)
    @test all(name -> isdefined(PVT, name) && name ∉ exported, qualified)
end
