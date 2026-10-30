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
    # Public qualified bindings also appear in names on Julia 1.11+.
    exported = filter(
        name -> Base.isexported(PVT, name) && name !== :PositionVelocityTime,
        names(PVT),
    )
    @test sort(bare) == sort(exported)
    @test all(name -> isdefined(PVT, name) && name ∉ exported, qualified)
    @static if VERSION >= v"1.11"
        @test all(name -> Base.ispublic(PVT, name), qualified)
        @test !Base.ispublic(PVT, :empty_keeping_capacity!)
        @test !Base.ispublic(PVT, :check_observation_keys)
    end
end
