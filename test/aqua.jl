using Aqua

@testset "Aqua.jl" begin
    # `persistent_tasks` resolves the package from the registry into a fresh
    # temporary project, so it cannot pass while this branch requires the
    # not-yet-released GNSSDecoder 4 — re-enable with that release.
    Aqua.test_all(PositionVelocityTime; persistent_tasks = false)
end
