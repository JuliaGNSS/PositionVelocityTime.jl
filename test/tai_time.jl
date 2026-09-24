@testset "TAITime" begin
    # The J2000 origin and a fraction normalised into [0, 1).
    @test TAITime(DateTime(2000, 1, 1, 12)) == TAITime(0, 0.0)
    @test TAITime(10, 1.25) == TAITime(11, 0.25)
    @test TAITime(10, -0.25) == TAITime(9, 0.75)

    # A calendar label round-trips at millisecond resolution, before J2000 as well.
    for datetime in (DateTime(2021, 5, 31, 12, 53, 14, 118), DateTime(1980, 1, 6, 0, 0, 19))
        @test DateTime(TAITime(datetime)) == datetime
    end

    t = TAITime(675737594, 0.75)
    @test t + 0.5 == TAITime(675737595, 0.25)
    @test t - 1.0 == TAITime(675737593, 0.75)
    @test (t + 1234.5) - t == 1234.5
    @test t < t + 1e-6
    @test t ≈ t + 1e-9
    @test sprint(show, TAITime(DateTime(2021, 5, 31, 12, 53, 14, 118))) ==
          "2021-05-31T12:53:14.118 TAI"

    @testset "AstroTime extension" begin
        epoch = TAIEpoch(2021, 5, 31, 12, 53, 14.1183385390904732)
        @test TAIEpoch(TAITime(epoch)) == epoch
        @test TAITime(epoch) == TAITime(epoch.second, epoch.fraction)
        # The same calendar label on both sides, so the two share their origin.
        @test TAITime(DateTime(2021, 5, 31, 12, 53, 14)) ==
              TAITime(TAIEpoch(2021, 5, 31, 12, 53, 14.0))
    end
end
