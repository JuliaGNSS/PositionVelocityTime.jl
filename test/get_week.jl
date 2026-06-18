using PositionVelocityTime: get_week
using GNSSDecoder: GNSSDecoderState, GPSL1Data
using GNSSSignals: GPSL1CA, GalileoE1B

# Build a minimal GNSSDecoderState{GPSL1Data} with the given broadcast
# 10-bit `trans_week`. `get_week` only reads `decoder.data.trans_week`
# so the rest of the state can stay at its default zero values.
function decoder_with_trans_week(trans_week)
    base = GNSSDecoderState(GPSL1CA(), 1)
    GNSSDecoderState(base; data = GPSL1Data(base.data; trans_week))
end

@testset "get_week resolves GPS L1 week-rollover ambiguity" begin
    # GPS week 0 begins 1980-01-06.
    # Cycle boundaries:  WN=1024 → 1999-08-22, WN=2048 → 2019-04-07,
    #                    WN=3072 → 2038-11-21, WN=4096 → 2058-07-08.

    # Cycle 0 (1980-01-06 .. 1999-08-22): broadcast WN ∈ [0, 1023]
    # ION RTL-SDR sample data was recorded 2017-09-10 in cycle 1,
    # broadcast trans_week ≈ 941 → absolute week 1024 + 941 = 1965.
    @testset "anchor in 2017 picks cycle 1 (post-1999 rollover)" begin
        decoder = decoder_with_trans_week(941)
        @test get_week(decoder; approximate_year = 2017) == 1024 + 941
        # ±9 years of slop still picks the correct cycle.
        @test get_week(decoder; approximate_year = 2009) == 1024 + 941
        @test get_week(decoder; approximate_year = 2018) == 1024 + 941
    end

    # Same trans_week but anchor in cycle 2 (post-2019 rollover):
    # absolute week 2048 + 941 = 2989 → ~April 2037. The picked cycle is
    # whichever places the absolute date closer to the anchor; with
    # trans_week = 941 (~mid-cycle), cycle 1 covers 2017 and cycle 2
    # covers 2037, so anchors near 2037 pick cycle 2 and anchors before
    # ~2027 still pick cycle 1.
    @testset "anchor in 2037 picks cycle 2 (post-2019 rollover)" begin
        decoder = decoder_with_trans_week(941)
        @test get_week(decoder; approximate_year = 2037) == 2048 + 941
        @test get_week(decoder; approximate_year = 2030) == 2048 + 941  # closer to 2037 than 2017
    end

    # Anchor in cycle 0 (between epoch and first rollover).
    @testset "anchor in 1995 picks cycle 0" begin
        decoder = decoder_with_trans_week(500)
        @test get_week(decoder; approximate_year = 1995) == 0 + 500
    end

    # Anchor in cycle 3 (after the upcoming 2038 rollover).
    @testset "anchor in 2050 picks cycle 3" begin
        decoder = decoder_with_trans_week(100)
        @test get_week(decoder; approximate_year = 2050) == 3072 + 100
    end

    # Boundary case: trans_week = 0 with anchor near the rollover.
    # 2019-04-07 (rollover to cycle 2) — anchor at 2019 should pick the
    # cycle whose week 0 is closest. Mid-2019 is closer to start-of-2019
    # (cycle 1, end) than start-of-2020 (cycle 2, beginning), so cycle 1
    # is selected.
    @testset "trans_week = 0 near rollover" begin
        decoder = decoder_with_trans_week(0)
        # Mid-2019 anchor ≈ week 2050, distance to cycle 1 base (week 1024)
        # is 1026 weeks; distance to cycle 2 base (week 2048) is 2 weeks.
        # So cycle 2 wins.
        @test get_week(decoder; approximate_year = 2019) == 2048
    end

    # Sanity-check that the default anchor matches an explicit current-UTC-year
    # anchor — both should pick the same cycle (the default reads `now()` at
    # call time, and the explicit comparison value reads it on the same line,
    # so they agree by construction).
    @testset "default anchor uses current UTC year" begin
        decoder = decoder_with_trans_week(500)
        current_year = Dates.year(Dates.now(Dates.UTC))
        @test get_week(decoder) == get_week(decoder; approximate_year = current_year)
    end
end
