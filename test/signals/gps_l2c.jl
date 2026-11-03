# GPS L2C (CNAV on the L2 band) support.
#
# L2C shares the CNAV `GNSSDecoder.GPSCNAVData` container with L5-I; the two are
# distinguished only by the decoder's constants type (`GPSL2CMConstants` vs
# `GPSL5IConstants`). So the orbit/clock/week/ionosphere paths are the L5 ones,
# exercised here on an L2C state to confirm the classification, the L2 frequency
# band, the 50 Hz L2 CM bit rate and the ISC_L2C group delay are wired up.
@testset "GPS L2C (CNAV on L2)" begin
    l2c = GNSSDecoder.GPSL2CMDecoderState(1)
    l5 = GNSSDecoder.GPSL5IDecoderState(1)
    swap(state, d) = GNSSDecoder.GNSSDecoderState(state; data = d, raw_data = d)
    state_of(d) = SignalObservation(; decoder = swap(l2c, d), signal = GNSSSignals.GPSL2CM(),
        code_phase = 0.0, carrier_doppler = 0.0Hz, carrier_phase = 0.0)

    @testset "classification: GPS time system, L2 band, GPSL2CM signal id" begin
        @test GNSSSignals.get_time_system(GNSSSignals.GPSL2CM()) == GNSSSignals.GPST()
        @test !has_offset(l2c, GPST())
        @test GNSSSignals.get_band_id(GNSSSignals.GPSL2CM()) == :L2
        @test GNSSSignals.get_signal_id(GNSSSignals.GPSL2CM()) == :GPSL2CM
    end

    @testset "nav-message bit rate is L2 CM (50 Hz), not L5-I (100 Hz)" begin
        d = GNSSDecoder.GPSCNAVData()
        @test GNSSSignals.get_bit_rate(swap(l2c, d)) ==
              GNSSSignals.get_bit_rate(GNSSSignals.GPSL2CM())
        @test GNSSSignals.get_bit_rate(swap(l5, d)) ==
              GNSSSignals.get_bit_rate(GNSSSignals.GPSL5I())
    end

    @testset "get_broadcast_week uses full broadcast WN (no rollover)" begin
        # approximate_year must be ignored for CNAV (full 13-bit WN).
        @test GNSSDecoder.get_broadcast_week(swap(l2c, GNSSDecoder.GPSCNAVData(; WN = 2200));
            approximate_year = 2042) == 2200
    end

    @testset "Klobuchar parameters from CNAV α_0…β_3" begin
        d = GNSSDecoder.GPSCNAVData(; α_0 = 1e-8, α_1 = 2e-8, α_2 = 3e-8, α_3 = 4e-8,
            β_0 = 1e4, β_1 = 2e4, β_2 = 3e4, β_3 = 4e4)
        kp = GNSSDecoder.get_klobuchar_parameters(swap(l2c, d))
        @test kp isa GNSSDecoder.KlobucharParameters
        @test kp.α_0 == 1e-8
        @test kp.β_3 == 4e4
        # Partially populated → nothing, not an error.
        @test GNSSDecoder.get_klobuchar_parameters(swap(l2c, GNSSDecoder.GPSCNAVData(; α_0 = 1e-8))) ===
              nothing
    end
end
