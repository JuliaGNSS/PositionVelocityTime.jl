# Opt-in real-data PVT integration test.
#
# Computes a full Position/Velocity/Time fix end to end from a real-sky GNSS
# recording: Acquisition.jl acquires GPS L1 C/A satellites, Tracking.jl tracks
# them, GNSSDecoder.jl decodes the navigation message, and `calc_pvt` turns the
# resulting `SatelliteState`s into a position — which we check lands at the
# recording site. This is the real-signal counterpart to the synthetic
# `pvt.jl` / `pvt_iono_tropo.jl` unit tests, exercising the whole JuliaGNSS
# stack the way a receiver actually uses it.
#
# WHY THIS CAPTURE (and not the Tracking.jl / GNSSDecoder.jl "III-7a" one):
#   Tracking.jl and GNSSDecoder.jl share the ~15.7 s Hanoi "III-7a_short"
#   capture. That is long enough to validate a single parity-checked subframe
#   (GNSSDecoder asserts only on `raw_data.TOW`), but NOT long enough to decode
#   a complete GPS ephemeris — subframes 1+2+3 with a consistent IODC span
#   ≥18 s. `calc_pvt` needs that full ephemeris (`decoder.data`) plus a healthy
#   flag for ≥4 satellites to place them in orbit, so a real fix is impossible
#   from III-7a. We therefore use the Fraunhofer IIS "L125" static roof
#   recording (Nuremberg, 2014-09-23, 210 s, CC BY-NC), of which we use the
#   first ~40 s — comfortably past the ephemeris span with margin for
#   acquisition + tracking lock-in.
#
# WHY ONLY GPS L1 C/A:
#   `calc_pvt` needs ≥4 healthy satellites of the *same* GNSS. The open-sky
#   roof recording has plenty of GPS L1 C/A satellites; the other systems/bands
#   present in the wider L125 set have fewer than 4 in view, so GPS L1 C/A is
#   the only constellation that can yield a fix here.
#
# CAPTURE FORMAT (`*_L1.bin`, the de-multiplexed L1 band of Flexiband config
# III-1b — parameters taken from Fraunhofer's own GNSS-SDR config
# `conf/gnss-sdr_Hybrid_byte.conf`, which points directly at this file):
#   • zero-IF complex baseband (`InputFilter.IF=0`), L1 RF centre 1575.42 MHz
#   • fs = 20 MHz (`SignalSource.sampling_frequency=20000000`)
#   • interleaved signed 8-bit I/Q, two bytes per complex sample
#     (`item_type=byte` + `Ibyte_To_Complex`): byte 2k-1 = I, byte 2k = Q.
#   The `.bin` is wrapped in a single-entry ZIP64 archive, DEFLATE-compressed.
#
# DOWNLOAD + CACHE + STREAMING DEMUX (we never materialise the decoded signal):
#   The full archive is ~3.45 GB (≈8.4 GB uncompressed at 2 B/sample). We
#   range-download only a bounded compressed *prefix* of the zip — enough to
#   inflate past `L125_SECONDS` — and cache that prefix as-is (set
#   PVT_TESTDATA_DIR to the cache dir, so CI skips the download on a hit). The
#   prefix stays compressed on disk; the test never writes the multi-GB decoded
#   signal. Instead, mirroring Tracking.jl's streaming demux
#   (JuliaGNSS/Tracking.jl#158), it raw-inflates the prefix and demuxes one
#   ~0.2 s window of int8 I/Q into ComplexF32 at a time, feeding each window
#   straight into `track!`. Peak memory is a single window (~16 MB), not the
#   ~16 GB a whole-span decode would need.

if get(ENV, "PVT_RUN_INTEGRATION_TEST", "false") != "true"
    @info "Skipping L125 real-data PVT integration test (downloads a large capture). " *
          "Set ENV[\"PVT_RUN_INTEGRATION_TEST\"] = \"true\" to run it."
else
    using Downloads: Downloads
    using CodecZlib: DeflateDecompressorStream
    using Unitful: Hz, ustrip
    using GNSSSignals: GPSL1CA
    using GNSSDecoder: GPSL1CADecoderState, decode, is_sat_healthy
    using Acquisition: acquire, is_detected
    using Tracking: TrackState, add_satellite!, track!, get_soft_bits, get_sat_state
    using Geodesy: LLA, ECEF, ECEFfromLLA, wgs84, euclidean_distance

    const L125_URL =
        "https://www2.iis.fraunhofer.de/flexiband/reference-data/" *
        "20140923_20-24-17_L125_roof_210s_L1.bin.zip"
    const L125_ZIP_SIZE = 3_453_820_212        # exact archive byte count
    const L125_FS = 20.0e6Hz                    # L1 band sampling rate
    const L125_IF = 0.0Hz                       # zero-IF complex baseband
    const L125_SECONDS = 100                    # how much signal we extract + track

    # 100 s spans ~3.3 GPS LNAV frames (a frame is 30 s, and a complete ephemeris
    # is subframes 1+2+3). From a cold acquisition that random-phase start means
    # a single satellite can need up to ~45 s to deliver a consecutive, IODC-
    # consistent 1-2-3 set, so 40 s yields only a handful of complete ephemerides;
    # 100 s clears the ≥4 healthy satellites `calc_pvt` needs with margin to spare
    # (so the fix is well-conditioned rather than balanced on the 4-satellite edge).

    # Bounded compressed prefix to fetch (and cache). The archive compresses at a
    # measured ~0.41 ratio (≈16 MB of zip per second of signal), so 100 s of
    # signal is ~1.6 GB compressed; 2.0 GB inflates to ~123 s, a safe margin.
    const L125_PREFIX_BYTES = min(2_000_000_000, L125_ZIP_SIZE)

    # Recording site: Fraunhofer IIS roof in north-east Nuremberg. The exact
    # surveyed antenna position is not published, so the reference is the position
    # this pipeline decodes. The fix is fully deterministic — identical to every
    # printed digit whether tracking runs single- or multi-threaded — and the
    # satellite set is pinned (`L125_PRNS`), so the geometry is fixed too. The
    # 100 m tolerance therefore only has to absorb numeric drift from future
    # dependency patch updates, while still pinning the solve to this rooftop (a
    # broken solve lands kilometres away or non-finite).
    const L125_REF_LLA = LLA(49.4872661, 11.1253752, 421.07)
    const L125_REF_TOL_M = 100.0
    const L125_GPS_WEEK = 1811                  # GPS week of 2014-09-23

    # Satellites that decode a complete, healthy ephemeris within L125_SECONDS of
    # this capture; the other acquired PRNs are too weak to yield a clean subframe
    # 1-2-3 set. Pinning the solve to them keeps the geometry — and the fix —
    # reproducible, which is what lets us assert the position to 100 m.
    const L125_PRNS = [13, 17, 20, 23]

    # Ensure the bounded compressed prefix of the archive is present in the
    # test-data directory (cached on CI), returning its path. On a miss this
    # range-downloads `L125_PREFIX_BYTES` of the zip in chunks, with retries
    # since the server can drop long connections. The prefix is kept compressed —
    # it is never expanded to a multi-GB file on disk.
    function l125_fetch()
        dir = get(ENV, "PVT_TESTDATA_DIR", tempdir())
        mkpath(dir)
        path = joinpath(dir, "L125_roof_210s_L1.prefix.zip")
        (isfile(path) && filesize(path) >= L125_PREFIX_BYTES) && return path
        tmp = path * ".part"
        open(tmp, "w") do out
            start = 0
            while start < L125_PREFIX_BYTES
                stop = min(start + 200_000_000 - 1, L125_PREFIX_BYTES - 1)
                want = stop - start + 1
                tries = 0
                while true
                    (tries += 1) > 30 && error("L125 download: chunk at $start failed")
                    part = tempname()
                    try
                        Downloads.download(
                            L125_URL,
                            part;
                            headers = ["Range" => "bytes=$start-$stop"],
                        )
                        if filesize(part) == want
                            write(out, read(part))
                            rm(part; force = true)
                            break
                        end
                    catch err
                        @warn "L125 chunk retry" start tries err
                    end
                    rm(part; force = true)
                end
                start = stop + 1
            end
        end
        mv(tmp, path; force = true)
        path
    end

    # Open the cached prefix and return a raw-DEFLATE stream positioned at the
    # start of the `.bin` payload, by parsing the ZIP local file header (30-byte
    # fixed part + filename + extra field; method 8 = DEFLATE). Reading from the
    # returned stream inflates on demand, so each window is decompressed only as
    # it is consumed.
    function l125_open(path)
        io = open(path)
        read(io, UInt32) == 0x04034b50 || error("L125: not a ZIP local file header")
        skip(io, 4)                       # version-needed + general-purpose flags
        read(io, UInt16) == 0x0008 || error("L125: entry is not DEFLATE-compressed")
        skip(io, 16)                      # mod time/date + CRC-32 + comp/uncomp sizes
        namelen = read(io, UInt16)
        extralen = read(io, UInt16)
        skip(io, Int(namelen) + Int(extralen))
        DeflateDecompressorStream(io)
    end

    # Inflate and demux the next `nsamples` complex samples from the stream
    # (signed 8-bit I/Q, two bytes each: byte 2k-1 = I, byte 2k = Q). Returns
    # fewer than requested once the prefix is exhausted.
    function l125_read_window(stream, nsamples)
        buf = read(stream, 2 * nsamples)
        n = length(buf) - (length(buf) % 2)
        out = Vector{ComplexF32}(undef, n ÷ 2)
        @inbounds for k in eachindex(out)
            out[k] = ComplexF32(reinterpret(Int8, buf[2k-1]), reinterpret(Int8, buf[2k]))
        end
        out
    end

    @testset "Flexiband L125 real-data PVT fix (GPS L1 C/A, Nuremberg roof)" begin
        path = l125_fetch()
        fs = ustrip(Hz, L125_FS)

        # Acquire GPS L1 C/A on the first 40 ms (10 ms coherent-x-noncoherent is
        # plenty for the strong open-sky roof signals), inflated + demuxed from
        # just that leading window.
        stream = l125_open(path)
        acq_block = l125_read_window(stream, round(Int, fs * 0.040))
        close(stream)
        acq = filter(
            r -> is_detected(r; pfa = 1e-8),
            acquire(
                GPSL1CA(),
                acq_block,
                L125_FS,
                1:32;
                interm_freq = L125_IF,
                num_coherently_integrated_code_periods = 1,
                num_noncoherent_accumulations = 10,
            ),
        )
        prns = sort!([a.prn for a in acq])
        @info "Acquired GPS L1 C/A satellites" prns
        @test length(prns) >= 4

        # Re-open the prefix and stream it from the start (so the acquisition code
        # phase and Doppler line up with sample 0), inflating + demuxing one 0.2 s
        # window at a time and feeding each straight into `track!`. Only a single
        # window of samples is ever live, so peak memory stays flat regardless of
        # `L125_SECONDS`.
        ts = TrackState(; signals = (gps = (GPSL1CA(),),))
        for a in acq
            ts = add_satellite!(ts, a; group = :gps)
        end
        decoders = Dict(p => GPSL1CADecoderState(p) for p in prns)

        stream = l125_open(path)
        window_len = round(Int, fs * 0.2)
        total = round(Int, fs * L125_SECONDS)
        processed = 0
        while processed < total
            window = l125_read_window(stream, min(window_len, total - processed))
            isempty(window) && break
            processed += length(window)
            ts = track!(window, ts, L125_FS; intermediate_frequency = L125_IF)
            for p in prns
                s = get_soft_bits(ts, :gps, p)
                isempty(s) || (decoders[p] = decode(decoders[p], s, length(s)))
            end
        end
        close(stream)
        @test processed >= round(Int, fs * (L125_SECONDS - 5))   # streamed the prefix

        # Pin the solve (via the Tracking extension) to the satellites that
        # reliably decode a complete, healthy ephemeris in this capture, so the
        # geometry — and the fix — is reproducible and can be checked tightly.
        healthy = [p for p in prns if is_sat_healthy(decoders[p])]
        @info "Satellites with a decoded healthy ephemeris" healthy
        @test issubset(L125_PRNS, healthy)
        states = [
            SatelliteState(decoders[p], GPSL1CA(), get_sat_state(ts, :gps, p)) for
            p in L125_PRNS
        ]

        # Solve PVT. `approximate_year` resolves the GPS 1024-week rollover; the
        # recording is from 2014, not "now".
        pvt = calc_pvt(states; approximate_year = 2014)
        @test get_num_used_sats(pvt) >= 4

        lla = get_LLA(pvt)
        @info "PVT fix" lat = lla.lat lon = lla.lon alt = lla.alt used =
            get_num_used_sats(pvt) gdop = get_gdop(pvt)
        ref_ecef = ECEFfromLLA(wgs84)(L125_REF_LLA)
        @test euclidean_distance(pvt.position, ref_ecef) < L125_REF_TOL_M
        @test isfinite(get_gdop(pvt)) && get_gdop(pvt) > 0

        # The decoded navigation time must place the recording on 2014-09-23
        # (GPS week 1811) — a genuine, parity-/ephemeris-validated end-to-end
        # decode, not just a plausible-looking geometry.
        @test all(
            s -> PositionVelocityTime.get_week(s.decoder; approximate_year = 2014) ==
                L125_GPS_WEEK,
            states,
        )
        @test !isnothing(pvt.time)
    end
end
