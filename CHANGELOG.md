# Changelog

# [7.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v6.0.1...v7.0.0) (2026-09-29)

## [6.0.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v6.0.0...v6.0.1) (2026-09-29)


### Bug Fixes

* accept GNSSDecoder 5 ([a0f37f9](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/a0f37f9247e1edf8a50eb0940cc09c49e1c66416))

# [6.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.4.0...v6.0.0) (2026-09-28)


* feat!: compile the PVT solve into a trimmed juliac executable ([6b820d6](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/6b820d669dfd053e080f65e0af3756cee199e0eb))
* feat!: solve the PVT without allocating, into an explicit output ([f1dd2ab](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/f1dd2ab9b7b21bcda0f6d3ae5576b1da8a2ea80d))
* feat!: take signal groups, and solve behind a flat-row barrier ([c06b5d7](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/c06b5d77957acdb111c1841fd715e75fa9b08b2e))


### Bug Fixes

* refuse a duplicate (signal, PRN) pair before writing the solution ([506857d](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/506857deab1041116864f12cf3031eb33cd2971d))
* refuse groups of bare state vectors with the migration directions ([7cc0811](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/7cc0811c9a35e0c81284b29edaacad90db3bb0f5))


### Performance Improvements

* count distinct satellites without boxing a key per satellite ([49fb434](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/49fb4347ad1eac2fcee00902fbe61b654e4a919a))


### BREAKING CHANGES

* `PVTSolution` is a mutable struct.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
* `PVTSolution.time` is a `TAITime` rather than an AstroTime
`TAIEpoch`, and AstroTime is no longer loaded with the package —
`using AstroTime` and `TAIEpoch(pvt.time)` converts it. `PVTSolution`'s
`position` and `velocity` and `SatInfo.position` are `ECEF{Float64}`, and
`inter_system_biases` is keyed (and `reference_system` typed) by
`SupportedTimeSystem` rather than `GNSSSignals.TimeSystem`.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
* `calc_pvt` takes signal groups, not a vector of
`SatelliteState`s. `calc_pvt(PositionVelocityTime.signal_groups(states))` is
the mechanical translation and the documented bridge, though it is
inference-blind by construction; a receiver that already keeps its satellites
per signal should name its groups instead. With Tracking loaded,
`signal_groups(track_state, decoders)` converts a whole `TrackState`.

The documented Measurement-Model Surface is re-specified on the flat row:
`decide_bias_layout`, `predict_atmospheric_delays`, `calc_hub_range_offsets`,
`calc_user_velocity_and_clock_drift`, `calc_time_scale_offsets` and
`select_ionospheric_correction` take measurements (or groups) rather than
satellite states plus parallel classification vectors; `ionospheric_delay`
takes a carrier frequency rather than a signal; `day_of_year` takes a system
start date; and `calc_steering_offset` takes a `BroadcastTimeOffset`
(`broadcast_time_offset(decoder, target)` builds one).

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>

# [5.4.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.3.0...v5.4.0) (2026-09-05)


### Bug Fixes

* fold the week wrap out of the pseudorange differencing ([dac3b5f](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/dac3b5f88d72978823bd33215a264647c56e1cdc)), closes [GNSSDecoder.jl#89](https://github.com/GNSSDecoder.jl/issues/89)
* survive week rollovers, spurious warm-start roots, and iono order ([9e48961](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/9e489611b3631724402bc8c18089efa236790ceb))


### Features

* collapse scarce epochs onto any hub system, not just GPS Time ([edc8464](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/edc8464b4ea9a1ee437d96657fcad5ad071b8ebd))
* implement BDGIM, the BDS-3 broadcast ionospheric model ([a77ec0b](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/a77ec0bf991f46b8d9d71b3b9916b26bfba0023c))
* support the Galileo E5b/E6-B and BeiDou signals ([022b051](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/022b0511c0be4748046fc88a4fce76e2dd99af18))

# [5.3.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.2.0...v5.3.0) (2026-09-03)


### Features

* precompile every navigation-data type the solver dispatches on ([d44c8f1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/d44c8f18f6ab11492f389cd829f466f7f21aa625))
* precompile the Galileo and mixed-constellation solves as well ([851a0b2](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/851a0b21ac404427471731adef98cbe7649d6ada))
* precompile the PVT solve ([5d06f1c](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/5d06f1c1e223e02834cf98bc7c6125a448ff9308)), closes [GNSSReceiver.jl#107](https://github.com/GNSSReceiver.jl/issues/107)

# [5.2.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.1.0...v5.2.0) (2026-08-20)


### Features

* support Tracking 8 ([9fc9e4f](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/9fc9e4f1623fcb3f5258835dee78d20eb7492683)), closes [JuliaGNSS/Tracking.jl#229](https://github.com/JuliaGNSS/Tracking.jl/issues/229)

# [5.1.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.0.3...v5.1.0) (2026-08-18)


### Features

* support Tracking 7 ([91069f6](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/91069f68fbe5e3f5b16216c1f3adce32c685641c)), closes [JuliaGNSS/Tracking.jl#223](https://github.com/JuliaGNSS/Tracking.jl/issues/223)

## [5.0.3](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.0.2...v5.0.3) (2026-08-06)


### Bug Fixes

* map the tropospheric delay with the Niell mapping functions ([4153e1c](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/4153e1cddbeba9d8439ff5d638d6a862c59d787c)), closes [#62](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/62)

## [5.0.2](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.0.1...v5.0.2) (2026-08-05)


### Bug Fixes

* report residuals as observed minus computed ([1bb2fe8](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/1bb2fe8fe13e478f82f094abe6cd6b00cfb69646))

## [5.0.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v5.0.0...v5.0.1) (2026-08-05)


### Bug Fixes

* run the position solve to its optimum ([cc487ce](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/cc487ce57641703e01f57faef428bebb02d7f3f1))

# [5.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.3.0...v5.0.0) (2026-08-05)


* feat!: report post-fit range-rate residuals per satellite ([65aba95](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/65aba951deb6269289046f73c8e921f0bd5a25ab))


### BREAKING CHANGES

* `SatInfo` gains a `rate_residual` field, so the exported struct's
positional-constructor arity/order and field layout change. Keyword construction
and field access are unaffected; positional construction and reflection-based code
must be updated.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>

# [4.3.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.2.2...v4.3.0) (2026-08-05)


### Features

* support Tracking.jl v6 ([79f8a2a](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/79f8a2ab52f5a5655d62bed3752a71920c96656c))

## [4.2.2](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.2.1...v4.2.2) (2026-08-04)


### Bug Fixes

* correct clock drift epoch, velocity singularity and carrier-phase unit ([81f785d](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/81f785da09a99eb99514492515c7a78ba7748dd4)), closes [#41](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/41) [#42](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/42) [#59](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/59) [#41](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/41) [#42](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/42) [#59](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/59)

## [4.2.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.2.0...v4.2.1) (2026-08-04)


### Bug Fixes

* bound satellite elevation in tropospheric model ([678c429](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/678c4294968eac2bf75498a652e907b78d6b52f2))

# [4.2.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.1.2...v4.2.0) (2026-08-03)


### Features

* support Tracking 5 ([010ffa3](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/010ffa3cd097f1ec52d11dea6e82162008eae3a4))

## [4.1.2](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.1.1...v4.1.2) (2026-07-30)


### Bug Fixes

* skip unsolvable PVT epochs instead of throwing ([8af96c5](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/8af96c516ff0bc4b339ec91a2c68d7b8e6c09626))

## [4.1.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.1.0...v4.1.1) (2026-07-29)


### Bug Fixes

* reject degenerate PVT geometries instead of throwing ([c83ff9a](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/c83ff9a4a44af818d496127901123981533d29a8))

# [4.1.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.0.1...v4.1.0) (2026-07-26)


### Features

* support Tracking 4 ([d35afec](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/d35afec))

## [4.0.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v4.0.0...v4.0.1) (2026-07-07)


### Bug Fixes

* gate PVT satellites on decoding completeness ([8fe32d4](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/8fe32d4f29449896e916877be910901b212a659c))

# [4.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v3.1.0...v4.0.0) (2026-07-07)


* feat!: enrich PVT solution with IFB reference bands, units, and course over ground ([6a528ce](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/6a528ce93978a5c6ab1b829ba3eec364188e1d70))


### BREAKING CHANGES

* `PVTSolution` gains a `course_over_ground` field, so the exported
struct's positional-constructor arity/order and field layout change. Several fields
are now Unitful quantities rather than bare `Float64` metres: `time_correction`, the
`inter_system_biases` values and `SatInfo.residual` are `typeof(1.0m)`. And
`inter_frequency_biases` is now `Dict{Symbol,InterFrequencyBias}` — read `.value`
(a `typeof(1.0m)`) and `.reference` (the anchor band) instead of a bare `Float64`.
Consumers doing arithmetic must handle the units (e.g. `ustrip(u"m", x)`), and
positional construction and reflection-based code must be updated. Reading via
keyword construction and the accessor fields is otherwise unaffected.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>

# [3.1.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v3.0.0...v3.1.0) (2026-07-06)


### Features

* support Tracking 3 ([984d893](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/984d893))

# [3.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v2.2.0...v3.0.0) (2026-07-05)


* feat!: multi-GNSS PVT with GPS L2C/L5/L1C and Galileo E5a ([a968fff](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/a968ffff38c41132bef208abbe50a7047d79204f))


### BREAKING CHANGES

* `PVTSolution.sats` is now a
`Dictionary{Tuple{Symbol,Int},SatInfo}` (was `Dict{Int,SatInfo}`); index it
with `get_sat_info(pvt, signal, prn)`. `SatInfo` gains a `residual` field and
`calc_DOP` takes the user position and primary clock index. The
`get_gdop`/`get_pdop`/`get_hdop`/`get_vdop`/`get_tdop` accessors are removed;
read DOP from the `dop` field instead (e.g. `pvt.dop.GDOP`).
`get_num_used_sats` is removed: with `sats` keyed by (signal, PRN),
`length(pvt.sats)` counts measurements, not satellites.
`get_frequency_offset` is removed: compute it inline as
`pvt.relative_clock_drift * base_frequency`.
`PVTSolution.reference_system` and the `inter_system_biases` keys are now
`GNSSSignals.TimeSystem` values (GPST()/GST()), not :GPS/:Galileo symbols.
Requires GNSSSignals 3.3 and GNSSDecoder 3.6.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>

# [2.2.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v2.1.0...v2.2.0) (2026-06-23)


### Features

* support GNSSDecoder 2 ([da00512](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/da00512909b0219e2575f94133b94924a786e6cf))

# [2.1.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v2.0.0...v2.1.0) (2026-06-22)


### Features

* ionospheric and tropospheric corrections in calc_pvt ([#38](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/38)) ([c2a74e0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/c2a74e008d96fafedc247609cff8022119afd0eb))

# [2.0.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.6...v2.0.0) (2026-06-19)


* feat!: migrate to Tracking 2 (GNSSSignals 2, GNSSDecoder 1.3) ([b4bec52](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/b4bec520b5680e57b91aaed737356f6322e215b1))


### Bug Fixes

* **benchmark:** pick GPS L1 type by GNSSSignals version ([57e5520](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/57e552088bc572d8caca42728a043f3707b17738))
* **benchmark:** update fixtures for GNSSDecoder 1.3 / GNSSSignals 2 ([0387b36](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/0387b3627a3f63f603cfdbdbaebc2e523d27a5fc))


### BREAKING CHANGES

* drops support for GNSSSignals 1 and Tracking 1; the core
API now requires the v2 ecosystem.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>

## [1.0.6](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.5...v1.0.6) (2026-05-07)


### Bug Fixes

* resolve GPS L1 week-rollover ambiguity via approximate_year ([e7e2191](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/e7e219128346626da70432b618ca4a86ef9914c6))

## [1.0.5](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.4...v1.0.5) (2026-05-07)


### Performance Improvements

* avoid materializing healthy_states via findall + view ([263801d](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/263801d469c5f0e2207e46cdb989103899d53e00))

## [1.0.4](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.3...v1.0.4) (2026-05-07)


### Performance Improvements

* only apply geodesic acceleration on cold start (iszero prev_ξ) ([d73e3fa](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/d73e3fad9f2eed486c34fb0797198b6931181793))
* use geodesic acceleration in LM solve for user_position ([ac14723](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/ac14723e722a4ad66adcf5261bafc2c1fab94393))

## [1.0.3](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.2...v1.0.3) (2026-05-07)


### Performance Improvements

* use in-place LM model and Jacobian in user_position ([2c8fefe](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/2c8fefe1c432314ae8c3e9481dfc553694fbc195))

## [1.0.2](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.1...v1.0.2) (2026-05-07)


### Performance Improvements

* stack-allocate calc_DOP and reuse times in velocity solve ([84d968c](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/84d968c4d1060e96e85f4ceb6ecd62ca475023a5)), closes [#26](https://github.com/JuliaGNSS/PositionVelocityTime.jl/issues/26)

## [1.0.1](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v1.0.0...v1.0.1) (2026-05-07)


### Performance Improvements

* parameterize SatelliteState on decoder and system types ([1fc631f](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/1fc631f3cbefe886a4893cff6b29961408f13a10))

# [0.3.0](https://github.com/JuliaGNSS/PositionVelocityTime.jl/compare/v0.2.2...v0.3.0) (2026-03-24)


### Features

* add docstrings, Documenter.jl docs, and Aqua.jl tests ([1c13f31](https://github.com/JuliaGNSS/PositionVelocityTime.jl/commit/1c13f31a8eae33db951bb8355a441867fb8451bf))
