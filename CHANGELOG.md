# Changelog

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
