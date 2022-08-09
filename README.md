
# PositionVelocityTime.jl
Calculates position and time by using GNSS data

## Features
* User Position calculation
* Satellite Position calculation
* Precision estimation (GDOP)

## Preparing

### Install
```julia
julia> ]
pkg> add git@github.com:JuliaGNSS/PositionVelocityTime.jl.git
```

Decoded data and code phase of satellite must be combined in the provided `SatelliteState` struct. 
```julia
using PositionVelocityTime, GNSSSignals, GNSSDecoder
decoder = GPSL1DecoderState()
gpsl1 = GPSL1()
sat_state = SatelliteState(
    decoder = decoder,
    system = gpsl1,
    code_phase = code_phase,
    carrier_phase = carrier_phase # optional
)
```
The declaration of `carrier_phase` is optional due to its small effect on the user position.

For user position computation at least 4 decoded satellites must be provided. 

## Usage

### User position calculation
The function 
```julia
# You need at least 4 satellite states
calc_PVT(sat_states)
``` 
provides a complete position calculation.