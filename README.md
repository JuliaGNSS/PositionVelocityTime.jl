
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

### Initializing
```julia
julia> using PositionVelocityTime, GNSSDecoder
julia> #decode Signals here
```


Decoded data and code phase of satellite must be combined in the provided `SatelliteState` struct. 
```julia
julia> sat_state = SatelliteState(
    decoder = decoder, 
    code_phase = code_phase,
    carrier_phase = carrier_phase # optional
)
```
The declaration of `carrier_phase` is optional due to its small effect on the user position.

For user position computation at least 4 decoded satellites must be provided. 

## Usage

### User position Calculation
The function 
```
calc_PVT(system, sat_states)
``` 
provides a complete position calculation.

Exemplary output:
```julia
julia> gpsl1 = GPSL1()
julia> pvt = calc_pvt(gpsl1, sat_states)
PositionVelocityTime.PVTSolution
  position: ECEF{Float64}
  time_correction: Float64 -2.392890916479146e7
  time: AstroTime.Epochs.TAIEpoch{Float64}
  dop: PositionVelocityTime.DOP
  used_sats: Array{Int64}((5,)) [2, 4, 11, 25, 30]
  sat_positions: Array{ECEF}((5,))
```