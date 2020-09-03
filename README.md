
# PVT.jl
Calculates Positions by using GPS Data

## Features
* User Position calculation
* Satellite Position calculation
* Precision estimation (GDOP)

## Usage

### Install
```julia
julia> ]add PVT
```

### Preparing
```julia
julia> using PVT, GNSSDecoder
julia> #decode Signals here
```


The decoding using the GNSSDecoder module must be completed before beginning. For user position computation at least 4 decoded satellites must be handed over. 


If using single functions, please check `can_get_sat_position` prior to computing.

```julia
#´dc´: Decoder
julia> can_get_sat_position(dc)
true
```

For automated position calculating, the function 

`calc_PVT(dcs::Vector{GNSSDecoderState}, code_phases::Vector{Float64}, carrier_phases = -1)`

is provided. Since this function checks the input arguments for usability, all data can be passed. Note that the input arguments needs to have the same size to prevent assignment errors. The input argument `carrier_phases` is optional due to its small effect on position calculation.

```julia
julia> calc_PVT(dcs, code_phases, carrier_phases)
```

Output:
```julia
([4.0188794844854493e6, 426955.64428302745, 4.918459570283906e6, -2.0419758225928288e7], 1.7019567876997732)
```

The first 3 values represent the user position in ECEF coordinates, the fourth the calculated travel time correction. The last value represents the GDOP (Geometric Dilution of Precision).  



