
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


If using single functions, please check `is_sat_healthy_and_decodable` prior to computing.

```julia
#´dc´: Decoder
julia> is_sat_healthy_and_decodable(dc)
true
```

### Position Calculation
The function `calc_PVT(dcs::Vector{GNSSDecoderState}, code_phases::Vector{Float64}, carrier_phases = zeros(length(dcs)))` provides a complete position calculation. Since this function checks the input arguments for usability, all data can be passed. Note that the input arguments needs to have the same size to prevent assignment errors. The input argument `carrier_phases` is optional due to its small effect on position calculation.

```julia
#´dcs´: Array of decoder
#´code_phases´: Array of code phases
#´carrier phases´: Array of carrier phases
julia> calc_PVT(dcs, code_phases, carrier_phases)
```

Output:
```julia
PVTSolution(
    [4.0186749839887144e6, 427051.1942215096, 4.918252576909532e6],
    -2.3409334780245904e7, 
    17.961339568218765)
```

The first member `pos` represents the user position in ECEF coordinates, the second `receiver_time_correction` the calculated travel time correction. The third member `GDOP` represents the Geometric Dilution of Precision.  



